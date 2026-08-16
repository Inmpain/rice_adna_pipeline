#!/usr/bin/env python3
"""Fail-closed, receipt-based controller for the rice ecotype PCA v2 workflow.

This controller does not invent commands or infer that a scientific result
"looks right".  It executes only commands versioned in workflow.json, in the
declared order.  A stage becomes complete only after exit code 0 (command
stage) or an explicit evidence-backed PASS (manual review stage).

Every receipt is bound to the stage definition, its tracked implementation
files, the frozen config, and prerequisite receipts.  Editing any of those
makes the receipt stale and locks downstream stages until the changed stage is
run again.  Failed attempts are never overwritten and are automatically
packaged as a small debug tarball for the user to return to the code author.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import re
import subprocess
import sys
import tarfile
from pathlib import Path


HERE = Path(__file__).resolve().parent
DEFAULT_PLAN = HERE / "workflow.json"
REPO_ROOT = HERE.parents[2]
DEFAULT_CONFIG = REPO_ROOT / "scripts/ecotype_pca_v2/config/ecotype_pca_v2.yaml"
VALID_KINDS = {"command", "manual_gate", "blocked"}
VALID_RESOURCES = {"local", "login", "slurm"}


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def attempt_stamp() -> str:
    return dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def canonical_sha256(value) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return sha256_bytes(payload.encode("utf-8"))


def read_json(path: Path):
    with path.open() as handle:
        return json.load(handle)


def atomic_write_json(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with tmp.open("w") as handle:
        json.dump(value, handle, indent=2, sort_keys=True, ensure_ascii=False)
        handle.write("\n")
    os.replace(tmp, path)


def resolve_repo_file(rel_path: str) -> Path:
    path = (REPO_ROOT / rel_path).resolve()
    try:
        path.relative_to(REPO_ROOT.resolve())
    except ValueError as exc:
        raise ValueError(f"tracked path escapes repository: {rel_path}") from exc
    return path


def validate_plan(plan: dict) -> None:
    required_top = {"schema_version", "engine_contract_version", "workflow_id", "global_policy", "stages"}
    missing_top = required_top - set(plan)
    if missing_top:
        raise ValueError(f"workflow plan missing top-level fields: {sorted(missing_top)}")
    if not isinstance(plan["stages"], list) or not plan["stages"]:
        raise ValueError("workflow plan must contain at least one stage")

    seen = set()
    previous_id = None
    for index, stage in enumerate(plan["stages"]):
        for field in ("id", "title", "kind", "prerequisites", "tracked_files"):
            if field not in stage:
                raise ValueError(f"stage #{index} missing {field}")
        stage_id = stage["id"]
        if not re.fullmatch(r"[0-9]{2,3}_[a-z0-9_]+", stage_id):
            raise ValueError(f"invalid stage id: {stage_id!r}")
        if stage_id in seen:
            raise ValueError(f"duplicate stage id: {stage_id}")
        if stage["kind"] not in VALID_KINDS:
            raise ValueError(f"{stage_id}: invalid kind {stage['kind']!r}")
        if index == 0 and stage["prerequisites"]:
            raise ValueError(f"{stage_id}: first stage must have no prerequisites")
        if previous_id is not None and previous_id not in stage["prerequisites"]:
            raise ValueError(
                f"{stage_id}: must depend on immediately preceding stage {previous_id}; "
                "the workflow is deliberately linear and may not skip a gate"
            )
        unknown_prereqs = set(stage["prerequisites"]) - seen
        if unknown_prereqs:
            raise ValueError(f"{stage_id}: prerequisites are not earlier stages: {sorted(unknown_prereqs)}")

        for rel_path in stage["tracked_files"]:
            path = resolve_repo_file(rel_path)
            if not path.is_file():
                raise ValueError(f"{stage_id}: tracked file missing: {rel_path}")

        if stage["kind"] == "command":
            if not isinstance(stage.get("command"), list) or not stage["command"]:
                raise ValueError(f"{stage_id}: command stage requires a non-empty argv list")
            if stage.get("resource") not in VALID_RESOURCES:
                raise ValueError(f"{stage_id}: command stage needs resource in {sorted(VALID_RESOURCES)}")
        elif stage["kind"] == "manual_gate":
            if not stage.get("required_evidence"):
                raise ValueError(f"{stage_id}: manual gate must describe required evidence")
        elif not stage.get("blockers"):
            raise ValueError(f"{stage_id}: blocked stage must list concrete blockers")

        seen.add(stage_id)
        previous_id = stage_id


def load_plan(path: Path) -> dict:
    plan = read_json(path)
    validate_plan(plan)
    return plan


def source_revision() -> dict:
    marker = REPO_ROOT / ".source_revision"
    if marker.is_file():
        value = marker.read_text().strip()
        return {"revision": value or "UNKNOWN", "source": ".source_revision", "dirty": None}
    try:
        rev = subprocess.run(
            ["git", "-C", str(REPO_ROOT), "rev-parse", "HEAD"],
            capture_output=True, text=True, check=True,
        ).stdout.strip()
        dirty = bool(subprocess.run(
            ["git", "-C", str(REPO_ROOT), "status", "--porcelain"],
            capture_output=True, text=True, check=True,
        ).stdout.strip())
        return {"revision": rev, "source": "git", "dirty": dirty}
    except (OSError, subprocess.CalledProcessError):
        return {"revision": "UNKNOWN", "source": "unavailable", "dirty": None}


def stage_digest(plan: dict, stage: dict, config_path: Path) -> tuple[str, dict]:
    tracked = {}
    for rel_path in stage["tracked_files"]:
        tracked[rel_path] = sha256_file(resolve_repo_file(rel_path))
    ingredients = {
        "schema_version": plan["schema_version"],
        "engine_contract_version": plan["engine_contract_version"],
        "workflow_id": plan["workflow_id"],
        "global_policy": plan["global_policy"],
        "stage": stage,
        "tracked_file_sha256": tracked,
        "config_sha256": sha256_file(config_path),
    }
    return canonical_sha256(ingredients), ingredients


def receipt_path(state_dir: Path, stage_id: str) -> Path:
    return state_dir / "receipts" / f"{stage_id}.json"


def current_prerequisite_hashes(state_dir: Path, stage: dict) -> dict:
    hashes = {}
    for prereq in stage["prerequisites"]:
        path = receipt_path(state_dir, prereq)
        if not path.is_file():
            raise ValueError(f"missing prerequisite receipt: {path}")
        hashes[prereq] = sha256_file(path)
    return hashes


def validate_receipt(plan: dict, stage: dict, state_dir: Path, config_path: Path) -> tuple[bool, str]:
    path = receipt_path(state_dir, stage["id"])
    if not path.is_file():
        return False, "no receipt"
    try:
        receipt = read_json(path)
        if receipt.get("workflow_id") != plan["workflow_id"]:
            return False, "workflow ID changed"
        if receipt.get("stage_id") != stage["id"]:
            return False, "receipt stage ID mismatch"
        expected_digest, _ = stage_digest(plan, stage, config_path)
        if receipt.get("stage_digest") != expected_digest:
            return False, "stage/config/implementation digest changed"
        expected_prereqs = current_prerequisite_hashes(state_dir, stage)
        if receipt.get("prerequisite_receipt_sha256") != expected_prereqs:
            return False, "prerequisite receipt changed"
        if receipt.get("result") != "PASS":
            return False, "receipt result is not PASS"
        if stage["kind"] == "command":
            attempt_dir = Path(receipt.get("attempt_dir", "")).resolve()
            outputs = receipt.get("outputs")
            if not attempt_dir.is_dir() or not isinstance(outputs, list):
                return False, "attempt directory/output manifest missing"
            for output in outputs:
                output_path = (attempt_dir / output["path"]).resolve()
                try:
                    output_path.relative_to(attempt_dir)
                except ValueError:
                    return False, "output manifest path escapes attempt directory"
                if not output_path.is_file():
                    return False, f"recorded output missing: {output['path']}"
                if output_path.stat().st_size != output["size"]:
                    return False, f"recorded output size changed: {output['path']}"
                if sha256_file(output_path) != output["sha256"]:
                    return False, f"recorded output hash changed: {output['path']}"
        elif stage["kind"] == "manual_gate":
            evidence = receipt.get("evidence")
            if not isinstance(evidence, list) or not evidence:
                return False, "manual evidence manifest missing"
            for item in evidence:
                evidence_path = Path(item["path"])
                if not evidence_path.is_file():
                    return False, f"manual evidence missing: {evidence_path}"
                if evidence_path.stat().st_size != item["size"]:
                    return False, f"manual evidence size changed: {evidence_path}"
                if sha256_file(evidence_path) != item["sha256"]:
                    return False, f"manual evidence hash changed: {evidence_path}"
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        return False, f"invalid receipt: {exc}"
    return True, "valid receipt"


def stage_statuses(plan: dict, state_dir: Path, config_path: Path) -> list[dict]:
    statuses = []
    all_previous_complete = True
    for stage in plan["stages"]:
        valid, detail = validate_receipt(plan, stage, state_dir, config_path)
        has_receipt = receipt_path(state_dir, stage["id"]).is_file()
        if valid:
            status = "COMPLETE"
        elif not all_previous_complete:
            status = "LOCKED"
            detail = "an earlier stage is incomplete or stale"
        elif has_receipt:
            status = "STALE"
        elif stage["kind"] == "blocked":
            status = "BLOCKED"
            detail = "; ".join(stage["blockers"])
        else:
            status = "READY"
        statuses.append({"stage": stage, "status": status, "detail": detail})
        all_previous_complete = all_previous_complete and status == "COMPLETE"
    return statuses


def next_incomplete(statuses: list[dict]) -> dict | None:
    for entry in statuses:
        if entry["status"] != "COMPLETE":
            return entry
    return None


def output_manifest(attempt_dir: Path) -> list[dict]:
    rows = []
    for path in sorted(p for p in attempt_dir.rglob("*") if p.is_file()):
        stat = path.stat()
        rows.append({
            "path": str(path.relative_to(attempt_dir)),
            "size": stat.st_size,
            "sha256": sha256_file(path),
        })
    return rows


def make_debug_bundle(plan_path: Path, config_path: Path, state_dir: Path,
                      stage_id: str, attempt_dir: Path | None) -> Path:
    bundle_dir = state_dir / "debug_bundles"
    bundle_dir.mkdir(parents=True, exist_ok=True)
    bundle = bundle_dir / f"{attempt_stamp()}.{stage_id}.debug.tar.gz"
    with tarfile.open(bundle, "w:gz") as archive:
        archive.add(plan_path, arcname="workflow/workflow.json")
        archive.add(config_path, arcname="workflow/ecotype_pca_v2.yaml")
        receipts = state_dir / "receipts"
        if receipts.is_dir():
            archive.add(receipts, arcname="state/receipts")
        if attempt_dir is not None and attempt_dir.is_dir():
            archive.add(attempt_dir, arcname=f"attempts/{stage_id}/{attempt_dir.name}")
        marker = REPO_ROOT / ".source_revision"
        if marker.is_file():
            archive.add(marker, arcname="workflow/.source_revision")
    return bundle


def run_stage(args, plan: dict, stage: dict) -> int:
    if stage["kind"] != "command":
        raise ValueError(f"{stage['id']} is {stage['kind']}, not a runnable command stage")
    if stage["resource"] == "slurm" and not os.environ.get("SLURM_JOB_ID"):
        hint = stage.get("slurm_hint", "submit this controller command through sbatch")
        raise ValueError(
            f"{stage['id']} requires a SLURM allocation and will not run on a login node. {hint}"
        )

    attempt_dir = args.state_dir / "attempts" / stage["id"] / attempt_stamp()
    attempt_dir.mkdir(parents=True, exist_ok=False)
    digest, ingredients = stage_digest(plan, stage, args.config)
    prereq_hashes = current_prerequisite_hashes(args.state_dir, stage)
    revision = source_revision()

    command = [
        token.replace("{repo_root}", str(REPO_ROOT))
             .replace("{config}", str(args.config))
             .replace("{attempt_dir}", str(attempt_dir))
        for token in stage["command"]
    ]
    metadata = {
        "stage_id": stage["id"],
        "started_at": utc_now(),
        "command": command,
        "cwd": str(REPO_ROOT),
        "attempt_dir": str(attempt_dir),
        "source_revision": revision,
        "stage_digest": digest,
        "digest_ingredients": ingredients,
        "prerequisite_receipt_sha256": prereq_hashes,
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
    }
    atomic_write_json(attempt_dir / "attempt.json", metadata)

    env = os.environ.copy()
    env.update({
        "RICE_PCA_REPO_ROOT": str(REPO_ROOT),
        "RICE_PCA_CONFIG": str(args.config),
        "RICE_PCA_ATTEMPT_DIR": str(attempt_dir),
        "RICE_PCA_STATE_DIR": str(args.state_dir),
        "RICE_PCA_STAGE_ID": stage["id"],
        "RICE_PCA_STAGE_DIGEST": digest,
    })

    log_path = attempt_dir / "runner.log"
    print(f"RUN {stage['id']}: {' '.join(command)}")
    print(f"attempt: {attempt_dir}")
    with log_path.open("w") as log_handle:
        proc = subprocess.Popen(
            command, cwd=REPO_ROOT, env=env, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, text=True, bufsize=1,
        )
        assert proc.stdout is not None
        for line in proc.stdout:
            sys.stdout.write(line)
            sys.stdout.flush()
            log_handle.write(line)
            log_handle.flush()
        returncode = proc.wait()

    finished = {**metadata, "finished_at": utc_now(), "returncode": returncode}
    if returncode != 0:
        atomic_write_json(attempt_dir / "failure.json", finished)
        bundle = make_debug_bundle(args.plan, args.config, args.state_dir, stage["id"], attempt_dir)
        print(f"FAIL: {stage['id']} exited {returncode}", file=sys.stderr)
        print(f"DEBUG_BUNDLE={bundle}", file=sys.stderr)
        return returncode

    receipt = {
        "workflow_id": plan["workflow_id"],
        "stage_id": stage["id"],
        "stage_title": stage["title"],
        "kind": stage["kind"],
        "result": "PASS",
        "completed_at": utc_now(),
        "source_revision": revision,
        "stage_digest": digest,
        "prerequisite_receipt_sha256": prereq_hashes,
        "attempt_dir": str(attempt_dir),
        "command": command,
        "outputs": output_manifest(attempt_dir),
    }
    atomic_write_json(receipt_path(args.state_dir, stage["id"]), receipt)
    print(f"PASS: {stage['id']}")
    print(f"RECEIPT={receipt_path(args.state_dir, stage['id'])}")
    return 0


def accept_manual_gate(args, plan: dict, stage: dict) -> int:
    if stage["kind"] != "manual_gate":
        raise ValueError(f"{stage['id']} is not a manual review gate")
    if args.decision != "PASS":
        raise ValueError("manual gate accepts only the exact token --decision PASS")
    if not args.evidence:
        raise ValueError("manual gate requires at least one --evidence file")

    evidence = []
    for value in args.evidence:
        path = Path(value).expanduser().resolve()
        if not path.is_file():
            raise ValueError(f"evidence file not found: {path}")
        evidence.append({"path": str(path), "size": path.stat().st_size, "sha256": sha256_file(path)})

    digest, _ = stage_digest(plan, stage, args.config)
    receipt = {
        "workflow_id": plan["workflow_id"],
        "stage_id": stage["id"],
        "stage_title": stage["title"],
        "kind": stage["kind"],
        "result": "PASS",
        "completed_at": utc_now(),
        "source_revision": source_revision(),
        "stage_digest": digest,
        "prerequisite_receipt_sha256": current_prerequisite_hashes(args.state_dir, stage),
        "decision": args.decision,
        "review_note": args.note,
        "evidence": evidence,
    }
    atomic_write_json(receipt_path(args.state_dir, stage["id"]), receipt)
    print(f"PASS: manual gate {stage['id']}")
    print(f"RECEIPT={receipt_path(args.state_dir, stage['id'])}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, default=DEFAULT_PLAN)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--state-dir", type=Path, default=Path("ecotype_pca_v2_workflow_state"))
    sub = parser.add_subparsers(dest="action", required=True)
    sub.add_parser("validate-plan", help="validate workflow schema and tracked files")
    sub.add_parser("status", help="show every stage and its receipt/gate status")
    sub.add_parser("next", help="show the only stage currently allowed to advance")
    run = sub.add_parser("run", help="run the current command stage; cannot skip ahead")
    run.add_argument("stage_id", nargs="?", help="optional expected stage ID; defaults to current stage")
    accept = sub.add_parser("accept", help="accept the current manual review gate with evidence")
    accept.add_argument("stage_id")
    accept.add_argument("--decision", required=True)
    accept.add_argument("--evidence", action="append", default=[])
    accept.add_argument("--note", required=True)
    debug = sub.add_parser("debug-bundle", help="package the latest attempt and receipts")
    debug.add_argument("stage_id", nargs="?", help="defaults to the current stage")
    return parser


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    args.plan = args.plan.expanduser().resolve()
    args.config = args.config.expanduser().resolve()
    args.state_dir = args.state_dir.expanduser().resolve()

    try:
        if not args.config.is_file():
            raise ValueError(f"config not found: {args.config}")
        plan = load_plan(args.plan)
        statuses = stage_statuses(plan, args.state_dir, args.config)
        current = next_incomplete(statuses)

        if args.action == "validate-plan":
            print(f"PASS: {args.plan}")
            print(f"workflow_id={plan['workflow_id']} stages={len(plan['stages'])}")
            return 0
        if args.action == "status":
            for entry in statuses:
                print(f"{entry['stage']['id']}\t{entry['status']}\t{entry['detail']}")
            return 0
        if args.action == "next":
            if current is None:
                print("WORKFLOW_COMPLETE")
                return 0
            stage = current["stage"]
            print(f"NEXT={stage['id']} STATUS={current['status']} KIND={stage['kind']}")
            print(stage["title"])
            if stage["kind"] == "blocked":
                for blocker in stage["blockers"]:
                    print(f"BLOCKER: {blocker}")
                return 2
            return 0
        if current is None:
            raise ValueError("workflow is already complete")

        args.state_dir.mkdir(parents=True, exist_ok=True)

        requested_id = getattr(args, "stage_id", None) or current["stage"]["id"]
        if requested_id != current["stage"]["id"]:
            raise ValueError(
                f"out-of-order request refused: requested {requested_id}, "
                f"current stage is {current['stage']['id']} ({current['status']})"
            )
        stage = current["stage"]
        if current["status"] not in {"READY", "STALE"}:
            raise ValueError(f"current stage {stage['id']} is {current['status']}: {current['detail']}")

        if args.action == "run":
            return run_stage(args, plan, stage)
        if args.action == "accept":
            return accept_manual_gate(args, plan, stage)
        if args.action == "debug-bundle":
            attempts_root = args.state_dir / "attempts" / stage["id"]
            attempts = sorted(p for p in attempts_root.iterdir() if p.is_dir()) if attempts_root.is_dir() else []
            latest = attempts[-1] if attempts else None
            bundle = make_debug_bundle(args.plan, args.config, args.state_dir, stage["id"], latest)
            print(f"DEBUG_BUNDLE={bundle}")
            return 0
        raise AssertionError(args.action)
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        print(f"FATAL: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
