#!/usr/bin/env python3
"""Pure-Python checks for ordering, receipts, and config invalidation."""

import importlib.util
import json
import os
import tempfile
from pathlib import Path


HERE = Path(__file__).resolve().parent
WORKFLOW_DIR = HERE.parent
SCRIPT = WORKFLOW_DIR / "ecotype_pca_workflow.py"
SPEC = importlib.util.spec_from_file_location("workflow_controller", SCRIPT)
MOD = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MOD)


def check(name, condition):
    if not condition:
        raise AssertionError(name)
    print(f"PASS: {name}")


def write_passing_receipt(plan, stage, state_dir, config_path):
    digest, _ = MOD.stage_digest(plan, stage, config_path)
    attempt_dir = state_dir / f"fake_attempt_{stage['id']}"
    attempt_dir.mkdir(parents=True, exist_ok=True)
    artifact = attempt_dir / "artifact.txt"
    artifact.write_text("ok\n")
    receipt = {
        "workflow_id": plan["workflow_id"],
        "stage_id": stage["id"],
        "result": "PASS",
        "stage_digest": digest,
        "prerequisite_receipt_sha256": MOD.current_prerequisite_hashes(state_dir, stage),
        "attempt_dir": str(attempt_dir),
        "outputs": [{
            "path": "artifact.txt",
            "size": artifact.stat().st_size,
            "sha256": MOD.sha256_file(artifact),
        }],
    }
    MOD.atomic_write_json(MOD.receipt_path(state_dir, stage["id"]), receipt)
    return artifact


def write_manual_receipt(plan, stage, state_dir, config_path):
    digest, _ = MOD.stage_digest(plan, stage, config_path)
    evidence_path = state_dir / f"evidence_{stage['id']}.txt"
    evidence_path.write_text("reviewed and accepted\n")
    receipt = {
        "workflow_id": plan["workflow_id"],
        "stage_id": stage["id"],
        "result": "PASS",
        "stage_digest": digest,
        "prerequisite_receipt_sha256": MOD.current_prerequisite_hashes(state_dir, stage),
        "evidence": [{
            "path": str(evidence_path),
            "size": evidence_path.stat().st_size,
            "sha256": MOD.sha256_file(evidence_path),
        }],
    }
    MOD.atomic_write_json(MOD.receipt_path(state_dir, stage["id"]), receipt)


def main():
    plan = MOD.load_plan(WORKFLOW_DIR / "workflow.json")
    config_source = MOD.DEFAULT_CONFIG
    check("default_plan_valid", len(plan["stages"]) == 9)
    check("server_preflight_resource_is_slurm", plan["stages"][1]["resource"] == "slurm")

    saved_slurm_job_id = os.environ.pop("SLURM_JOB_ID", None)
    try:
        try:
            MOD.run_stage(None, plan, plan["stages"][1])
        except ValueError as exc:
            check("server_preflight_refuses_login_node", "requires a SLURM allocation" in str(exc))
        else:
            raise AssertionError("server_preflight_refuses_login_node")
    finally:
        if saved_slurm_job_id is not None:
            os.environ["SLURM_JOB_ID"] = saved_slurm_job_id

    with tempfile.TemporaryDirectory() as temp:
        root = Path(temp)
        state = root / "state"
        config = root / "config.yaml"
        config.write_bytes(config_source.read_bytes())

        statuses = MOD.stage_statuses(plan, state, config)
        check("first_stage_ready", statuses[0]["status"] == "READY")
        check("second_stage_locked", statuses[1]["status"] == "LOCKED")

        first = plan["stages"][0]
        artifact = write_passing_receipt(plan, first, state, config)
        statuses = MOD.stage_statuses(plan, state, config)
        check("valid_receipt_completes_first", statuses[0]["status"] == "COMPLETE")
        check("receipt_unlocks_only_next", statuses[1]["status"] == "READY")
        check("later_stage_remains_locked", statuses[2]["status"] == "LOCKED")

        artifact.write_text("tampered\n")
        statuses = MOD.stage_statuses(plan, state, config)
        check("output_change_stales_receipt", statuses[0]["status"] == "STALE")
        check("changed_output_relocks_downstream", statuses[1]["status"] == "LOCKED")
        artifact.write_text("ok\n")

        config.write_text(config.read_text() + "\n# test-only digest change\n")
        statuses = MOD.stage_statuses(plan, state, config)
        check("config_change_stales_receipt", statuses[0]["status"] == "STALE")
        check("stale_receipt_relocks_downstream", statuses[1]["status"] == "LOCKED")

    broken = json.loads(json.dumps(plan))
    broken["stages"][1]["prerequisites"] = []
    try:
        MOD.validate_plan(broken)
    except ValueError:
        print("PASS: skipped_immediate_prerequisite_rejected")
    else:
        raise AssertionError("skipped_immediate_prerequisite_rejected")

    with tempfile.TemporaryDirectory() as temp:
        root = Path(temp)
        state = root / "state"
        config = root / "config.yaml"
        config.write_bytes(config_source.read_bytes())
        for stage in plan["stages"][:3]:
            write_passing_receipt(plan, stage, state, config)
        write_manual_receipt(plan, plan["stages"][3], state, config)
        statuses = MOD.stage_statuses(plan, state, config)
        check("manual_evidence_receipt_valid", statuses[3]["status"] == "COMPLETE")
    check("implemented_stage_is_available", statuses[4]["status"] in {"READY", "COMPLETE"})
        check("post_blocker_stage_locked", statuses[5]["status"] == "LOCKED")

    print("\nALL WORKFLOW CONTROLLER TESTS PASSED")


if __name__ == "__main__":
    main()
