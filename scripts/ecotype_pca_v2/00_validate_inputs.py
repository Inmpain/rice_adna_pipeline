#!/usr/bin/env python3
"""Phase 1 / Batch 1, script 00.

Read-only. Verifies every input, tool, and label assumption the rest of
ecotype_pca_v2 depends on, before any panel manifest or audit is built.
Never modifies data. Never touches results/ecotype_pca/ (v1).

Exit codes:
  0  all SHOTGUN-track checks passed (capture-track status reported but
     does not block shotgun-track work unless --require-capture is given)
  2  an unexpected FAIL (missing file, tool, or format problem)
  3  a BLOCKED item required by --track is unresolved (e.g. capture bait BED)
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_ecotype_v2 import (base_argparser, setup_logger, load_config, tool_version,
                             read_eigenstrat_ind, write_manifest_tsv, check_output_not_present)


def check_panel_files(cfg, panel_key, logger):
    """Returns (ok: bool, found: dict of path->exists)."""
    pinfo = cfg["inputs"][panel_key]
    d = Path(pinfo["dir"])
    prefix = pinfo["prefix"]
    suffix = pinfo["filtered_suffix"]
    candidates = {
        "raw_ind": d / f"{prefix}.ind",
        "raw_snp": d / f"{prefix}.snp",
        "raw_eigenstratgeno": d / f"{prefix}.eigenstratgeno",
        "raw_geno": d / f"{prefix}.geno",
        "filtered_ind": d / f"{prefix}{suffix}.ind",
        "filtered_eigenstratgeno": d / f"{prefix}{suffix}.eigenstratgeno",
        "filtered_geno": d / f"{prefix}{suffix}.geno",
    }
    found = {k: v.is_file() for k, v in candidates.items()}
    ok_required = found["raw_ind"] and found["raw_snp"] and found["filtered_ind"]
    ok_geno = (found["raw_eigenstratgeno"] or found["raw_geno"])
    for k, exists in found.items():
        logger.info(f"[{panel_key}] {k}: {'FOUND' if exists else 'missing'} -> {candidates[k]}")
    if not ok_required:
        logger.error(f"[{panel_key}] REQUIRED files missing (raw .ind/.snp or filtered .ind)")
    if not ok_geno:
        logger.error(f"[{panel_key}] no genotype matrix found (.eigenstratgeno or .geno)")
    return (ok_required and ok_geno), found, candidates


def check_axis_labels(cfg, panel_key, candidates, logger):
    """For panel_A_3k / panel_C_civan: verify configured axis_labels exist in the
    filtered .ind with n>0 each, and (panel C only) that the total equals
    expected_axis_builder_n EXACTLY. Never auto-corrects a mismatch."""
    pcfg = cfg[panel_key]
    ind_path = candidates["filtered_ind"]
    if not ind_path.is_file():
        logger.error(f"[{panel_key}] cannot check axis labels, filtered .ind missing")
        return False
    rows = read_eigenstrat_ind(ind_path)
    counts = {}
    for r in rows:
        counts[r["label"]] = counts.get(r["label"], 0) + 1
    logger.info(f"[{panel_key}] all distinct labels in {ind_path.name}: {counts}")

    axis_labels = pcfg["axis_labels"]
    missing_labels = [lb for lb in axis_labels if counts.get(lb, 0) == 0]
    if missing_labels:
        logger.error(
            f"[{panel_key}] BLOCKED: config axis_labels not found (0 individuals) in "
            f"{ind_path.name}: {missing_labels}. Actual distinct label strings present: "
            f"{sorted(counts)}. This usually means the config's label spelling does not "
            f"match what the labeling script actually wrote -- do not guess a fix, confirm "
            f"the real string against the file and update config accordingly."
        )
        return False

    axis_n = sum(counts.get(lb, 0) for lb in axis_labels)
    expected_n = pcfg.get("expected_axis_builder_n")
    if expected_n is not None and axis_n != expected_n:
        logger.error(
            f"[{panel_key}] BLOCKED: axis-builder sample count = {axis_n}, "
            f"expected exactly {expected_n} per frozen spec. Per-label breakdown: "
            f"{ {lb: counts.get(lb, 0) for lb in axis_labels} }. STOP -- do not adjust "
            f"labels or filters to force this number; report and wait for instruction."
        )
        return False

    logger.info(f"[{panel_key}] axis-builder labels OK: n={axis_n}"
                + (f" (matches expected {expected_n})" if expected_n is not None else ""))
    return True


def check_panel_b(cfg, candidates, logger):
    ind_path = candidates["filtered_ind"]
    if not ind_path.is_file():
        ind_path = candidates["raw_ind"]
        logger.warning("[panel_B_720] no filtered .ind found, reporting on raw .ind instead")
    if not ind_path.is_file():
        logger.error("[panel_B_720] no .ind file found at all")
        return False
    rows = read_eigenstrat_ind(ind_path)
    counts = {}
    malformed = 0
    for r in rows:
        if not r["id"] or not r["label"]:
            malformed += 1
            continue
        counts[r["label"]] = counts.get(r["label"], 0) + 1
    logger.info(f"[panel_B_720] N samples={len(rows)}, malformed rows={malformed}, "
                f"per-label counts={counts}")
    logger.info("[panel_B_720] axis_mode=all_modern: no whitelist check applies -- "
                "only technically-malformed rows may be excluded later, never by "
                "population identity (spec section 5).")
    return True


def check_capture_bait_bed(cfg, logger):
    bed = cfg["inputs"].get("capture_bait_bed")
    if not bed:
        logger.error(
            "[capture] BLOCKED: inputs.capture_bait_bed is not set in config. No rice "
            "capture bait BED was found anywhere in the repo (all 4 branches searched "
            "2026-08-14) or in any server path referenced by existing docs. Per spec "
            "section 3.1/6, capture marker universe construction (CAPTURE.TV / "
            "CAPTURE.ALL for all three panels) cannot start until this file's server "
            "path is confirmed and filled into config.inputs.capture_bait_bed. "
            "SHOTGUN-track work is not blocked by this."
        )
        return False
    p = Path(bed)
    if not p.is_file():
        logger.error(f"[capture] BLOCKED: inputs.capture_bait_bed is set to {p} but that "
                      f"file does not exist on this server.")
        return False
    logger.info(f"[capture] bait BED found: {p}")
    return True


def check_tools(logger):
    tools = {
        "plink2": ["plink2", "--version"],
        "convertf": ["convertf"],
        "smartpca": ["smartpca"],
        "bedtools": ["bedtools", "--version"],
        "samtools": ["samtools", "--version"],
    }
    ok = True
    for name, cmd in tools.items():
        v = tool_version(cmd, logger)
        if v is None:
            logger.error(f"[tools] {name}: NOT FOUND on PATH")
            ok = False
        else:
            logger.info(f"[tools] {name}: {v}")
    for mod in ("yaml", "pandas", "numpy"):
        try:
            __import__(mod)
            logger.info(f"[python] {mod}: importable")
        except ImportError:
            logger.error(f"[python] {mod}: NOT importable")
            ok = False
    return ok


def main():
    ap = base_argparser(__doc__)
    ap.add_argument("--track", choices=["shotgun", "capture", "both"], default="shotgun",
                     help="which track's readiness determines the exit code")
    args = ap.parse_args()
    logger, log_path = setup_logger("00_validate_inputs", args.out_dir)
    cfg = load_config(args.config)

    report_path = Path(args.out_dir) / "00_validate_inputs.report.tsv"
    check_output_not_present([report_path], args.overwrite, logger)

    results = {}
    logger.info("=== tool/software checks ===")
    results["tools"] = check_tools(logger)

    logger.info("=== panel_A_3k file checks ===")
    ok_a, found_a, cand_a = check_panel_files(cfg, "panel_A_3k", logger)
    results["panel_A_files"] = ok_a
    if ok_a:
        results["panel_A_labels"] = check_axis_labels(cfg, "panel_A_3k", cand_a, logger)

    logger.info("=== panel_B_720 file checks ===")
    ok_b, found_b, cand_b = check_panel_files(cfg, "panel_B_720", logger)
    results["panel_B_files"] = ok_b
    if ok_b:
        results["panel_B_summary"] = check_panel_b(cfg, cand_b, logger)

    logger.info("=== panel_C_civan file checks ===")
    ok_c, found_c, cand_c = check_panel_files(cfg, "panel_C_civan", logger)
    results["panel_C_files"] = ok_c
    if ok_c:
        results["panel_C_labels"] = check_axis_labels(cfg, "panel_C_civan", cand_c, logger)

    logger.info("=== ancient BAM directory ===")
    bam_dir = Path(cfg["inputs"]["ancient_bam_dir"])
    results["ancient_bam_dir"] = bam_dir.is_dir()
    logger.info(f"[ancient] {bam_dir}: {'FOUND' if bam_dir.is_dir() else 'missing'}")
    if bam_dir.is_dir():
        bams = sorted(bam_dir.glob("*.bam"))
        logger.info(f"[ancient] {len(bams)} .bam files present: "
                    f"{[b.name for b in bams[:20]]}{' ...' if len(bams) > 20 else ''}")

    logger.info("=== capture bait BED ===")
    results["capture_bait_bed"] = check_capture_bait_bed(cfg, logger)

    logger.info("=== results_v2_root ===")
    v2_root = Path(cfg["results_v2_root"])
    results["results_v2_root_parent_writable"] = v2_root.parent.is_dir()
    logger.info(f"[layout] results_v2_root parent {v2_root.parent}: "
                f"{'exists' if v2_root.parent.is_dir() else 'MISSING'} "
                f"(script does not create it -- infrastructure setup only)")

    rows = [{"check": k, "status": "PASS" if v else "FAIL"} for k, v in results.items()]
    write_manifest_tsv(report_path, rows, ["check", "status"])
    logger.info(f"wrote {report_path}")

    shotgun_required = ["tools", "panel_A_files", "panel_A_labels", "panel_B_files",
                         "panel_B_summary", "panel_C_files", "panel_C_labels",
                         "ancient_bam_dir"]
    shotgun_ok = all(results.get(k, False) for k in shotgun_required)
    capture_ok = results.get("capture_bait_bed", False)

    logger.info(f"SHOTGUN-track readiness: {'PASS' if shotgun_ok else 'FAIL'}")
    logger.info(f"CAPTURE-track readiness: {'PASS' if capture_ok else 'BLOCKED'}")

    if args.track == "capture" and not capture_ok:
        sys.exit(3)
    if args.track == "both" and not (shotgun_ok and capture_ok):
        sys.exit(3 if shotgun_ok else 2)
    if not shotgun_ok:
        sys.exit(2)
    sys.exit(0)


if __name__ == "__main__":
    main()
