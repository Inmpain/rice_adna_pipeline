#!/usr/bin/env python3
"""Phase 1 / Batch 1, script 05.

Capture marker universe step: panel SNP positions ∩ rice capture bait BED
(spec section 3.1/6 -- this intersection must happen BEFORE any capture-
track reference marker selection, and capture reference marker selection
must never proceed without it).

BLOCKED as of 2026-08-14: no capture bait BED has been found anywhere in
this repo (all 4 branches searched) or in any server path documented in
file_path.md / ECOTYPE_PCA_PANEL.md. This script refuses to run without an
explicit --bait-bed (or config.inputs.capture_bait_bed) pointing at a real
file -- it will not guess a path or fall back to running shotgun-only.

Written and ready so capture-track work can start the moment the bait BED's
server path is confirmed; do not treat its existence in this batch as
capture track being unblocked.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_ecotype_v2 import (base_argparser, setup_logger, load_config,
                             run_cmd_stream_to_file, run_cmd_discard_stdout,
                             check_output_not_present, is_transversion)


def main():
    ap = base_argparser(__doc__)
    ap.add_argument("--label", required=True)
    ap.add_argument("--bim", required=True, help="panel .bim file (from script 02)")
    ap.add_argument("--bait-bed", default=None,
                     help="rice capture bait BED path; overrides config.inputs.capture_bait_bed")
    args = ap.parse_args()
    logger, _ = setup_logger(f"05_intersect_panel_baits.{args.label}", args.out_dir)
    cfg = load_config(args.config)

    bait_bed = args.bait_bed or cfg["inputs"].get("capture_bait_bed")
    if not bait_bed:
        logger.error(
            "BLOCKED: no bait BED given (--bait-bed) and config.inputs.capture_bait_bed "
            "is not set. Per spec section 3.1, capture reference marker selection must "
            "never proceed without first intersecting against the bait BED -- refusing "
            "to run rather than silently falling back to shotgun-equivalent behavior. "
            "Confirm the bait BED's server path and either pass --bait-bed or fill in "
            "config.inputs.capture_bait_bed, then re-run."
        )
        sys.exit(3)
    bait_path = Path(bait_bed)
    if not bait_path.is_file():
        logger.error(f"BLOCKED: --bait-bed / config path does not exist on disk: {bait_path}")
        sys.exit(3)

    out_dir = Path(args.out_dir)
    snp_bed = out_dir / f"{args.label}.panel_snps.bed"
    hit_bed = out_dir / f"{args.label}.panel_snps.in_bait.bed"
    compat_snp = out_dir / f"{args.label}.capture_compatible.snp"
    check_output_not_present([snp_bed, hit_bed, compat_snp], args.overwrite, logger)

    bim_rows = []
    with open(args.bim) as fh:
        for line in fh:
            parts = line.split()
            chrom, snpid, _cm, pos, a1, a2 = parts[:6]
            bim_rows.append((chrom, int(pos), snpid, a1, a2))
    logger.info(f"panel total SNP (from {args.bim}): {len(bim_rows)}")

    bim_rows.sort(key=lambda r: (r[0], r[1]))
    with open(snp_bed, "w") as fh:
        for chrom, pos, snpid, a1, a2 in bim_rows:
            fh.write(f"{chrom}\t{pos - 1}\t{pos}\t{snpid}\n")

    # Validation-only: confirms snp_bed parses as well-formed BED. Its stdout IS
    # the (large) re-sorted file, which we don't need -- discard rather than
    # capture-then-log (capturing multi-million-line stdout into the log was
    # both a memory risk and made the log file unusable).
    run_cmd_discard_stdout(["bedtools", "sort", "-i", str(snp_bed)], logger)

    # Same reasoning: intersect's stdout is the actual result data, streamed
    # straight to hit_bed rather than captured into a Python string first.
    cmd = ["bedtools", "intersect", "-u", "-a", str(snp_bed), "-b", str(bait_path)]
    run_cmd_stream_to_file(cmd, hit_bed, logger)

    hit_ids = set()
    with open(hit_bed) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) >= 4:
                hit_ids.add(parts[3])

    by_id = {snpid: (a1, a2) for _c, _p, snpid, a1, a2 in bim_rows}
    tv_hit_ids = set()
    unresolvable_tv = 0
    for snpid in hit_ids:
        a1, a2 = by_id[snpid]
        tv = is_transversion(a1, a2)
        if tv is True:
            tv_hit_ids.add(snpid)
        elif tv is None:
            unresolvable_tv += 1

    with open(compat_snp, "w") as fh:
        for snpid in sorted(hit_ids):
            fh.write(snpid + "\n")

    logger.info(f"SNP inside bait BED (ALL): {len(hit_ids)}")
    logger.info(f"capture-compatible TV SNP: {len(tv_hit_ids)}")
    if unresolvable_tv:
        logger.warning(f"{unresolvable_tv} in-bait SNPs had non-single-base or identical "
                        f"alleles, could not classify as TV/transition -- excluded from "
                        f"the TV count, kept in the ALL count")
    logger.info(f"wrote {compat_snp} ({len(hit_ids)} SNP IDs, ALL track)")
    logger.info("(TV-restricted capture_compatible list is a subset filter applied "
                "downstream by 07_make_fixed_markers.sh, not written separately here)")


if __name__ == "__main__":
    main()
