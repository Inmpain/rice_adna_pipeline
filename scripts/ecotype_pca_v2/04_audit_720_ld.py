#!/usr/bin/env python3
"""Phase 1 / Batch 1, script 04.

LD decay by physical-distance bin, on a handful of chromosomes (default
1,3,6,12 per docs/ECOTYPE_PCA_PANEL_QC_DESIGN.md section 3 step 5). Written
for panel B (6.7M_720) since that panel's LD structure is the least
understood of the three, but usable generically via --bfile/--keep/--label.

Audit only -- this determines what LD-pruning window is appropriate, it
does not itself prune anything. Read-only against panel inputs.

NOTE: plink2 LD-report flag syntax (--r2-unphased / --ld-window-kb /
--ld-window-r2 / --ld-window) has changed across plink2 alpha builds. This
script logs the exact command and full stdout/stderr; if the installed
build (server has v2.0.0-a.6.9LM per prior sessions) rejects a flag, paste
the error back rather than guessing a workaround -- this is a software-
compatibility question, not something to silently route around.
"""
import sys
from pathlib import Path
from collections import defaultdict

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_ecotype_v2 import base_argparser, setup_logger, load_config, run_cmd, check_output_not_present

DIST_BINS_BP = [
    (0, 1_000, "0-1kb"), (1_000, 2_000, "1-2kb"), (2_000, 5_000, "2-5kb"),
    (5_000, 10_000, "5-10kb"), (10_000, 20_000, "10-20kb"), (20_000, 50_000, "20-50kb"),
    (50_000, 100_000, "50-100kb"), (100_000, 200_000, "100-200kb"), (200_000, 500_000, "200-500kb"),
]


def bin_label(dist):
    for lo, hi, lab in DIST_BINS_BP:
        if lo <= dist < hi:
            return lab
    return None  # beyond 500kb: not part of this decay report


def main():
    ap = base_argparser(__doc__)
    ap.add_argument("--label", required=True)
    ap.add_argument("--bfile", required=True)
    ap.add_argument("--keep", default=None,
                     help="plink2 --keep sample list; for panel B this should be the "
                          "technically-valid-samples list, not a population subset")
    ap.add_argument("--chroms", default="1,3,6,12")
    ap.add_argument("--max-window-kb", type=int, default=500,
                     help="upper edge of the LD-decay report window")
    args = ap.parse_args()
    logger, _ = setup_logger(f"04_audit_720_ld.{args.label}", args.out_dir)
    load_config(args.config)

    out_dir = Path(args.out_dir)
    summary_path = out_dir / f"{args.label}.ld_decay.summary.tsv"
    check_output_not_present([summary_path], args.overwrite, logger)

    chroms = [c.strip() for c in args.chroms.split(",") if c.strip()]
    binned = defaultdict(lambda: defaultdict(list))  # chrom -> bin_label -> [r2,...]

    for chrom in chroms:
        chr_out = out_dir / f"{args.label}.ld_decay.chr{chrom}"
        cmd = ["plink2", "--bfile", args.bfile, "--chr", chrom,
               "--r2-unphased", "--ld-window-kb", str(args.max_window_kb),
               "--ld-window-r2", "0", "--ld-window", "999999",
               "--out", str(chr_out)]
        if args.keep:
            cmd[3:3] = ["--keep", args.keep]
        proc = run_cmd(cmd, logger, check=False)
        if proc.returncode != 0:
            logger.error(f"chr{chrom}: plink2 LD command failed (exit {proc.returncode}). "
                         f"See stdout/stderr above -- likely a flag-syntax mismatch for this "
                         f"plink2 build, not a data problem. Do not substitute a different "
                         f"LD tool without reporting this first.")
            continue

        vcor_path = None
        for cand in (chr_out.with_suffix(".vcor"), Path(str(chr_out) + ".vcor")):
            if cand.is_file():
                vcor_path = cand
                break
        if vcor_path is None:
            logger.error(f"chr{chrom}: expected .vcor output not found next to {chr_out} -- "
                         f"list {out_dir} and check what plink2 actually named its output")
            continue

        with open(vcor_path) as fh:
            header = fh.readline().lstrip("#").split()
            try:
                pos_a_idx = header.index("POS_A")
                pos_b_idx = header.index("POS_B")
                r2_idx = header.index("UNPHASED_R2")
            except ValueError:
                logger.error(f"chr{chrom}: unexpected .vcor header {header} -- column names "
                             f"may differ on this plink2 build, adjust before re-running")
                continue
            n_pairs = 0
            for line in fh:
                parts = line.split()
                if len(parts) <= max(pos_a_idx, pos_b_idx, r2_idx):
                    continue
                dist = abs(int(parts[pos_b_idx]) - int(parts[pos_a_idx]))
                lab = bin_label(dist)
                if lab is None:
                    continue
                binned[chrom][lab].append(float(parts[r2_idx]))
                n_pairs += 1
        logger.info(f"chr{chrom}: {n_pairs} SNP pairs binned within {args.max_window_kb}kb")

    rows = []
    for chrom in chroms:
        for lo, hi, lab in DIST_BINS_BP:
            vals = binned[chrom].get(lab, [])
            mean_r2 = sum(vals) / len(vals) if vals else None
            rows.append((chrom, lab, len(vals), mean_r2))

    with open(summary_path, "w") as fh:
        fh.write("chrom\tdistance_bin\tn_pairs\tmean_r2\n")
        for chrom, lab, n, mean_r2 in rows:
            fh.write(f"{chrom}\t{lab}\t{n}\t{'' if mean_r2 is None else f'{mean_r2:.5f}'}\n")
    logger.info(f"wrote {summary_path}")
    for chrom, lab, n, mean_r2 in rows:
        logger.info(f"chr{chrom} {lab}: n={n} mean_r2={mean_r2}")


if __name__ == "__main__":
    main()
