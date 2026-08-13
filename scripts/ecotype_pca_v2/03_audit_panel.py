#!/usr/bin/env python3
"""Phase 1 / Batch 1, script 03.

Generic per-panel MAF / missingness / spacing audit. Does not filter or
delete anything -- audit only, per spec section 5's "Panel B: audit first,
do not drop data or run PCA yet", generalized to all three panels since
Phase 1 needs the same MAF/missingness picture for A and C too (computed
on the reference/axis-builder sample set only for A and C; on all
technically-valid samples for B, per --keep or its absence).

Produces, prefixed by --label (e.g. "720", "29m3k", "civan"):
  {label}.audit.samples.tsv
  {label}.audit.maf.tsv
  {label}.audit.missingness.tsv
  {label}.audit.spacing.tsv
  {label}.audit.summary.txt

Requires PLINK2 bed/bim/fam (output of 02_convert_eigenstrat_for_plink.sh).
"""
import sys
from pathlib import Path
from collections import Counter

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_ecotype_v2 import (base_argparser, setup_logger, load_config, run_cmd,
                             read_eigenstrat_ind, check_output_not_present)

MAF_BINS = [(0, 0), (0, 0.001), (0.001, 0.005), (0.005, 0.01),
            (0.01, 0.05), (0.05, 0.10), (0.10, float("inf"))]
MAF_BIN_LABELS = ["0", "0-0.001", "0.001-0.005", "0.005-0.01",
                   "0.01-0.05", "0.05-0.10", ">0.10"]

MISS_BINS = [(0, 0.01), (0.01, 0.05), (0.05, 0.10), (0.10, 0.20), (0.20, float("inf"))]
MISS_BIN_LABELS = ["0-0.01", "0.01-0.05", "0.05-0.10", "0.10-0.20", ">0.20"]


def bin_value(v, bins, labels):
    if v == 0:
        return labels[0] if bins[0] == (0, 0) else _bin_nonzero(v, bins, labels)
    return _bin_nonzero(v, bins, labels)


def _bin_nonzero(v, bins, labels):
    for (lo, hi), lab in zip(bins, labels):
        if lo < v <= hi or (lo == 0 and v == 0):
            return lab
    return labels[-1]


def percentile(sorted_vals, p):
    if not sorted_vals:
        return None
    k = (len(sorted_vals) - 1) * p
    f = int(k)
    c = min(f + 1, len(sorted_vals) - 1)
    if f == c:
        return sorted_vals[f]
    return sorted_vals[f] + (sorted_vals[c] - sorted_vals[f]) * (k - f)


def main():
    ap = base_argparser(__doc__)
    ap.add_argument("--label", required=True, help="output filename prefix, e.g. 720, 29m3k, civan")
    ap.add_argument("--bfile", required=True, help="PLINK bfile prefix (bed/bim/fam), from script 02")
    ap.add_argument("--ind-file", required=True,
                     help="EIGENSTRAT .ind file to source per-sample population labels from "
                          "(for per-pop sample counts; not used for genotype computation)")
    ap.add_argument("--keep", default=None,
                     help="plink2 --keep sample list (FID IID); omit to use every sample in --bfile "
                          "(panel B / axis_mode=all_modern). For panel A/C this MUST be the "
                          "reference/axis-builder list from 06_build_reference_sample_set.py -- "
                          "MAF/missingness/LD must never be computed on ancient or non-reference "
                          "modern individuals.")
    args = ap.parse_args()
    logger, _ = setup_logger(f"03_audit_panel.{args.label}", args.out_dir)
    load_config(args.config)  # validated for side-effect (version check); no params used here

    out_dir = Path(args.out_dir)
    out_paths = [out_dir / f"{args.label}.audit.{suf}"
                 for suf in ("samples.tsv", "maf.tsv", "missingness.tsv", "spacing.tsv", "summary.txt")]
    check_output_not_present(out_paths, args.overwrite, logger)

    freq_out = out_dir / f"{args.label}.audit.freq"
    miss_out = out_dir / f"{args.label}.audit.miss"

    cmd_freq = ["plink2", "--bfile", args.bfile, "--freq", "--out", str(freq_out)]
    cmd_miss = ["plink2", "--bfile", args.bfile, "--missing", "--out", str(miss_out)]
    if args.keep:
        cmd_freq[3:3] = ["--keep", args.keep]
        cmd_miss[3:3] = ["--keep", args.keep]
        logger.info(f"restricting MAF/missingness computation to --keep {args.keep} "
                    f"(reference/axis-builder samples only)")
    else:
        logger.warning("no --keep given: MAF/missingness computed on ALL samples in --bfile. "
                        "This is only correct for panel B (axis_mode=all_modern); for panel "
                        "A/C this must never be used for the frozen reference marker set.")

    run_cmd(cmd_freq, logger)
    run_cmd(cmd_miss, logger)

    afreq_path = Path(str(freq_out) + ".afreq")
    vmiss_path = Path(str(miss_out) + ".vmiss")
    smiss_path = Path(str(miss_out) + ".smiss")
    for p in (afreq_path, vmiss_path, smiss_path):
        if not p.is_file():
            logger.error(f"expected plink2 output missing: {p} (plink2 output filenames can "
                          f"vary by version -- check {out_dir} and adjust if needed)")
            sys.exit(2)

    # --- MAF bins ---
    maf_counts = Counter()
    n_sites = 0
    with open(afreq_path) as fh:
        header = fh.readline().split()
        alt_freq_idx = header.index("ALT_FREQS")
        for line in fh:
            parts = line.split()
            freq = float(parts[alt_freq_idx])
            maf = min(freq, 1 - freq)
            maf_counts[bin_value(maf, MAF_BINS, MAF_BIN_LABELS)] += 1
            n_sites += 1
    with open(out_dir / f"{args.label}.audit.maf.tsv", "w") as fh:
        fh.write("maf_bin\tn_sites\tfraction\n")
        for lab in MAF_BIN_LABELS:
            n = maf_counts.get(lab, 0)
            fh.write(f"{lab}\t{n}\t{n / n_sites if n_sites else 0:.6f}\n")
    logger.info(f"MAF bins: {dict(maf_counts)} (n_sites={n_sites})")

    # --- site missingness bins ---
    site_miss_counts = Counter()
    with open(vmiss_path) as fh:
        header = fh.readline().split()
        fmiss_idx = header.index("F_MISS")
        for line in fh:
            parts = line.split()
            fmiss = float(parts[fmiss_idx])
            site_miss_counts[bin_value(fmiss, MISS_BINS, MISS_BIN_LABELS)] += 1

    # --- sample missingness ---
    sample_miss = []
    with open(smiss_path) as fh:
        header = fh.readline().split()
        fmiss_idx = header.index("F_MISS")
        iid_idx = header.index("IID")
        for line in fh:
            parts = line.split()
            sample_miss.append((parts[iid_idx], float(parts[fmiss_idx])))

    with open(out_dir / f"{args.label}.audit.missingness.tsv", "w") as fh:
        fh.write("metric\tbin_or_id\tvalue\n")
        for lab in MISS_BIN_LABELS:
            fh.write(f"site_missingness_bin\t{lab}\t{site_miss_counts.get(lab, 0)}\n")
        for iid, fmiss in sample_miss:
            fh.write(f"sample_missingness\t{iid}\t{fmiss:.6f}\n")
    logger.info(f"site missingness bins: {dict(site_miss_counts)}")
    logger.info(f"sample missingness range: min={min(v for _, v in sample_miss):.4f} "
                f"max={max(v for _, v in sample_miss):.4f}" if sample_miss else "no samples")

    # --- per-pop sample counts (from .ind, not from plink; --keep may subset this) ---
    ind_rows = read_eigenstrat_ind(args.ind_file)
    keep_ids = None
    if args.keep:
        keep_ids = set()
        with open(args.keep) as fh:
            for line in fh:
                parts = line.split()
                if parts:
                    keep_ids.add(parts[-1])  # IID is last column in a 2-col FID/IID keep file
    pop_counts = Counter()
    for r in ind_rows:
        if keep_ids is not None and r["id"] not in keep_ids:
            continue
        pop_counts[r["label"]] += 1
    with open(out_dir / f"{args.label}.audit.samples.tsv", "w") as fh:
        fh.write("label\tn\n")
        for lab, n in sorted(pop_counts.items()):
            fh.write(f"{lab}\t{n}\n")
    logger.info(f"per-pop sample counts (post --keep if given): {dict(pop_counts)}")

    # --- adjacent SNP spacing, pooled within-chromosome diffs ---
    bim_path = Path(args.bfile + ".bim")
    by_chrom = {}
    with open(bim_path) as fh:
        for line in fh:
            parts = line.split()
            chrom, pos = parts[0], int(parts[3])
            by_chrom.setdefault(chrom, []).append(pos)
    diffs = []
    for chrom, positions in by_chrom.items():
        positions.sort()
        diffs.extend(b - a for a, b in zip(positions, positions[1:]))
    diffs.sort()
    n_diffs = len(diffs)
    stats = {
        "median": percentile(diffs, 0.5),
        "P10": percentile(diffs, 0.10),
        "P25": percentile(diffs, 0.25),
        "P75": percentile(diffs, 0.75),
        "P90": percentile(diffs, 0.90),
        "fraction_lt_1kb": sum(1 for d in diffs if d < 1000) / n_diffs if n_diffs else None,
        "fraction_lt_5kb": sum(1 for d in diffs if d < 5000) / n_diffs if n_diffs else None,
        "fraction_lt_10kb": sum(1 for d in diffs if d < 10000) / n_diffs if n_diffs else None,
    }
    with open(out_dir / f"{args.label}.audit.spacing.tsv", "w") as fh:
        fh.write("metric\tvalue\n")
        for k, v in stats.items():
            fh.write(f"{k}\t{v}\n")
    logger.info(f"adjacent SNP spacing (bp, within-chromosome, n_diffs={n_diffs}): {stats}")

    with open(out_dir / f"{args.label}.audit.summary.txt", "w") as fh:
        fh.write(f"panel label: {args.label}\n")
        fh.write(f"bfile: {args.bfile}\n")
        fh.write(f"keep list: {args.keep or '(none -- all samples in bfile)'}\n")
        fh.write(f"n_sites: {n_sites}\n")
        fh.write(f"n_samples_in_pop_table: {sum(pop_counts.values())}\n")
        fh.write(f"maf_bins: {dict(maf_counts)}\n")
        fh.write(f"site_missingness_bins: {dict(site_miss_counts)}\n")
        fh.write(f"spacing_stats: {stats}\n")
    logger.info(f"wrote audit outputs to {out_dir} with prefix {args.label}.audit.*")


if __name__ == "__main__":
    main()
