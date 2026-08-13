#!/usr/bin/env python3
"""Phase 1 / Batch 1, script 08.

Panel B "paperlike_5kb" route (spec section 5, B-paperlike): from a
site-missingness/MAF-filtered (but NOT LD-pruned) SNP set, keep at most one
SNP per non-overlapping 5000bp window, deterministically, seed=20260814.

Input must be the .geno_maf_filtered.bim produced by
07_make_fixed_markers.sh --sensitivity thinning_only (i.e. geno<=0.10,
maf>=0.01 already applied, no LD pruning). This script does not re-apply
those filters itself -- it only thins.

Output is named *.paperlike_5kb.fixed.snplist -- never Wang2017_exact, per
spec's explicit naming prohibition (this is a reprocessing of a different,
higher-density matrix, not a reproduction of the published analysis).
"""
import random
import sys
from pathlib import Path
from collections import defaultdict

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_ecotype_v2 import base_argparser, setup_logger, load_config, md5_file, check_output_not_present


def main():
    ap = base_argparser(__doc__)
    ap.add_argument("--label", required=True)
    ap.add_argument("--geno-maf-bim", required=True,
                     help="output of 07_make_fixed_markers.sh --sensitivity thinning_only")
    args = ap.parse_args()
    logger, _ = setup_logger(f"08_make_5kb_thinned_markers.{args.label}", args.out_dir)
    cfg = load_config(args.config)

    window_bp = cfg["panel_B_720"]["paperlike_5kb"]["window_bp"]
    seed = cfg["panel_B_720"]["paperlike_5kb"]["seed"]
    logger.info(f"window_bp={window_bp} seed={seed} (from config.panel_B_720.paperlike_5kb)")

    out_dir = Path(args.out_dir)
    out_snplist = out_dir / f"{args.label}.paperlike_5kb.fixed.snplist"
    out_manifest = out_dir / f"{args.label}.paperlike_5kb.marker_manifest.tsv"
    check_output_not_present([out_snplist, out_manifest], args.overwrite, logger)

    by_chrom = defaultdict(list)
    n_in = 0
    with open(args.geno_maf_bim) as fh:
        for line in fh:
            parts = line.split()
            chrom, snpid, _cm, pos = parts[0], parts[1], parts[2], int(parts[3])
            by_chrom[chrom].append((pos, snpid))
            n_in += 1
    logger.info(f"input SNPs (post geno/MAF, pre-LD): {n_in}, chroms: {sorted(by_chrom)}")

    rng = random.Random(seed)
    kept = []
    for chrom in sorted(by_chrom):
        positions = sorted(by_chrom[chrom])
        windows = defaultdict(list)
        for pos, snpid in positions:
            windows[pos // window_bp].append((pos, snpid))
        for win_idx in sorted(windows):
            candidates = sorted(windows[win_idx])  # stable order before RNG draw
            chosen = candidates[rng.randrange(len(candidates))]
            kept.append((chrom, win_idx, chosen[1]))

    with open(out_snplist, "w") as fh:
        for _chrom, _win, snpid in kept:
            fh.write(snpid + "\n")
    md5 = md5_file(out_snplist)

    with open(out_manifest, "w") as fh:
        fh.write("panel\troute\tinput_snps\twindow_bp\tseed\tn_windows_with_data\t"
                  "n_snps_selected\tmd5\n")
        fh.write(f"B\tpaperlike_5kb\t{n_in}\t{window_bp}\t{seed}\t{len(kept)}\t{len(kept)}\t{md5}\n")

    logger.info(f"selected {len(kept)} SNPs (1 per populated {window_bp}bp window) "
                f"from {n_in} input SNPs")
    logger.info(f"wrote {out_snplist} (md5={md5}) and {out_manifest}")


if __name__ == "__main__":
    main()
