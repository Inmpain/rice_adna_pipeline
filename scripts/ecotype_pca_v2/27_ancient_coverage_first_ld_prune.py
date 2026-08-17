#!/usr/bin/env python3
"""Script 27: ancient-coverage-first LD pruning.

07_make_fixed_markers.sh's default route LD-prunes across the FULL geno/MAF-
filtered panel first, then intersects the survivors with ancient BAM
coverage. LD pruning and ancient coverage are drawn independently of each
other, so that intersection can collapse almost to nothing: Panel C ALL
track went 17,708 LD-survivors x 3,687 ancient-covered sites -> 47 in
common, because for any LD block plink2 may have kept the one representative
SNP ancient DNA never happened to sequence, while dropping the neighbor(s)
it did.

This route instead restricts LD pruning's candidate pool to the ancient-
covered sites FIRST, so every survivor is guaranteed ancient-covered by
construction -- LD pruning here only removes redundant ancient-covered SNPs
that are themselves in tight LD with each other (to avoid overweighting one
haplotype block), never a SNP that was simply the "wrong" representative of
a block ancient DNA never touched.

Input must be the .geno_maf_filtered.{bed,bim,fam} produced by
07_make_fixed_markers.sh --stage geno_maf_only (same intermediate 08's Panel
B paperlike_5kb route consumes) plus a candidate snplist that is already the
intersection of that .bim's SNP IDs with the ancient union coverage snplist
(25_intersect_snplists.py output) -- this script does not compute that
intersection itself, only the restricted LD pruning on top of it.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_ecotype_v2 import (base_argparser, setup_logger, load_config, md5_file,
                             check_output_not_present, run_cmd, resolve_marker_params,
                             tool_version)


def main():
    ap = base_argparser(__doc__)
    ap.add_argument("--panel", required=True, choices=("A", "B", "C"))
    ap.add_argument("--sensitivity", required=True,
                     choices=("primary", "S1", "S2", "S3", "S4"))
    ap.add_argument("--library-type", required=True)
    ap.add_argument("--track", required=True, choices=("TV", "ALL"))
    ap.add_argument("--geno-maf-bfile", required=True,
                     help="prefix of 07_make_fixed_markers.sh --stage geno_maf_only bed/bim/fam")
    ap.add_argument("--restrict-to", required=True,
                     help="candidate snplist: geno/MAF-passing SNP IDs already intersected "
                          "with the ancient union coverage snplist (25_intersect_snplists.py)")
    ap.add_argument("--label", required=True)
    ap.add_argument("--threads", type=int, default=2)
    args = ap.parse_args()
    logger, _ = setup_logger(f"27_ancient_coverage_first_ld_prune.{args.label}", args.out_dir)
    cfg = load_config(args.config)

    params = resolve_marker_params(cfg, args.panel, args.sensitivity)
    window_kb, r2 = params["ld_window_kb"], params["ld_r2"]
    logger.info(f"ld_window_kb={window_kb} ld_r2={r2} (panel={params['panel_key']}, "
                f"sensitivity={args.sensitivity}, route=ancient_coverage_first)")

    n_candidate = sum(1 for line in open(args.restrict_to) if line.strip())
    logger.info(f"candidate (ancient-covered, geno/MAF-passing) SNPs: {n_candidate}")
    if n_candidate == 0:
        logger.error("FATAL: --restrict-to is empty -- zero ancient-covered candidates")
        sys.exit(3)

    out_dir = Path(args.out_dir)
    ld_prefix = out_dir / f"{args.label}.{args.track}.{args.sensitivity}.ancient_first.ld"
    pruned_prefix = out_dir / f"{args.label}.{args.track}.{args.sensitivity}.ancient_first.pruned"
    out_snplist = out_dir / f"{args.label}.{args.track}.{args.sensitivity}.ancient_first.fixed.snplist"
    out_manifest = out_dir / f"{args.label}.{args.track}.{args.sensitivity}.ancient_first.marker_manifest.tsv"
    check_output_not_present([out_snplist, out_manifest], args.overwrite, logger)

    logger.info(f"plink2 version: {tool_version(['plink2', '--version'], logger)}")

    run_cmd(["plink2", "--bfile", args.geno_maf_bfile, "--extract", args.restrict_to,
             "--indep-pairwise", f"{window_kb}kb", str(r2),
             "--threads", str(args.threads), "--out", str(ld_prefix)], logger)
    run_cmd(["plink2", "--bfile", args.geno_maf_bfile, "--extract", f"{ld_prefix}.prune.in",
             "--threads", str(args.threads), "--make-bed", "--out", str(pruned_prefix)], logger)

    n_kept = 0
    with open(f"{pruned_prefix}.bim") as fh, open(out_snplist, "w") as out:
        for line in fh:
            out.write(line.split()[1] + "\n")
            n_kept += 1
    md5 = md5_file(out_snplist)

    with open(out_manifest, "w") as fh:
        fh.write("panel\tlibrary_type\ttrack\tsensitivity\troute\tcandidate_n\t"
                  "after_LD_or_thinning\tld_window_kb\tld_r2\tmd5\n")
        fh.write(f"{args.panel}\t{args.library_type}\t{args.track}\t{args.sensitivity}\t"
                 f"ancient_coverage_first\t{n_candidate}\t{n_kept}\t{window_kb}\t{r2}\t{md5}\n")

    logger.info(f"kept {n_kept}/{n_candidate} ancient-covered candidates after LD pruning "
                f"(window={window_kb}kb r2={r2})")
    logger.info(f"wrote {out_snplist} (md5={md5}) and {out_manifest}")


if __name__ == "__main__":
    main()
