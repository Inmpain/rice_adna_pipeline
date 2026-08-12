#!/usr/bin/env python3
"""
Merge one or more ancient samples' pseudo-haploid genotype calls
(pseudo_haploid_call.py output) into a panel's own EIGENSTRAT genotype
file, as extra individual columns.

WHY THIS SHAPE, NOT SEPARATE FILES: this is the standard smartpca
-lsqproject setup (see ECOTYPE_PCA_PANEL.md section 2 / 5.2). smartpca
takes ONE genotype/snp/ind triple. Eigenvectors are computed only from
individuals whose population label (last column of .ind) appears in the
poplistname file passed to smartpca; any individual NOT in that list is
projected onto the resulting fixed axes without influencing them. So the
ancient samples must live in the SAME genotype file as the modern panel,
just given a population label that poplistname deliberately excludes --
there is no other mechanism to keep the modern eigenvectors fixed while
projecting new samples.

EFFICIENCY: merges ALL ancient samples for a given panel in ONE pass
over the panel's .eigenstratgeno file, not one pass per sample -- for
29M_3k that file has ~29.6 million lines, so re-reading it once per
sample (16 times) would be wasteful compared to once total.

Usage:
  python3 merge_ancient_into_panel.py \\
    --panel-geno panel.eigenstratgeno --panel-ind panel.ind \\
    --calls SAMPLE1=sample1.calls.txt SAMPLE2=sample2.calls.txt ... \\
    --ancient-poplabel Ancient \\
    --out-geno merged.eigenstratgeno --out-ind merged.ind

Each calls entry's file must be pseudo_haploid_call.py's output run
against the SAME panel .snp file as --panel-geno/--panel-ind, i.e. same
SNP count and order -- this script checks call-file lengths match each
other but cannot independently verify they match the panel's SNP count
beyond a line-count comparison (see the WARNING it prints if they
don't).
"""
import argparse
import sys


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--panel-geno", required=True)
    p.add_argument("--panel-ind", required=True)
    p.add_argument(
        "--calls", nargs="+", required=True,
        help="SAMPLE_ID=path/to/calls.txt, one per ancient sample (pseudo_haploid_call.py output)",
    )
    p.add_argument(
        "--ancient-poplabel", default="Ancient",
        help="population label written to .ind for every ancient sample added here -- "
             "this label must NOT be listed in the poplistname file given to smartpca, "
             "so these samples get lsqproject'd instead of used to build eigenvectors",
    )
    p.add_argument("--out-geno", required=True)
    p.add_argument("--out-ind", required=True)
    return p.parse_args()


def main():
    args = parse_args()

    samples, call_strings = [], []
    for spec in args.calls:
        if "=" not in spec:
            sys.exit(f"--calls entries must be SAMPLE_ID=path, got: {spec}")
        sample_id, path = spec.split("=", 1)
        with open(path) as f:
            s = f.read().replace("\n", "")
        samples.append(sample_id)
        call_strings.append(s)

    n_snps = len(call_strings[0])
    for sid, s in zip(samples, call_strings):
        if len(s) != n_snps:
            sys.exit(
                f"call file length mismatch: {sid} has {len(s)} calls, "
                f"expected {n_snps} (from {samples[0]}) -- these must all "
                f"be pseudo_haploid_call.py runs against the SAME panel "
                f".snp file"
            )
    sys.stderr.write(
        f"[merge] {len(samples)} ancient samples ({', '.join(samples)}), "
        f"{n_snps} SNP calls expected per sample\n"
    )

    n_lines = 0
    with open(args.panel_geno) as fin, open(args.out_geno, "w") as fout:
        for i, line in enumerate(fin):
            line = line.rstrip("\n")
            extra = "".join(s[i] for s in call_strings)
            fout.write(line + extra + "\n")
            n_lines += 1

    if n_lines != n_snps:
        sys.stderr.write(
            f"[merge] WARNING: panel .eigenstratgeno has {n_lines} lines "
            f"but call files had {n_snps} entries each -- these should "
            f"match exactly (both are supposed to follow the same panel "
            f".snp order). Output was still written using positional "
            f"row alignment, but do not trust it until this is resolved.\n"
        )
    else:
        sys.stderr.write(f"[merge] OK: {n_lines} SNP rows, counts matched\n")

    with open(args.panel_ind) as fin, open(args.out_ind, "w") as fout:
        for line in fin:
            fout.write(line if line.endswith("\n") else line + "\n")
        for sid in samples:
            fout.write(f"{sid}\tU\t{args.ancient_poplabel}\n")

    sys.stderr.write(f"[merge] wrote {args.out_geno} and {args.out_ind}\n")


if __name__ == "__main__":
    main()
