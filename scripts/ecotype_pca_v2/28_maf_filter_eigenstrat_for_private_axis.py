#!/usr/bin/env python3
"""Script 28: MAF-filter a raw panel's EIGENSTRAT triple for v1's per-sample
private-axis pipeline (scripts/ecotype_pca/run_sample_panel_pca.sh).

That script takes --panel-snp/--panel-geno/--panel-ind as plain file paths
and does no filtering of its own by design (each ancient sample uses
whatever of those sites it happens to cover). Left unfiltered, most sites
are low-frequency/uninformative and dilute population-structure signal --
this produces a MAF-filtered EIGENSTRAT triple (same row order, fewer rows)
so run_sample_panel_pca.sh can be pointed at *that* instead of the raw
panel, with zero changes to that script itself.

MAF is computed only over the reference/axis-building samples (the same
--keep list 06_build_reference_sample_set.py already produces for the
shared-track runners), not the whole panel, so wild/outgroup samples in the
panel don't skew which sites look common. Output keeps every sample column
-- only site selection is restricted by the reference subset, not who's in
the final file.

.snp and .eigenstratgeno are read in strict lockstep (line N of one
corresponds to line N of the other) -- unlike fixed_projection_lib.iter_panel_snp,
this does NOT silently skip malformed .snp lines, because doing so here
would desync the two files' row alignment.
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_ecotype_v2 import load_config


def read_ind_ids(path):
    ids = []
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            if not parts:
                continue
            ids.append(parts[0])
    return ids


def read_keep_ids(path):
    ids = set()
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 2:
                continue
            ids.add(parts[1])  # FID\tIID -- IID matches .ind column 1
    return ids


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--panel", required=True, choices=("A", "B", "C"))
    ap.add_argument("--panel-snp", required=True)
    ap.add_argument("--panel-geno", required=True)
    ap.add_argument("--panel-ind", required=True)
    ap.add_argument("--keep", required=True,
                     help="reference_samples.keep (FID\\tIID per line, from "
                          "06_build_reference_sample_set.py) -- restricts which "
                          "columns count toward each site's MAF, not which "
                          "columns are kept in the output")
    ap.add_argument("--label", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()

    panel_key = {"A": "panel_A_3k", "B": "panel_B_720", "C": "panel_C_civan"}[args.panel]
    cfg = load_config(args.config)
    maf_threshold = float(cfg[panel_key]["maf"])

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_snp = out_dir / f"{args.label}.maf{maf_threshold}.snp"
    out_geno = out_dir / f"{args.label}.maf{maf_threshold}.eigenstratgeno"
    out_manifest = out_dir / f"{args.label}.maf{maf_threshold}.manifest.tsv"
    if not args.overwrite:
        for p in (out_snp, out_geno, out_manifest):
            if p.exists():
                sys.exit(f"FATAL: refusing to overwrite existing output (pass --overwrite): {p}")

    sample_ids = read_ind_ids(args.panel_ind)
    keep_ids = read_keep_ids(args.keep)
    mask = [sid in keep_ids for sid in sample_ids]
    n_ref = sum(mask)
    if n_ref == 0:
        sys.exit(f"FATAL: zero of {len(sample_ids)} panel samples matched --keep {args.keep}")
    print(f"MAF computed over {n_ref}/{len(sample_ids)} reference samples "
          f"(panel={panel_key}, maf>={maf_threshold})", file=sys.stderr)

    n_in = n_kept = 0
    with open(args.panel_snp) as snp_fh, open(args.panel_geno) as geno_fh, \
         open(out_snp, "w") as snp_out, open(out_geno, "w") as geno_out:
        for snp_line, geno_line in zip(snp_fh, geno_fh):
            n_in += 1
            geno_row = geno_line.rstrip("\n")
            if len(geno_row) != len(sample_ids):
                raise ValueError(f".eigenstratgeno row {n_in}: {len(geno_row)} genotypes "
                                  f"!= {len(sample_ids)} samples in {args.panel_ind}")
            alt_alleles = 0
            called_alleles = 0
            for ref_here, code in zip(mask, geno_row):
                if not ref_here or code == "9":
                    continue
                alt_alleles += int(code)
                called_alleles += 2
            if called_alleles == 0:
                continue
            alt_freq = alt_alleles / called_alleles
            maf = min(alt_freq, 1.0 - alt_freq)
            if maf < maf_threshold:
                continue
            snp_out.write(snp_line if snp_line.endswith("\n") else snp_line + "\n")
            geno_out.write(geno_line if geno_line.endswith("\n") else geno_line + "\n")
            n_kept += 1

    snp_lines = sum(1 for _ in open(args.panel_snp))
    geno_lines = sum(1 for _ in open(args.panel_geno))
    if snp_lines != geno_lines:
        raise ValueError(f".snp has {snp_lines} lines but .eigenstratgeno has {geno_lines} -- "
                          f"not row-aligned, cannot trust this filter's output")

    with open(out_manifest, "w") as fh:
        fh.write("panel\tmaf_threshold\treference_samples_n\ttotal_samples_n\t"
                  "sites_in\tsites_kept\n")
        fh.write(f"{panel_key}\t{maf_threshold}\t{n_ref}\t{len(sample_ids)}\t{n_in}\t{n_kept}\n")

    print(f"PASS: {n_kept}/{n_in} sites pass MAF>={maf_threshold} (over {n_ref} reference "
          f"samples). Wrote {out_snp} and {out_geno}", file=sys.stderr)
    print(f"Use --panel-ind {args.panel_ind} unchanged (samples don't change, only sites do)",
          file=sys.stderr)


if __name__ == "__main__":
    main()
