#!/usr/bin/env python3
"""Phase 1 / Batch 1, script 06.

Writes the authoritative per-panel reference/axis-builder sample list.
Everything downstream (MAF, missingness, LD, fixed marker freeze, smartpca
poplistname) must read this file, not re-derive its own notion of
"reference" from the .ind labels directly.

Panel A (29M_3k): reference = IND ∪ AUS ∪ ARO ∪ TRJ ∪ TEJ exactly.
  ADM is reported (project_labels) but never in the keep list.
Panel B (6.7M_720): reference = all samples in --ind-file except rows that
  are structurally malformed (empty id/label) or, if --smiss-file is given,
  samples with genotype-missing fraction >= --all-missing-threshold. Never
  excluded for "population looks strange" reasons.
Panel C (civan): reference = the 6 configured domesticated labels exactly.
  HARD-FAILS if the total does not equal config's expected_axis_builder_n
  (595) -- does not proceed with a "close enough" count, per spec section 5.

Read-only against the .ind file; writes only into --out-dir.
"""
import sys
from pathlib import Path
from collections import Counter

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_ecotype_v2 import (base_argparser, setup_logger, load_config,
                             read_eigenstrat_ind, write_manifest_tsv, check_output_not_present)


def build_panel_A(cfg, rows, logger):
    axis_labels = set(cfg["panel_A_3k"]["axis_labels"])
    project_labels = set(cfg["panel_A_3k"].get("project_labels", []))
    keep, other = [], Counter()
    for r in rows:
        if r["label"] in axis_labels:
            keep.append(r)
        elif r["label"] in project_labels:
            other[f"project:{r['label']}"] += 1
        else:
            other[f"unclassified:{r['label']}"] += 1
    return keep, other


def build_panel_B(cfg, rows, logger, smiss_file, threshold):
    excluded = Counter()
    smiss = {}
    if smiss_file:
        with open(smiss_file) as fh:
            header = fh.readline().split()
            iid_idx = header.index("IID")
            fmiss_idx = header.index("F_MISS")
            for line in fh:
                parts = line.split()
                smiss[parts[iid_idx]] = float(parts[fmiss_idx])
    keep = []
    for r in rows:
        if not r["id"].strip() or not r["label"].strip():
            excluded["malformed_row"] += 1
            continue
        if smiss_file and r["id"] in smiss and smiss[r["id"]] >= threshold:
            excluded["all_genotype_missing"] += 1
            continue
        keep.append(r)
    if not smiss_file:
        logger.warning("no --smiss-file given: 'all genotype missing' technical exclusions "
                        "were NOT checked for panel B. Run 03_audit_panel.py for panel B "
                        "first and pass its .smiss output here if any all-missing samples "
                        "need excluding (spec section 5: technical failures only, never "
                        "population-based exclusion).")
    return keep, excluded


def build_panel_C(cfg, rows, logger):
    axis_labels = set(cfg["panel_C_civan"]["axis_labels"])
    expected_n = cfg["panel_C_civan"]["expected_axis_builder_n"]
    keep, other = [], Counter()
    for r in rows:
        if r["label"] in axis_labels:
            keep.append(r)
        else:
            other[r["label"]] += 1
    counts = Counter(r["label"] for r in keep)
    missing = [lb for lb in axis_labels if counts.get(lb, 0) == 0]
    if missing:
        logger.error(
            f"BLOCKED: config axis_labels not present with n>0 in the .ind file: {missing}. "
            f"Distinct labels actually present: {sorted(set(r['label'] for r in rows))}. "
            f"Do not guess a corrected spelling -- confirm the real label string on the "
            f"server and update config.panel_C_civan.axis_labels only on explicit instruction."
        )
        sys.exit(3)
    if len(keep) != expected_n:
        logger.error(
            f"BLOCKED: axis-builder total = {len(keep)}, expected exactly {expected_n}. "
            f"Per-label breakdown: {dict(counts)}. STOP per spec section 17 acceptance "
            f"condition 'exactly 595 domesticated samples define axes' -- do not adjust "
            f"the label set or filters to force this number."
        )
        sys.exit(3)
    return keep, other


def main():
    ap = base_argparser(__doc__)
    ap.add_argument("--panel", required=True, choices=["A", "B", "C"])
    ap.add_argument("--label", required=True)
    ap.add_argument("--ind-file", required=True)
    ap.add_argument("--smiss-file", default=None,
                     help="panel B only: 03_audit_panel.py's .smiss output, for detecting "
                          "all-genotype-missing technical exclusions")
    ap.add_argument("--all-missing-threshold", type=float, default=0.999)
    args = ap.parse_args()
    logger, _ = setup_logger(f"06_build_reference_sample_set.{args.label}", args.out_dir)
    cfg = load_config(args.config)

    out_dir = Path(args.out_dir)
    keep_path = out_dir / f"{args.label}.reference_samples.keep"
    manifest_path = out_dir / f"{args.label}.reference_samples_manifest.tsv"
    check_output_not_present([keep_path, manifest_path], args.overwrite, logger)

    rows = read_eigenstrat_ind(args.ind_file)
    logger.info(f"panel {args.panel} ({args.label}): {len(rows)} total rows in {args.ind_file}")

    if args.panel == "A":
        keep, other = build_panel_A(cfg, rows, logger)
    elif args.panel == "B":
        keep, other = build_panel_B(cfg, rows, logger, args.smiss_file, args.all_missing_threshold)
    else:
        keep, other = build_panel_C(cfg, rows, logger)

    with open(keep_path, "w") as fh:
        for r in keep:
            fh.write(f"{r['id']}\t{r['id']}\n")  # FID=IID, plink2 --keep format

    counts = Counter(r["label"] for r in keep)
    manifest_rows = [{"label": lb, "n": n, "role": "reference"} for lb, n in sorted(counts.items())]
    manifest_rows += [{"label": k, "n": n, "role": "excluded_or_other"} for k, n in sorted(other.items())]
    manifest_rows.append({"label": "TOTAL_REFERENCE", "n": len(keep), "role": "reference"})
    write_manifest_tsv(manifest_path, manifest_rows, ["label", "n", "role"])

    logger.info(f"reference/axis-builder N = {len(keep)}, per-label: {dict(counts)}")
    logger.info(f"other (not reference): {dict(other)}")
    logger.info(f"wrote {keep_path} and {manifest_path}")


if __name__ == "__main__":
    main()
