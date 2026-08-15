#!/usr/bin/env python3
"""Phase 1 / Batch 1, script 01.

Panel-level manifest: N samples (raw vs filtered), N SNPs, per-label sample
counts, chromosome list, biallelic/REF-ALT column presence. This is the
lightweight "what do we actually have" snapshot -- NOT the fixed-marker
manifest (that is 07_make_fixed_markers.sh's job, after MAF/LD/TV-ALL).

Read-only against panel inputs. Writes only to --out-dir.
"""
import sys
from pathlib import Path
from collections import Counter

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_ecotype_v2 import (base_argparser, setup_logger, load_config,
                             read_eigenstrat_ind, iter_eigenstrat_snp,
                             read_eigenstrat_snp_header_probe, write_manifest_tsv,
                             check_output_not_present, find_duplicate_ids)

PANEL_KEYS = ["panel_A_3k", "panel_B_720", "panel_C_civan"]


def manifest_for_panel(cfg, panel_key, logger):
    pinfo = cfg["inputs"][panel_key]
    d = Path(pinfo["dir"])
    prefix = pinfo["prefix"]
    suffix = pinfo["filtered_suffix"]

    raw_ind = d / f"{prefix}.ind"
    raw_snp = d / f"{prefix}.snp"
    filtered_ind = d / f"{prefix}{suffix}.ind"

    if not raw_ind.is_file() or not raw_snp.is_file():
        logger.error(f"[{panel_key}] missing raw .ind/.snp, cannot build manifest -- "
                      f"run 00_validate_inputs.py first")
        return None, False

    ok = True
    raw_rows = read_eigenstrat_ind(raw_ind)
    raw_n = len(raw_rows)
    raw_labels = Counter(r["label"] for r in raw_rows)

    dup_sample_ids = find_duplicate_ids(r["id"] for r in raw_rows)
    if dup_sample_ids:
        sample = list(dup_sample_ids.items())[:10]
        logger.error(f"[{panel_key}] ERROR: {len(dup_sample_ids)} duplicate sample IDs "
                      f"in {raw_ind}, e.g. {sample}")
        ok = False

    filtered_n = None
    filtered_labels = {}
    if filtered_ind.is_file():
        filtered_rows = read_eigenstrat_ind(filtered_ind)
        filtered_n = len(filtered_rows)
        filtered_labels = Counter(r["label"] for r in filtered_rows)
    else:
        logger.warning(f"[{panel_key}] no filtered .ind found at {filtered_ind}")

    col_counts = read_eigenstrat_snp_header_probe(raw_snp, n=10)
    ncol = col_counts[0] if col_counts else None
    has_ref_alt = ncol == 6
    if len(set(col_counts)) > 1:
        logger.error(f"[{panel_key}] ERROR: inconsistent .snp column count across first "
                      f"{len(col_counts)} lines: {col_counts} -- format problem, "
                      f"do not proceed with this panel until resolved")
        ok = False

    chroms = Counter()
    n_snps = 0
    non_biallelic_alpha = 0
    snp_ids = []
    try:
        for rec in iter_eigenstrat_snp(raw_snp):
            n_snps += 1
            chroms[rec["chrom"]] += 1
            snp_ids.append(rec["snpid"])
            if rec["ref"] is not None and rec["alt"] is not None:
                if len(rec["ref"]) != 1 or len(rec["alt"]) != 1:
                    non_biallelic_alpha += 1
    except ValueError as e:
        logger.error(f"[{panel_key}] ERROR: {raw_snp} failed to parse: {e}")
        ok = False

    dup_snp_ids = find_duplicate_ids(snp_ids)
    if dup_snp_ids:
        sample = list(dup_snp_ids.items())[:10]
        logger.error(f"[{panel_key}] ERROR: {len(dup_snp_ids)} duplicate SNP IDs "
                      f"in {raw_snp}, e.g. {sample}")
        ok = False

    logger.info(f"[{panel_key}] raw N samples={raw_n}, raw N snps={n_snps}, "
                f"has_ref_alt_columns={has_ref_alt}, chroms={sorted(chroms, key=str)}")
    logger.info(f"[{panel_key}] raw per-label counts: {dict(raw_labels)}")
    if filtered_n is not None:
        logger.info(f"[{panel_key}] filtered N samples={filtered_n}, "
                    f"filtered per-label counts: {dict(filtered_labels)}")

    row = {
        "panel": panel_key,
        "raw_n_samples": raw_n,
        "filtered_n_samples": filtered_n if filtered_n is not None else "NA",
        "n_snps": n_snps,
        "n_chroms": len(chroms),
        "chrom_list": ";".join(sorted(chroms, key=lambda x: (len(x), x))),
        "snp_columns": ncol,
        "has_ref_alt_columns": has_ref_alt,
        "non_biallelic_alpha_alleles": non_biallelic_alpha,
        "raw_label_counts": ";".join(f"{k}={v}" for k, v in sorted(raw_labels.items())),
        "filtered_label_counts": ";".join(f"{k}={v}" for k, v in sorted(filtered_labels.items())),
    }
    return row, ok


def main():
    ap = base_argparser(__doc__)
    args = ap.parse_args()
    logger, _ = setup_logger("01_make_panel_manifest", args.out_dir)
    cfg = load_config(args.config)

    out_path = Path(args.out_dir) / "panel_manifest.tsv"
    check_output_not_present([out_path], args.overwrite, logger)

    rows = []
    fail = False
    for panel_key in PANEL_KEYS:
        row, ok = manifest_for_panel(cfg, panel_key, logger)
        if not ok:
            fail = True
        if row is not None:
            rows.append(row)

    fieldnames = ["panel", "raw_n_samples", "filtered_n_samples", "n_snps", "n_chroms",
                  "chrom_list", "snp_columns", "has_ref_alt_columns",
                  "non_biallelic_alpha_alleles", "raw_label_counts", "filtered_label_counts"]
    write_manifest_tsv(out_path, rows, fieldnames)
    logger.info(f"wrote {out_path}")

    sys.exit(1 if fail else 0)


if __name__ == "__main__":
    main()
