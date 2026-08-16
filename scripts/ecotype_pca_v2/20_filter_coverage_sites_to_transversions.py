#!/usr/bin/env python3
"""Post-filter an already-computed ancient coverage survey (from
19_survey_ancient_coverage.py) down to transversion-only sites, without
re-scanning any BAM. TV/transition status is a property of a panel site's
REF/ALT pair, not of any individual read, so this only needs the panel .snp
plus the survey's own TSV outputs.

Transversion pairs: A/C, A/T, C/G, G/T. Transition pairs (excluded): A/G, C/T.
This matches the exact definition already enforced by
10_call_ancient_fixed_markers.py's hard FATAL check for transition SNPs in a
--track TV fixed list -- that check is not touched or bypassed here; this
script exists so a track-TV marker list handed to Stage 50 never trips it.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

from fixed_projection_lib import iter_panel_snp, read_tsv, refuse_existing, sha256_file, write_tsv
from lib_ecotype_v2 import is_transversion


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel-snp", required=True)
    parser.add_argument("--union-sites", required=True, help="ancient_union_sites.tsv from 19_survey_ancient_coverage.py")
    parser.add_argument("--core-sites", required=True, help="ancient_core_sites.tsv from 19_survey_ancient_coverage.py")
    parser.add_argument("--per-sample-summary", required=True, help="per_sample_coverage_summary.tsv from 19_survey_ancient_coverage.py")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def filter_tv_rows(rows: list[dict], tv_lookup: dict) -> list[dict]:
    kept = []
    for row in rows:
        key = (row["chrom"], int(row["pos"]))
        if tv_lookup.get(key) is True:
            kept.append(row)
    return kept


def main() -> int:
    args = parse_args()
    output_dir = Path(args.out_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    union_out = output_dir / "ancient_union_sites.TV.tsv"
    core_out = output_dir / "ancient_core_sites.TV.tsv"
    per_sample_out = output_dir / "per_sample_coverage_summary.TV.tsv"
    manifest_out = output_dir / "filter_manifest.json"

    try:
        refuse_existing([union_out, core_out, per_sample_out, manifest_out], args.overwrite)

        tv_lookup = {}
        panel_site_n = 0
        panel_tv_site_n = 0
        for record in iter_panel_snp(args.panel_snp):
            panel_site_n += 1
            key = (record["chrom"], record["pos"])
            verdict = is_transversion(record["ref"], record["alt"])
            tv_lookup[key] = verdict
            if verdict is True:
                panel_tv_site_n += 1
        if panel_site_n == 0:
            raise ValueError(f"{args.panel_snp}: no usable 6-column .snp rows found")

        union_rows = read_tsv(args.union_sites)
        core_rows = read_tsv(args.core_sites)
        union_tv = filter_tv_rows(union_rows, tv_lookup)
        core_tv = filter_tv_rows(core_rows, tv_lookup)
        fieldnames = ["snp_id", "chrom", "pos", "n_samples_covered", "samples_covered"]
        write_tsv(union_out, union_tv, fieldnames)
        write_tsv(core_out, core_tv, fieldnames)

        per_sample_tv_count = Counter()
        for row in union_tv:
            for sample in row["samples_covered"].split(","):
                if sample:
                    per_sample_tv_count[sample] += 1

        per_sample_rows = read_tsv(args.per_sample_summary)
        per_sample_tv_rows = []
        for row in per_sample_rows:
            sample = row["sample"]
            covered_tv = per_sample_tv_count.get(sample, 0)
            per_sample_tv_rows.append({
                "sample": sample, "bam_path": row["bam_path"],
                "qualifying_read_n": row["qualifying_read_n"],
                "panel_sites_total": panel_tv_site_n, "panel_sites_covered": covered_tv,
                "panel_sites_covered_fraction": f"{(covered_tv / panel_tv_site_n) if panel_tv_site_n else 0.0:.10f}",
            })
        write_tsv(per_sample_out, per_sample_tv_rows,
                  ["sample", "bam_path", "qualifying_read_n", "panel_sites_total",
                   "panel_sites_covered", "panel_sites_covered_fraction"])

        manifest = {
            "schema_version": 1, "script": "20_filter_coverage_sites_to_transversions.py",
            "panel_snp": str(args.panel_snp), "panel_snp_sha256": sha256_file(args.panel_snp),
            "panel_site_n": panel_site_n, "panel_tv_site_n": panel_tv_site_n,
            "transversion_pairs": ["A/C", "A/T", "C/G", "G/T"],
            "transition_pairs_excluded": ["A/G", "C/T"],
            "source_union_sites": str(args.union_sites), "source_core_sites": str(args.core_sites),
            "source_per_sample_summary": str(args.per_sample_summary),
            "union_site_n_before": len(union_rows), "union_site_n_after_tv": len(union_tv),
            "core_site_n_before": len(core_rows), "core_site_n_after_tv": len(core_tv),
            "note": ("BAMs were not re-scanned -- TV status is a property of the panel "
                     "site's REF/ALT pair, looked up per (chrom,pos) from --panel-snp. "
                     "10_call_ancient_fixed_markers.py's hard FATAL check for transition "
                     "SNPs in a --track TV list is unaffected and still applies."),
        }
        manifest_out.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    except (OSError, ValueError, KeyError) as exc:
        print(f"FATAL: {exc}", file=sys.stderr)
        return 3

    print(f"PASS: panel TV sites={panel_tv_site_n}/{panel_site_n}; "
          f"union {len(union_rows)}->{len(union_tv)}; core {len(core_rows)}->{len(core_tv)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
