#!/usr/bin/env python3
"""Extract the snp_id column from a coverage-survey sites TSV (from
19_survey_ancient_coverage.py, or its .TV.tsv variant from
20_filter_coverage_sites_to_transversions.py) into the plain newline-
separated ID list format 09_export_fixed_reference_eigenstrat.py's
--fixed-snplist argument expects."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from fixed_projection_lib import read_tsv, refuse_existing


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sites-tsv", required=True,
                         help="ancient_union_sites.tsv / ancient_core_sites.tsv (or .TV.tsv variant)")
    parser.add_argument("--out", required=True)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    out_path = Path(args.out)
    try:
        refuse_existing([out_path], args.overwrite)
        rows = read_tsv(args.sites_tsv)
        if not rows:
            raise ValueError(f"{args.sites_tsv}: no rows (empty coverage-survey sites file)")
        ids = [row["snp_id"] for row in rows]
        if len(ids) != len(set(ids)):
            raise ValueError(f"{args.sites_tsv}: duplicate snp_id values")
        out_path.write_text("\n".join(ids) + "\n")
    except (OSError, KeyError, ValueError) as exc:
        print(f"FATAL: {exc}", file=sys.stderr)
        return 3
    print(f"PASS: wrote {len(ids)} marker IDs to {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
