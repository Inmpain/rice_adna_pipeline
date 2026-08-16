#!/usr/bin/env python3
"""Build the frozen ancient_callability.tsv information report."""

import argparse
import csv
import sys
from pathlib import Path

from fixed_projection_lib import POOLED_LIBRARY_TYPE, parse_sample_paths, read_calls, refuse_existing, write_tsv
from lib_ecotype_v2 import load_config


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--fixed-snp", required=True)
    parser.add_argument("--calls", action="append", required=True, help="repeat SAMPLE=CALLS_PATH")
    parser.add_argument("--metadata", help="optional TSV with sample, age, depth")
    parser.add_argument("--panel", choices=("A", "B", "C"), required=True)
    parser.add_argument("--library-type", choices=(POOLED_LIBRARY_TYPE,), required=True)
    parser.add_argument("--track", choices=("TV", "ALL"), required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        refuse_existing([args.out], args.overwrite)
        cfg = load_config(args.config)
        very_low = int(cfg["information_flags"]["very_low"])
        low = int(cfg["information_flags"]["low"])
        moderate = int(cfg["information_flags"]["moderate"])
        fixed_n = sum(1 for line in open(args.fixed_snp) if line.strip())
        if fixed_n == 0:
            raise ValueError("fixed .snp is empty")
        metadata = {}
        if args.metadata:
            with open(args.metadata, newline="") as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    metadata[row["sample"]] = row
        rows = []
        for sample, path in parse_sample_paths(args.calls):
            calls = read_calls(path)
            if len(calls) != fixed_n:
                raise ValueError(f"{sample}: {len(calls)} calls != {fixed_n} fixed markers")
            if "1" in calls:
                raise ValueError(f"{sample}: pseudo-haploid calls may not contain genotype 1")
            callable_n = sum(call in "02" for call in calls)
            if callable_n < very_low:
                flag = "VERY_LOW"
            elif callable_n < low:
                flag = "LOW"
            elif callable_n < moderate:
                flag = "MODERATE"
            else:
                flag = "HIGHER"
            item = metadata.get(sample, {})
            rows.append({
                "sample": sample, "age": item.get("age", "NA"), "depth": item.get("depth", "NA"),
                "panel": args.panel, "library_type": args.library_type, "track": args.track,
                "fixed_marker_n": fixed_n, "callable_n": callable_n,
                "callable_fraction": f"{callable_n / fixed_n:.10f}", "information_flag": flag,
            })
        write_tsv(args.out, rows, list(rows[0]))
    except (OSError, KeyError, ValueError) as exc:
        print(f"FATAL: {exc}", file=sys.stderr)
        return 3
    print(f"PASS: wrote {args.out} for {len(rows)} sample(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
