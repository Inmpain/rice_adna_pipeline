#!/usr/bin/env python3
"""Write descriptive callable-count and Jaccard matrices without filtering loci."""

import argparse
import sys
from pathlib import Path

from fixed_projection_lib import parse_sample_paths, read_calls, refuse_existing


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--calls", action="append", required=True, help="repeat SAMPLE=CALLS_PATH")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    count_path = out_dir / "callable_count_matrix.tsv"
    jaccard_path = out_dir / "jaccard_matrix.tsv"
    try:
        refuse_existing([count_path, jaccard_path], args.overwrite)
        entries = parse_sample_paths(args.calls)
        calls = {sample: read_calls(path) for sample, path in entries}
        lengths = {len(value) for value in calls.values()}
        if len(lengths) != 1:
            raise ValueError("call files have different fixed-marker lengths")
        sets = {sample: {i for i, call in enumerate(value) if call in "02"} for sample, value in calls.items()}
        samples = [sample for sample, _ in entries]
        with open(count_path, "w") as count_handle, open(jaccard_path, "w") as jaccard_handle:
            header = "sample\t" + "\t".join(samples) + "\n"
            count_handle.write(header)
            jaccard_handle.write(header)
            for left in samples:
                intersections, scores = [], []
                for right in samples:
                    intersection = len(sets[left] & sets[right])
                    union = len(sets[left] | sets[right])
                    intersections.append(str(intersection))
                    scores.append(f"{(intersection / union) if union else 1.0:.10f}")
                count_handle.write(left + "\t" + "\t".join(intersections) + "\n")
                jaccard_handle.write(left + "\t" + "\t".join(scores) + "\n")
    except (OSError, ValueError) as exc:
        print(f"FATAL: {exc}", file=sys.stderr)
        return 3
    print("PASS: overlap matrices are descriptive only; no marker intersection was created")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
