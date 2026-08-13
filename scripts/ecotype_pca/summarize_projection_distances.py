#!/usr/bin/env python3
"""
Rank one projected sample's nearest modern population centroids in PC
space, from a smartpca .evec file -- docs/ECOTYPE_PCA_EXECUTION_PLAN.md
section 5's required report format ("最近现代群体: TEJ 76% / 次近: ...",
not a bare "falls near TEJ" statement).

WHY min-pop-size EXISTS: a population label with very few individuals
(e.g. this project's Civáň wild-outgroup singletons -- O. barthii,
O. glaberrima etc, n=1 each, see build_civan_population_labels.py) has a
"centroid" that is really just one individual's own coordinates, not a
population estimate -- computing a meaningful nearest-population ranking
requires excluding these from the distance ranking, not silently letting
a n=1 group win just because a single real ancient sample's noisy
projection happened to land near that one point. They stay computable
via a lower --min-pop-size if specifically wanted, default 5.

Usage:
  python3 summarize_projection_distances.py \\
    --evec sample.panel.TV.evec \\
    --sample SAMPLE_ID \\
    --out sample.panel.TV.nearest.tsv
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from collections import defaultdict
from pathlib import Path


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--evec", required=True)
    parser.add_argument("--sample", required=True, help="individual ID (first column of .evec) to rank")
    parser.add_argument("--num-pcs", type=int, default=2, help="how many leading PCs to use for distance (default 2)")
    parser.add_argument(
        "--exclude-label",
        action="append",
        default=None,
        help="population label(s) to exclude from ranking (repeatable); defaults to just 'Ancient'",
    )
    parser.add_argument("--min-pop-size", type=int, default=5, help="ignore centroids built from fewer individuals than this")
    parser.add_argument("--out", default=None, help="optional TSV of the full ranking")
    return parser.parse_args(argv)


def validate_args(args: argparse.Namespace) -> None:
    if not Path(args.evec).is_file():
        raise FileNotFoundError(f".evec file not found: {args.evec}")
    if args.num_pcs < 1:
        raise ValueError("--num-pcs must be >= 1")
    if args.min_pop_size < 1:
        raise ValueError("--min-pop-size must be >= 1")


def load_evec(path: Path, num_pcs: int) -> tuple[dict[str, tuple[float, ...]], dict[str, list[tuple[float, ...]]]]:
    """Return (sample_id -> coords, label -> [coords, ...]) -- last field of each row is the label."""
    by_sample: dict[str, tuple[float, ...]] = {}
    by_label: dict[str, list[tuple[float, ...]]] = defaultdict(list)
    with path.open() as handle:
        first = handle.readline()
        if not first.startswith("#eigvals") and not first.strip().startswith("#"):
            raise ValueError(f"{path}: expected first line to start with '#eigvals', got {first[:40]!r}")
        for line_no, line in enumerate(handle, start=2):
            fields = line.split()
            if len(fields) < 2 + num_pcs:
                raise ValueError(f"{path}:{line_no}: expected at least {2 + num_pcs} fields, got {len(fields)}: {line!r}")
            sample_id = fields[0]
            label = fields[-1]
            coords = tuple(float(x) for x in fields[1 : 1 + num_pcs])
            if sample_id in by_sample:
                raise ValueError(f"{path}:{line_no}: duplicate individual ID {sample_id!r}")
            by_sample[sample_id] = coords
            by_label[label].append(coords)
    return by_sample, by_label


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    exclude_labels = set(args.exclude_label) if args.exclude_label else {"Ancient"}
    try:
        validate_args(args)
        by_sample, by_label = load_evec(Path(args.evec), args.num_pcs)

        if args.sample not in by_sample:
            raise ValueError(f"sample {args.sample!r} not found in {args.evec} (found {len(by_sample)} individuals)")
        target = by_sample[args.sample]

        ranked = []
        for label, pts in by_label.items():
            if label in exclude_labels:
                continue
            n = len(pts)
            if n < args.min_pop_size:
                continue
            centroid = tuple(sum(p[i] for p in pts) / n for i in range(args.num_pcs))
            dist = math.sqrt(sum((target[i] - centroid[i]) ** 2 for i in range(args.num_pcs)))
            ranked.append((label, n, dist, centroid))
        ranked.sort(key=lambda row: row[2])

        if not ranked:
            raise ValueError(
                f"no population label in {args.evec} has >= {args.min_pop_size} individuals "
                "(after excluding {sorted(exclude_labels)}) -- nothing to rank against"
            )

    except (OSError, ValueError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    coord_str = ", ".join(f"{c:.4f}" for c in target)
    print(f"[nearest] {args.sample}: PC1..PC{args.num_pcs} = ({coord_str})", file=sys.stderr)
    for rank, (label, n, dist, centroid) in enumerate(ranked, start=1):
        centroid_str = ", ".join(f"{c:.4f}" for c in centroid)
        print(f"[nearest]   #{rank} {label:22s} n={n:4d} dist={dist:.4f} centroid=({centroid_str})", file=sys.stderr)

    if args.out:
        out_path = Path(args.out)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with out_path.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(
                ["sample", "rank", "population", "n", "distance"] + [f"pc{i + 1}" for i in range(args.num_pcs)]
            )
            for rank, (label, n, dist, centroid) in enumerate(ranked, start=1):
                writer.writerow(
                    [args.sample, rank, label, n, f"{dist:.6f}"] + [f"{c:.6f}" for c in centroid]
                )
        print(f"[done] wrote {args.out}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
