#!/usr/bin/env python3
"""Audit smartpca sample membership while allowing all-missing projections.

smartpca may omit an ancient projection row when that sample is missing at
every fixed marker.  It must never omit a modern sample that contributes to or
is projected onto the fixed modern coordinate system.
"""

import argparse

from fixed_projection_lib import read_evec, read_ind, write_tsv


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evec", required=True)
    parser.add_argument("--ind", required=True)
    parser.add_argument("--expected-n", type=int)
    parser.add_argument("--projection-label", default="Ancient")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    evec_rows = read_evec(args.evec, 10)
    ind_rows = read_ind(args.ind)
    evec_ids = [row["id"] for row in evec_rows]
    ind_by_id = {row["id"]: row for row in ind_rows}

    if len(evec_ids) != len(set(evec_ids)):
        raise SystemExit("duplicate IDs in evec")
    unknown = sorted(set(evec_ids) - set(ind_by_id))
    if unknown:
        raise SystemExit("evec contains IDs absent from ind: " + ",".join(unknown[:10]))
    if args.expected_n is not None and len(ind_rows) != args.expected_n:
        raise SystemExit(
            f"ind has {len(ind_rows)} samples, expected {args.expected_n}"
        )

    missing = sorted(set(ind_by_id) - set(evec_ids))
    missing_nonprojection = [
        sample for sample in missing
        if ind_by_id[sample]["label"] != args.projection_label
    ]
    if missing_nonprojection:
        raise SystemExit(
            "smartpca omitted non-projection sample IDs: "
            + ",".join(missing_nonprojection[:10])
        )

    write_tsv(
        args.out,
        [
            {"metric": "evec_sample_n", "value": len(evec_ids)},
            {"metric": "ind_sample_n", "value": len(ind_rows)},
            {"metric": "evec_samples_missing_from_pca", "value": len(missing)},
            {"metric": "missing_sample_ids", "value": ",".join(missing)},
            {"metric": "pc_n", "value": 10},
        ],
        ["metric", "value"],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
