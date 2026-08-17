#!/usr/bin/env python3
"""Extract one ancient sample's own covered-and-callable marker IDs from
its .call_sites.tsv (10_call_ancient_fixed_markers.py output).

Becomes that sample's private marker subset for the v1-style per-sample
private-axis projection (run_sample_panel_pca.sh's original design --
see docs/ECOTYPE_PCA_PANEL.md and plot_pca_projection.py's own docstring
on why per-sample-restricted axes are not cross-comparable with each
other or with the shared-matrix axis). FATALs if the sample has zero
covered sites -- callers must catch this and skip that sample's private
axis gracefully, not treat it as a batch-fatal error.
"""

from __future__ import annotations

import argparse
import sys

from fixed_projection_lib import read_tsv, refuse_existing

CALLABLE_STATUSES = {"CALLED_REF", "CALLED_ALT"}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--call-sites", required=True,
                         help="{sample}.{panel}.{library_type}.{track}.call_sites.tsv")
    parser.add_argument("--out", required=True)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        refuse_existing([args.out], args.overwrite)
        rows = read_tsv(args.call_sites)
        ids = [row["snp_id"] for row in rows if row["status"] in CALLABLE_STATUSES]
        if not ids:
            raise ValueError(f"{args.call_sites}: zero CALLED_REF/CALLED_ALT sites for this sample")
        with open(args.out, "w") as out:
            out.write("\n".join(ids) + "\n")
    except (OSError, KeyError, ValueError) as exc:
        print(f"FATAL: {exc}", file=sys.stderr)
        return 3
    print(f"PASS: {len(ids)} covered marker IDs written to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
