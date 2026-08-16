#!/usr/bin/env python3
"""Export a fixed-marker EIGENSTRAT triple while retaining every modern sample."""

from __future__ import annotations

import argparse
import itertools
import json
import sys
from pathlib import Path

from fixed_projection_lib import (
    POOLED_LIBRARY_TYPE, iter_snp, read_ind, read_keep_ids, refuse_existing,
    sha256_file,
)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel", choices=("A", "B", "C"), required=True)
    parser.add_argument("--library-type", choices=(POOLED_LIBRARY_TYPE,), required=True)
    parser.add_argument("--track", choices=("TV", "ALL"), required=True)
    parser.add_argument("--panel-snp", required=True)
    parser.add_argument("--panel-geno", required=True)
    parser.add_argument("--panel-ind", required=True)
    parser.add_argument("--fixed-snplist", required=True)
    parser.add_argument("--reference-keep", required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    output_dir = Path(args.out_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = output_dir / f"{args.label}.{args.library_type}.{args.track}.fixed_reference"
    outputs = [Path(str(prefix) + suffix) for suffix in (".snp", ".eigenstratgeno", ".ind", ".poplistname", ".manifest.json")]
    try:
        refuse_existing(outputs, args.overwrite)
        ind_rows = read_ind(args.panel_ind)
        sample_index = {row["id"]: row for row in ind_rows}
        reference_ids = read_keep_ids(args.reference_keep)
        missing_reference = sorted(set(reference_ids) - set(sample_index))
        if missing_reference:
            raise ValueError(f"reference keep IDs absent from .ind: {missing_reference[:10]}")
        with open(args.fixed_snplist) as handle:
            fixed_ids = [line.strip() for line in handle if line.strip()]
        if not fixed_ids or len(fixed_ids) != len(set(fixed_ids)):
            raise ValueError("fixed SNP list is empty or contains duplicate IDs")
        fixed_set = set(fixed_ids)

        snp_out, geno_out = outputs[0], outputs[1]
        found = []
        with open(args.panel_geno) as genotype_handle, open(snp_out, "w") as snp_handle, open(geno_out, "w") as genotype_output:
            for snp, genotype_line in itertools.zip_longest(iter_snp(args.panel_snp), genotype_handle):
                if snp is None or genotype_line is None:
                    raise ValueError("panel .snp and .eigenstratgeno row counts differ")
                genotype = genotype_line.strip()
                if len(genotype) != len(ind_rows) or set(genotype) - set("0129"):
                    raise ValueError(f"invalid genotype row for SNP {snp['id']}: width/character mismatch")
                if snp["id"] in fixed_set:
                    if snp["ref"] is None:
                        raise ValueError(f"fixed SNP {snp['id']} lacks REF/ALT columns")
                    snp_handle.write(snp["line"] + "\n")
                    genotype_output.write(genotype + "\n")
                    found.append(snp["id"])
        missing_markers = sorted(fixed_set - set(found))
        if missing_markers or len(found) != len(fixed_ids):
            raise ValueError(f"fixed markers absent from panel or count mismatch: {missing_markers[:10]}")

        with open(outputs[2], "w") as handle:
            for row in ind_rows:
                handle.write(f"{row['id']}\t{row['sex']}\t{row['label']}\n")
        reference_set = set(reference_ids)
        axis_labels = []
        for row in ind_rows:
            if row["id"] in reference_set and row["label"] not in axis_labels:
                axis_labels.append(row["label"])
        with open(outputs[3], "w") as handle:
            handle.write("\n".join(axis_labels) + "\n")
        manifest = {
            "schema_version": 1, "script": "09_export_fixed_reference_eigenstrat.py",
            "panel": args.panel, "library_type": args.library_type, "track": args.track,
            "fixed_marker_n": len(found), "modern_sample_n": len(ind_rows),
            "axis_builder_n": len(reference_ids), "axis_labels": axis_labels,
            "marker_order": "source_panel_order", "ancient_in_axis": False,
            "inputs_sha256": {name: sha256_file(path) for name, path in {
                "panel_snp": args.panel_snp, "panel_geno": args.panel_geno,
                "panel_ind": args.panel_ind, "fixed_snplist": args.fixed_snplist,
                "reference_keep": args.reference_keep,
            }.items()},
        }
        outputs[4].write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    except (OSError, ValueError) as exc:
        print(f"FATAL: {exc}", file=sys.stderr)
        return 3
    print(f"PASS: exported {len(found)} fixed markers, {len(reference_ids)} axis builders, {len(ind_rows)} modern samples")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
