#!/usr/bin/env python3
"""Classify each ancient sample's fixed-marker calling result into
technical_execution (PASS/FAIL -- did 10_call_ancient_fixed_markers.py finish
and leave a readable .call_report.tsv) and scientific_projection
(formal_validation_candidate / exploratory_projection / descriptive_only,
based on callable_n), per the 2026-08-17 decisions_log.md entry: every
sample supplied is always classified and kept in the output. This script
never drops a sample for having a low or zero callable_n -- the tier label
IS the honest report, not a filter."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from fixed_projection_lib import read_tsv, refuse_existing, write_tsv


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--call-report", action="append", required=True,
                         help="repeat SAMPLE=PATH, each a .call_report.tsv from "
                              "10_call_ancient_fixed_markers.py (path may not exist, "
                              "in which case that sample is recorded as technical_execution=FAIL)")
    parser.add_argument("--formal-min", type=int, default=200,
                         help="callable_n >= this -> formal_validation_candidate")
    parser.add_argument("--exploratory-min", type=int, default=50,
                         help="callable_n >= this (and < --formal-min) -> exploratory_projection; "
                              "below this -> descriptive_only")
    parser.add_argument("--out", required=True)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def parse_entries(values: list[str]) -> list[tuple[str, Path]]:
    entries = []
    seen = set()
    for value in values:
        if "=" not in value:
            raise ValueError(f"expected SAMPLE=PATH, got {value!r}")
        sample, path_text = value.split("=", 1)
        if not sample or sample in seen:
            raise ValueError(f"empty or duplicate sample ID in {value!r}")
        seen.add(sample)
        entries.append((sample, Path(path_text)))
    if not entries:
        raise ValueError("at least one --call-report SAMPLE=PATH is required")
    return entries


def classify(callable_n: int, formal_min: int, exploratory_min: int) -> str:
    if callable_n >= formal_min:
        return "formal_validation_candidate"
    if callable_n >= exploratory_min:
        return "exploratory_projection"
    return "descriptive_only"


def main() -> int:
    args = parse_args()
    out_path = Path(args.out)
    rows = []
    try:
        refuse_existing([out_path], args.overwrite)
        entries = parse_entries(args.call_report)

        for sample, path in entries:
            technical_note = ""
            callable_n = None
            if not path.is_file():
                technical_execution = "FAIL"
                technical_note = "call_report file not found"
            else:
                try:
                    report = {r["metric"]: r["value"] for r in read_tsv(path)}
                    callable_n = int(report["called"])
                    technical_execution = "PASS"
                except (KeyError, ValueError) as exc:
                    technical_execution = "FAIL"
                    technical_note = f"unreadable call_report: {exc}"
            scientific_projection = (
                classify(callable_n, args.formal_min, args.exploratory_min)
                if callable_n is not None else "NOT_APPLICABLE_TECHNICAL_FAILURE"
            )
            rows.append({
                "sample": sample, "call_report": str(path),
                "callable_n": callable_n if callable_n is not None else "NA",
                "technical_execution": technical_execution, "technical_note": technical_note,
                "scientific_projection": scientific_projection,
            })
        write_tsv(out_path, rows,
                  ["sample", "call_report", "callable_n", "technical_execution",
                   "technical_note", "scientific_projection"])
    except (OSError, ValueError) as exc:
        print(f"FATAL: {exc}", file=sys.stderr)
        return 3

    tiers: dict[str, int] = {}
    for row in rows:
        tiers[row["scientific_projection"]] = tiers.get(row["scientific_projection"], 0) + 1
    print(f"PASS: classified {len(rows)} sample(s); " + ", ".join(f"{k}={v}" for k, v in sorted(tiers.items())))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
