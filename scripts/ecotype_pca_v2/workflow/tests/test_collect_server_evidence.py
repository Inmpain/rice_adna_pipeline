#!/usr/bin/env python3
"""Pure regression tests for parsing one-pass samtools flagstat evidence."""

import importlib.util
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "collect_server_evidence.py"
SPEC = importlib.util.spec_from_file_location("collect_server_evidence", SCRIPT)
MOD = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MOD)


def check(name, condition):
    if not condition:
        raise AssertionError(name)
    print(f"PASS: {name}")


def main():
    text = """\
100 + 2 in total (QC-passed reads + QC-failed reads)
80 + 1 paired in sequencing
60 + 1 properly paired (75.00% : 100.00%)
"""
    counts = MOD.parse_flagstat(text)
    check("flagstat_total_combines_qc_categories", counts["records"] == 102)
    check("flagstat_paired_combines_qc_categories", counts["paired_flag"] == 81)
    check("flagstat_proper_pair_combines_qc_categories", counts["proper_pair_flag"] == 61)

    try:
        MOD.parse_flagstat("10 + 0 in total (QC-passed reads + QC-failed reads)\n")
    except ValueError as exc:
        check("flagstat_missing_fields_hard_fail", "missing fields" in str(exc))
    else:
        raise AssertionError("flagstat_missing_fields_hard_fail")

    print("\nALL SERVER EVIDENCE TESTS PASSED")


if __name__ == "__main__":
    main()
