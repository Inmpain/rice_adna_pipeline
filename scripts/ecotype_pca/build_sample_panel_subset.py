#!/usr/bin/env python3
"""
Row-subset a panel's .snp/.eigenstratgeno down to only the SNPs one
ancient sample's pseudo_haploid_call.py output actually covers, dropping
everything else -- docs/ECOTYPE_PCA_EXECUTION_PLAN.md section 5.

WHY THIS SCRIPT EXISTS SEPARATELY FROM pseudo_haploid_call.py:
pseudo_haploid_call.py deliberately keeps its output the SAME length and
row order as the panel's .snp file, uncovered sites written as "9" --
that invariant is required for merge_ancient_into_panel.py to append it
as a same-length extra column. But for the PRIMARY per-sample analysis
(section 5's design, as opposed to the auxiliary fixed-404K-CoreSNP
projection which deliberately keeps the full panel), running smartpca on
a panel where one ancient column is 99%+ "9" wastes computation and, per
section 5.2's earlier discussion in ECOTYPE_PCA_PANEL.md, is not what
the sample-specific-subset design calls for -- it calls for actually
shrinking the SNP set down to what this one sample can speak to, for
both the modern panel and the ancient column together. This script does
exactly that shrink, in the same "stream .snp and .eigenstratgeno in
lockstep by row number, never rejoin by SNP ID" style already used by
pseudo_haploid_call.py (see that script's docstring point 1) -- keeping
the SAME invariant here rather than introducing a SNP-ID-keyed join.

WARNING (must appear in any downstream report, see execution plan
section 5): every ancient sample ends up using a DIFFERENT SNP subset,
so PC coordinates from one sample's subset run are NOT directly
comparable to another sample's -- do not plot multiple samples' subset
runs on one shared PC1/PC2 axis pair and compare positions.

Usage:
  python3 build_sample_panel_subset.py \\
    --panel-snp panel.snp \\
    --panel-geno panel.eigenstratgeno \\
    --ancient-calls sample.panel.pseudohap.txt \\
    --out-panel-snp sample.panel.subset.snp \\
    --out-panel-geno sample.panel.subset.eigenstratgeno \\
    --out-ancient-calls sample.panel.subset.calls.txt
"""

from __future__ import annotations

import argparse
import os
import sys
import tempfile
from pathlib import Path

VALID_CALLS = frozenset("0129")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--panel-snp", required=True)
    parser.add_argument("--panel-geno", required=True)
    parser.add_argument("--ancient-calls", required=True, help="pseudo_haploid_call.py --out file for this sample+panel")
    parser.add_argument("--out-panel-snp", required=True)
    parser.add_argument("--out-panel-geno", required=True)
    parser.add_argument("--out-ancient-calls", required=True)
    parser.add_argument("--report", default=None, help="optional TSV: kept/dropped SNP counts")
    return parser.parse_args(argv)


def validate_args(args: argparse.Namespace) -> None:
    for label, path in (
        ("--panel-snp", args.panel_snp),
        ("--panel-geno", args.panel_geno),
        ("--ancient-calls", args.ancient_calls),
    ):
        if not Path(path).is_file():
            raise FileNotFoundError(f"{label} file not found: {path}")


class AtomicWriter:
    """Write to a temp file in the target's own directory; rename in on success only."""

    def __init__(self, path: Path):
        self.path = path
        self.path.parent.mkdir(parents=True, exist_ok=True)
        fd, tmp_name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
        os.close(fd)
        self.tmp_path = Path(tmp_name)
        self.handle = self.tmp_path.open("w")

    def write(self, s: str) -> None:
        self.handle.write(s)

    def commit(self) -> None:
        self.handle.close()
        os.replace(self.tmp_path, self.path)

    def abort(self) -> None:
        self.handle.close()
        self.tmp_path.unlink(missing_ok=True)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        validate_args(args)

        writers = [
            AtomicWriter(Path(args.out_panel_snp)),
            AtomicWriter(Path(args.out_panel_geno)),
            AtomicWriter(Path(args.out_ancient_calls)),
        ]
        out_snp, out_geno, out_calls = writers

        n_total = n_kept = n_dropped = 0
        try:
            with open(args.panel_snp) as f_snp, open(args.panel_geno) as f_geno, open(args.ancient_calls) as f_calls:
                for line_no, (snp_line, geno_line, call_line) in enumerate(zip(f_snp, f_geno, f_calls), start=1):
                    call = call_line.strip()
                    if len(call) != 1 or call not in VALID_CALLS:
                        raise ValueError(
                            f"{args.ancient_calls}:{line_no}: expected a single 0/1/2/9 character, got {call_line!r}"
                        )
                    n_total += 1
                    if call == "9":
                        n_dropped += 1
                        continue
                    n_kept += 1
                    out_snp.write(snp_line if snp_line.endswith("\n") else snp_line + "\n")
                    out_geno.write(geno_line if geno_line.endswith("\n") else geno_line + "\n")
                    out_calls.write(call + "\n")

                # Confirm all three inputs were exactly the same length --
                # zip() silently stops at the shortest, so a length mismatch
                # would otherwise pass with no error and a truncated result.
                sentinel = object()
                leftover_snp = next(f_snp, sentinel)
                leftover_geno = next(f_geno, sentinel)
                leftover_calls = next(f_calls, sentinel)
                if not (leftover_snp is sentinel and leftover_geno is sentinel and leftover_calls is sentinel):
                    raise ValueError(
                        f"input files have different line counts (compared {n_total} rows via zip(), "
                        "at least one file has more) -- panel-snp/panel-geno/ancient-calls must all "
                        "have exactly the same number of lines"
                    )

            if n_kept == 0:
                raise ValueError(
                    f"{args.ancient_calls} has zero covered (non-9) sites against this panel -- "
                    "nothing to build a subset from"
                )

        except BaseException:
            for w in writers:
                w.abort()
            raise

        for w in writers:
            w.commit()

        if args.report:
            report_path = Path(args.report)
            report_path.parent.mkdir(parents=True, exist_ok=True)
            fd, tmp_name = tempfile.mkstemp(prefix=f".{report_path.name}.", suffix=".tmp", dir=report_path.parent)
            os.close(fd)
            tmp_path = Path(tmp_name)
            try:
                with tmp_path.open("w") as handle:
                    handle.write("metric\tvalue\n")
                    handle.write(f"panel_total_snps\t{n_total}\n")
                    handle.write(f"kept_snps\t{n_kept}\n")
                    handle.write(f"dropped_snps\t{n_dropped}\n")
                os.replace(tmp_path, report_path)
            except BaseException:
                tmp_path.unlink(missing_ok=True)
                raise

    except (OSError, ValueError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(
        f"[subset] {args.ancient_calls}: {n_kept}/{n_total} SNPs kept "
        f"({100 * n_kept / n_total:.3f}%), {n_dropped} dropped as uncovered/missing",
        file=sys.stderr,
    )
    print(f"[done] wrote {args.out_panel_snp}, {args.out_panel_geno}, {args.out_ancient_calls}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
