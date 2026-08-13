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

--mask-from, AND WHY IT EXISTS (found via a real leave-one-out smoke
test, not anticipated in advance): by default, which rows get kept is
decided by --ancient-calls' OWN non-9 positions, and that file's own
values are what get written out -- correct for a real ancient sample.
But simulate_leaveoneout_projection.py's output is a REAL modern
individual's genotype masked down to an ancient sample's coverage
pattern, and that modern individual has its own pre-existing missingness
too -- naively subsetting it by its OWN non-9 positions (default
behavior) drops additional rows beyond what the ancient sample's mask
already dropped, so the "held-out individual" subset panel ends up with
FEWER, misaligned rows versus the real ancient sample's subset panel
built from the same mask (confirmed on a real run: 147 masked-in sites,
but only 123 survived because 24 were independently missing in the
held-out individual's own genotype) -- merge_ancient_into_panel.py
requires every column's calls file to have the exact same row count and
order, so these two subsets are then NOT mergeable into one matrix. Pass
--mask-from <the real ancient sample's calls file> together with
--ancient-calls <the simulated individual's calls file> to keep rows
using the MASK's non-9 positions while still writing out whatever the
ancient-calls file's actual value is there (including a legitimate "9"
if the held-out individual happens to be missing at that particular
site) -- this guarantees the two subset panels share identical rows.
Omitting --mask-from (the default) makes mask-from == ancient-calls,
reproducing the original single-file behavior exactly.

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
    parser.add_argument("--ancient-calls", required=True, help="pseudo_haploid_call.py --out file for this sample+panel (or simulate_leaveoneout_projection.py's --out)")
    parser.add_argument(
        "--mask-from",
        default=None,
        help="optional: decide kept rows from THIS file's non-9 positions instead of --ancient-calls' own "
        "(leave-one-out: pass the real ancient sample's calls here so the held-out individual's subset "
        "aligns row-for-row with it); defaults to --ancient-calls itself",
    )
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
    if args.mask_from is not None and not Path(args.mask_from).is_file():
        raise FileNotFoundError(f"--mask-from file not found: {args.mask_from}")


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

        mask_path = args.mask_from or args.ancient_calls

        n_total = n_kept = n_dropped = n_kept_but_call_missing = 0
        try:
            with open(args.panel_snp) as f_snp, open(args.panel_geno) as f_geno, \
                    open(args.ancient_calls) as f_calls, open(mask_path) as f_mask:
                for line_no, (snp_line, geno_line, call_line, mask_line) in enumerate(
                    zip(f_snp, f_geno, f_calls, f_mask), start=1
                ):
                    call = call_line.strip()
                    if len(call) != 1 or call not in VALID_CALLS:
                        raise ValueError(
                            f"{args.ancient_calls}:{line_no}: expected a single 0/1/2/9 character, got {call_line!r}"
                        )
                    mask_call = mask_line.strip()
                    if len(mask_call) != 1 or mask_call not in VALID_CALLS:
                        raise ValueError(
                            f"{mask_path}:{line_no}: expected a single 0/1/2/9 character, got {mask_line!r}"
                        )
                    n_total += 1
                    if mask_call == "9":
                        n_dropped += 1
                        continue
                    n_kept += 1
                    if call == "9":
                        # Only possible when --mask-from differs from
                        # --ancient-calls: the mask says this row is
                        # "covered", but the file we're extracting values
                        # FROM (e.g. a held-out modern individual in a
                        # leave-one-out run) happens to be missing here in
                        # its own real genotype. Kept anyway, on purpose --
                        # see the --mask-from docstring section.
                        n_kept_but_call_missing += 1
                    out_snp.write(snp_line if snp_line.endswith("\n") else snp_line + "\n")
                    out_geno.write(geno_line if geno_line.endswith("\n") else geno_line + "\n")
                    out_calls.write(call + "\n")

                # Confirm all four inputs were exactly the same length --
                # zip() silently stops at the shortest, so a length mismatch
                # would otherwise pass with no error and a truncated result.
                sentinel = object()
                leftovers = [next(f, sentinel) for f in (f_snp, f_geno, f_calls, f_mask)]
                if not all(x is sentinel for x in leftovers):
                    raise ValueError(
                        f"input files have different line counts (compared {n_total} rows via zip(), "
                        "at least one file has more) -- panel-snp/panel-geno/ancient-calls/mask-from must "
                        "all have exactly the same number of lines"
                    )

            if n_kept == 0:
                raise ValueError(
                    f"{mask_path} has zero covered (non-9) sites against this panel -- "
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
                    handle.write(f"mask_from\t{mask_path}\n")
                    handle.write(f"ancient_calls\t{args.ancient_calls}\n")
                    handle.write(f"panel_total_snps\t{n_total}\n")
                    handle.write(f"kept_snps\t{n_kept}\n")
                    handle.write(f"dropped_snps\t{n_dropped}\n")
                    handle.write(f"kept_but_call_missing\t{n_kept_but_call_missing}\n")
                os.replace(tmp_path, report_path)
            except BaseException:
                tmp_path.unlink(missing_ok=True)
                raise

    except (OSError, ValueError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(
        f"[subset] mask={mask_path} calls={args.ancient_calls}: {n_kept}/{n_total} SNPs kept "
        f"({100 * n_kept / n_total:.3f}%), {n_dropped} dropped by the mask, "
        f"{n_kept_but_call_missing} of the kept rows are '9' in ancient-calls itself",
        file=sys.stderr,
    )
    print(f"[done] wrote {args.out_panel_snp}, {args.out_panel_geno}, {args.out_ancient_calls}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
