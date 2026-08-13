#!/usr/bin/env python3
"""
Drop individuals with a given population label entirely from a panel's
.ind file AND its .eigenstratgeno matrix (not just from a smartpca
poplistname).

WHY DELETE INSTEAD OF JUST OMITTING FROM poplistname: smartpca
-lsqproject projects onto the fixed axes ANY individual in the .ind file
that ISN'T listed in poplistname (that's the exact mechanism used to
project the ancient samples themselves, see merge_ancient_into_panel.py's
docstring). Leaving unwanted individuals in .ind but out of poplistname
does NOT make them disappear from the projected result -- they'd still
show up on the final plot as unlabeled projected points, indistinguishable
from real ancient samples in smartpca's own logic. Actually removing them
from both .ind and .eigenstratgeno is the only way to make them not
appear at all.

WHY cut -c, NOT A PYTHON PER-LINE LOOP: an EIGENSTRAT .eigenstratgeno
file has one row per SNP and one character per individual, no
separators, with a row count in the tens of millions for a panel like
NB_final_snp (29,635,224 SNPs) -- a Python loop touching every character
of every row would take far too long. Dropping only a handful of
individuals out of thousands leaves a small number of large contiguous
"keep" ranges (e.g. dropping individuals at columns 6, 401, 1002 out of
3024 gives keep-ranges "1-5,7-400,402-1001,1003-3024") -- cut -c is a
compiled tool built exactly for this and was already this project's own
established, debugged approach for extracting columns from this same
matrix family (docs/3krgp_integration_and_simulation_prep.md section 3:
"用cut -c而非awk内部getline(更可靠)").

This script only COMPUTES the keep-ranges and .ind filtering and RUNS
cut (via subprocess) if --geno-in/--geno-out are both given; on a ~90GB
file this takes real wall-clock time, so run it via sbatch, not
interactively in a login-node foreground shell.

Usage:
  python3 filter_panel_by_label.py \\
    --ind /home/scratch/yinmt202607/db/29M_3k/NB_final_snp.labeled.ind \\
    --drop-label UNK \\
    --ind-out /home/scratch/yinmt202607/db/29M_3k/NB_final_snp.filtered.ind \\
    --geno-in /home/scratch/yinmt202607/db/29M_3k/NB_final_snp.eigenstratgeno \\
    --geno-out /home/scratch/yinmt202607/db/29M_3k/NB_final_snp.filtered.eigenstratgeno
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--ind", required=True, help="path to the labeled .ind file")
    parser.add_argument(
        "--drop-label",
        action="append",
        required=True,
        help="population label to drop entirely (repeatable for multiple labels)",
    )
    parser.add_argument("--ind-out", required=True, help="path for the filtered .ind file")
    parser.add_argument("--geno-in", default=None, help="path to the panel's .eigenstratgeno file (optional)")
    parser.add_argument("--geno-out", default=None, help="path for the filtered .eigenstratgeno file (required if --geno-in given)")
    parser.add_argument(
        "--print-cut-cmd-only",
        action="store_true",
        help="print the cut command instead of running it (for manual sbatch submission)",
    )
    return parser.parse_args(argv)


def validate_args(args: argparse.Namespace) -> None:
    if not Path(args.ind).is_file():
        raise FileNotFoundError(f".ind file not found: {args.ind}")
    if args.geno_in is not None:
        if not Path(args.geno_in).is_file():
            raise FileNotFoundError(f".eigenstratgeno file not found: {args.geno_in}")
        if args.geno_out is None:
            raise ValueError("--geno-out is required when --geno-in is given")


def load_ind(ind_path: Path) -> list[tuple[str, str, str]]:
    rows: list[tuple[str, str, str]] = []
    with ind_path.open() as handle:
        for line_no, line in enumerate(handle, start=1):
            fields = line.split()
            if len(fields) != 3:
                raise ValueError(f"{ind_path}:{line_no}: expected 3 fields, got {len(fields)}: {line!r}")
            rows.append(tuple(fields))
    return rows


def compute_keep_ranges(total: int, drop_positions: set[int]) -> str:
    """drop_positions is 1-indexed. Returns a cut -c range spec covering everything else."""
    ranges: list[str] = []
    start = None
    for pos in range(1, total + 1):
        if pos in drop_positions:
            if start is not None:
                ranges.append(f"{start}-{pos - 1}" if start != pos - 1 else str(start))
                start = None
        else:
            if start is None:
                start = pos
    if start is not None:
        ranges.append(f"{start}-{total}" if start != total else str(start))
    if not ranges:
        raise ValueError("keep-ranges is empty -- every individual would be dropped, refusing to proceed")
    return ",".join(ranges)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        validate_args(args)
        ind_rows = load_ind(Path(args.ind))
        drop_labels = set(args.drop_label)

        drop_positions: set[int] = set()
        keep_rows: list[tuple[str, str, str]] = []
        for i, (sample_id, sex, label) in enumerate(ind_rows, start=1):
            if label in drop_labels:
                drop_positions.add(i)
            else:
                keep_rows.append((sample_id, sex, label))

        if not drop_positions:
            raise ValueError(f"none of {sorted(drop_labels)} matched any label in {args.ind} -- check spelling")

        out_ind = Path(args.ind_out)
        out_ind.parent.mkdir(parents=True, exist_ok=True)
        with out_ind.open("w") as handle:
            for sample_id, sex, label in keep_rows:
                handle.write(f"{sample_id:>30} {sex} {label}\n")

        keep_spec = compute_keep_ranges(len(ind_rows), drop_positions)

        print(f"[filter] {args.ind}: {len(ind_rows)} individuals, dropping {len(drop_positions)} with label(s) {sorted(drop_labels)}", file=sys.stderr)
        print(f"[filter] kept {len(keep_rows)} individuals -> wrote {out_ind}", file=sys.stderr)
        print(f"[filter] cut -c keep-spec ({keep_spec.count(',') + 1} ranges): {keep_spec}", file=sys.stderr)

        if args.geno_in is None:
            print("[filter] --geno-in not given, skipping the .eigenstratgeno filtering step", file=sys.stderr)
            print(f"[filter] to filter it later, run: cut -c{keep_spec} {args.ind}_placeholder > OUTPUT", file=sys.stderr)
            return 0

        cut_cmd = ["cut", f"-c{keep_spec}", args.geno_in]
        if args.print_cut_cmd_only:
            print("[filter] would run:", " ".join(cut_cmd), ">", args.geno_out, file=sys.stderr)
            return 0

        print(f"[filter] running cut on {args.geno_in} (this is a large file, expect real wall-clock time)...", file=sys.stderr)
        geno_out = Path(args.geno_out)
        geno_out.parent.mkdir(parents=True, exist_ok=True)
        with geno_out.open("w") as out_handle:
            result = subprocess.run(cut_cmd, stdout=out_handle, check=False)
        if result.returncode != 0:
            geno_out.unlink(missing_ok=True)
            raise RuntimeError(f"cut exited with code {result.returncode}")

        print(f"[filter] wrote {geno_out}, verifying row widths on the first line...", file=sys.stderr)
        with geno_out.open() as handle:
            first_line = handle.readline().rstrip("\n")
        if len(first_line) != len(keep_rows):
            geno_out.unlink(missing_ok=True)
            raise RuntimeError(
                f"first output row has {len(first_line)} characters, expected {len(keep_rows)} "
                "-- output removed, something is wrong with the keep-spec or input file"
            )
        print(f"[filter] first-row width check passed ({len(first_line)} characters)", file=sys.stderr)

    except (OSError, ValueError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print("[done]", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
