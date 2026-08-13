#!/usr/bin/env python3
"""
Simulate a "digital ancient sample" from one modern individual already in
a panel, by masking it down to another (real ancient) sample's covered-
site pattern and re-calling pseudo-haploid genotypes at those sites --
the leave-one-out positive control from
docs/ECOTYPE_PCA_EXECUTION_PLAN.md section 6 / pseudo_haploid_call.py
docstring point 4.

WHY THIS IS THE AUTHORITATIVE CHECK, NOT check_ref.py: a modern
individual already in the panel has a KNOWN population label. If this
script's output -- fed through the exact same
build_sample_panel_subset.py -> merge_ancient_into_panel.py -> smartpca
-lsqproject chain used for real ancient samples -- projects back near
that individual's own known population, the whole chain (including
pseudo_haploid_call.py's REF/ALT 0/2 convention) is validated end to
end. If the encoding were backwards, this simulated individual would
project systematically AWAY from its true population, not just noisier
-- a far stronger signal than any FASTA match-rate spot check.

OUTPUT FORMAT MATCHES pseudo_haploid_call.py's --out EXACTLY (one 0/1/2/9
call per line, same row order as the panel's .snp file) so it is a
drop-in input to build_sample_panel_subset.py -- that script does not
need to know or care whether a calls file came from a real BAM or this
simulation.

WHAT THIS SCRIPT DOES NOT DO -- READ BEFORE USING: it does not touch the
panel's .ind file. Whoever runs the downstream smartpca step MUST give
--held-out-sample a population label that poplistname excludes for that
specific run (e.g. a one-off "LOO_HELDOUT" label), while every OTHER
individual sharing its true population keeps their real label and stays
in poplistname -- otherwise this individual's own real genotype column
would still help build the axis it's being projected onto, making the
test circular. smartpca's poplistname only filters by population label,
not by individual, so this has to be a per-run .ind edit done by the
caller; it is deliberately out of scope here (same separation of
concerns as pseudo_haploid_call.py not touching .ind either).

PSEUDO-HAPLOID SIMULATION LOGIC, per masked (covered) site: the modern
individual's own EIGENSTRAT genotype at that site is homozygous-ref (2),
homozygous-alt (0), or missing (9) -- output unchanged, there is only
one possible base to "read". A heterozygous (1) site has two possible
alleles; output a uniform-random choice of 0 or 2, mirroring
pseudo_haploid_call.py's own "one read drawn at random becomes the call"
rule for a real ancient sample at a truly heterozygous underlying site.

Usage:
  python3 simulate_leaveoneout_projection.py \\
    --panel-geno panel.eigenstratgeno \\
    --panel-ind panel.ind \\
    --held-out-sample SOME_KNOWN_SAMPLE_ID \\
    --mask-from real_ancient_sample.panel.pseudohap.txt \\
    --seed 0 \\
    --out heldout_sample.simulated.pseudohap.txt \\
    --report heldout_sample.simulated.report.tsv
"""

from __future__ import annotations

import argparse
import os
import random
import sys
import tempfile
from pathlib import Path

VALID_GENO = frozenset("0129")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--panel-geno", required=True)
    parser.add_argument("--panel-ind", required=True)
    parser.add_argument("--held-out-sample", required=True, help="sample ID (column 1 of --panel-ind) to mask and re-call")
    parser.add_argument("--mask-from", required=True, help="pseudo_haploid_call.py --out file whose covered/missing PATTERN to replicate")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--out", required=True)
    parser.add_argument("--report", default=None)
    return parser.parse_args(argv)


def validate_args(args: argparse.Namespace) -> None:
    for label, path in (("--panel-geno", args.panel_geno), ("--panel-ind", args.panel_ind), ("--mask-from", args.mask_from)):
        if not Path(path).is_file():
            raise FileNotFoundError(f"{label} file not found: {path}")


def find_column_index(ind_path: Path, sample_id: str) -> int:
    """0-based column index of sample_id within the panel's .eigenstratgeno rows."""
    matches: list[int] = []
    with ind_path.open() as handle:
        for i, line in enumerate(handle):
            fields = line.split()
            if not fields:
                continue
            if fields[0] == sample_id:
                matches.append(i)
    if not matches:
        raise ValueError(f"sample {sample_id!r} not found in {ind_path}")
    if len(matches) > 1:
        raise ValueError(f"sample {sample_id!r} appears {len(matches)} times in {ind_path} -- expected a unique ID")
    return matches[0]


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        validate_args(args)
        random.seed(args.seed)

        col_index = find_column_index(Path(args.panel_ind), args.held_out_sample)
        print(f"[loo] {args.held_out_sample} is column {col_index + 1} (1-based) of {args.panel_geno}", file=sys.stderr)

        n_total = n_masked_out = n_modern_missing = n_het_simulated = n_homo_kept = 0

        out_path = Path(args.out)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fd, tmp_name = tempfile.mkstemp(prefix=f".{out_path.name}.", suffix=".tmp", dir=out_path.parent)
        os.close(fd)
        tmp_path = Path(tmp_name)

        try:
            with open(args.panel_geno) as f_geno, open(args.mask_from) as f_mask, tmp_path.open("w") as f_out:
                for line_no, (geno_line, mask_line) in enumerate(zip(f_geno, f_mask), start=1):
                    n_total += 1
                    mask_call = mask_line.strip()
                    if len(mask_call) != 1 or mask_call not in VALID_GENO:
                        raise ValueError(f"{args.mask_from}:{line_no}: expected a single 0/1/2/9 character, got {mask_line!r}")

                    if mask_call == "9":
                        n_masked_out += 1
                        f_out.write("9\n")
                        continue

                    geno_row = geno_line.rstrip("\n")
                    if col_index >= len(geno_row):
                        raise ValueError(
                            f"{args.panel_geno}:{line_no}: row has only {len(geno_row)} characters, "
                            f"need index {col_index} for {args.held_out_sample} -- .panel-ind/.panel-geno mismatch?"
                        )
                    modern_call = geno_row[col_index]
                    if modern_call not in VALID_GENO:
                        raise ValueError(f"{args.panel_geno}:{line_no}: unexpected genotype character {modern_call!r} at column {col_index + 1}")

                    if modern_call == "9":
                        n_modern_missing += 1
                        f_out.write("9\n")
                    elif modern_call == "1":
                        n_het_simulated += 1
                        f_out.write(random.choice(("0", "2")) + "\n")
                    else:
                        n_homo_kept += 1
                        f_out.write(modern_call + "\n")

                sentinel = object()
                leftover_geno = next(f_geno, sentinel)
                leftover_mask = next(f_mask, sentinel)
                if not (leftover_geno is sentinel and leftover_mask is sentinel):
                    raise ValueError(
                        f"{args.panel_geno} and {args.mask_from} have different line counts "
                        f"(compared {n_total} rows via zip(), at least one file has more)"
                    )

            os.replace(tmp_path, out_path)
        except BaseException:
            tmp_path.unlink(missing_ok=True)
            raise

        n_called = n_homo_kept + n_het_simulated
        if args.report:
            report_path = Path(args.report)
            report_path.parent.mkdir(parents=True, exist_ok=True)
            fd, tmp_name = tempfile.mkstemp(prefix=f".{report_path.name}.", suffix=".tmp", dir=report_path.parent)
            os.close(fd)
            rtmp_path = Path(tmp_name)
            try:
                with rtmp_path.open("w") as handle:
                    handle.write("metric\tvalue\n")
                    handle.write(f"held_out_sample\t{args.held_out_sample}\n")
                    handle.write(f"mask_from\t{args.mask_from}\n")
                    handle.write(f"seed\t{args.seed}\n")
                    handle.write(f"total_panel_snps\t{n_total}\n")
                    handle.write(f"masked_out_by_ancient_pattern\t{n_masked_out}\n")
                    handle.write(f"modern_missing_at_masked_sites\t{n_modern_missing}\n")
                    handle.write(f"het_simulated\t{n_het_simulated}\n")
                    handle.write(f"homo_kept\t{n_homo_kept}\n")
                    handle.write(f"called\t{n_called}\n")
                os.replace(rtmp_path, report_path)
            except BaseException:
                rtmp_path.unlink(missing_ok=True)
                raise

    except (OSError, ValueError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(
        f"[loo] {args.held_out_sample}: {n_total} panel SNPs, "
        f"{n_total - n_masked_out} masked-in by {args.mask_from}'s coverage pattern, "
        f"{n_called} called ({n_homo_kept} homozygous copied, {n_het_simulated} het->random draw), "
        f"{n_modern_missing} missing in the modern sample itself at masked-in sites",
        file=sys.stderr,
    )
    print(f"[done] wrote {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
