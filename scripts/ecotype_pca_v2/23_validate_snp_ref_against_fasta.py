#!/usr/bin/env python3
"""Validate that a .snp file's declared REF allele actually matches the
reference FASTA sequence at that position.

Adapted from the same-lab Snakefile.pseudohaploid.from_panel's
validate_bim_a2_against_reference rule (2026-08-17 review): our own
pipeline had no equivalent check before this. fixed_projection_lib.iter_snp
only validates that REF/ALT are ACGT and distinct from each other -- it
never checks REF against the actual reference genome, so a strand-flip or
an off-by-one coordinate bug in the source panel would pass silently
through every downstream script.

Writes a small summary report (checked_sites / mismatch_n / up to
--max-report examples), not one row per SNP -- civan_snp.snp alone has
2.36M rows, a per-site TSV would be unusable and is not needed: either the
panel's REF alleles agree with the reference genome or they don't.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from fixed_projection_lib import format_contig, iter_panel_snp, refuse_existing


def load_fasta(path: str) -> dict[str, str]:
    seqs: dict[str, str] = {}
    chrom = None
    chunks: list[str] = []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if chrom is not None:
                    seqs[chrom] = "".join(chunks).upper()
                chrom = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if chrom is not None:
            seqs[chrom] = "".join(chunks).upper()
    if not seqs:
        raise ValueError(f"{path}: no sequences parsed (empty or malformed FASTA)")
    return seqs


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snp", required=True, help="raw upstream panel .snp, e.g. civan_snp.snp (6-column)")
    parser.add_argument("--fasta", required=True, help="reference FASTA the ancient BAMs were mapped against")
    parser.add_argument("--contig-format", default="chr%02d")
    parser.add_argument("--max-report", type=int, default=30)
    parser.add_argument("--out", required=True, help="small summary report path")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        refuse_existing([args.out], args.overwrite)
        seqs = load_fasta(args.fasta)
        checked = 0
        mismatches: list[str] = []
        for record in iter_panel_snp(args.snp):
            checked += 1
            contig = format_contig(args.contig_format, record["chrom"])
            seq = seqs.get(contig)
            if seq is None:
                mismatches.append(f"{record['id']}: contig {contig!r} not found in FASTA")
                continue
            idx = record["pos"] - 1
            if idx < 0 or idx >= len(seq):
                mismatches.append(f"{record['id']}: position {record['pos']} outside contig {contig} (len={len(seq)})")
                continue
            fasta_base = seq[idx]
            if fasta_base != record["ref"]:
                mismatches.append(
                    f"{record['id']}: declared REF={record['ref']} but FASTA {contig}:{record['pos']}={fasta_base}"
                )
        with open(args.out, "w") as out:
            out.write(f"checked_sites\t{checked}\n")
            out.write(f"mismatch_n\t{len(mismatches)}\n")
            out.write(f"fasta\t{args.fasta}\n")
            out.write(f"snp\t{args.snp}\n")
            for line in mismatches[: args.max_report]:
                out.write(f"mismatch_example\t{line}\n")
        if mismatches:
            print(
                f"FATAL: {len(mismatches)}/{checked} sites have REF mismatched against {args.fasta}; "
                f"see {args.out} for the first {min(len(mismatches), args.max_report)} example(s)",
                file=sys.stderr,
            )
            return 3
    except (OSError, ValueError) as exc:
        print(f"FATAL: {exc}", file=sys.stderr)
        return 3
    print(f"PASS: {checked} sites checked against {args.fasta}, 0 REF mismatches")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
