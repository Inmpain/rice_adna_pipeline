#!/usr/bin/env python3
"""Tally match_ref / match_alt / mismatch by chromosome and by allele pair
for the 720 panel vs irgsp.fa. Read-only; stdlib only."""
import sys
from collections import defaultdict


def load_fasta(path):
    seqs = {}
    chrom = None
    chunks = []
    with open(path) as fh:
        for line in fh:
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
    return seqs


def main(snp_path, fasta_path):
    seqs = load_fasta(fasta_path)
    by_chrom = defaultdict(lambda: {"match_ref": 0, "match_alt": 0, "mismatch": 0})
    by_pair = defaultdict(lambda: {"match_ref": 0, "match_alt": 0, "mismatch": 0})
    total = 0
    with open(snp_path) as fh:
        for line in fh:
            f = line.split()
            if len(f) != 6:
                continue
            _id, chrom, _gp, pos_text, ref, alt = f
            ref = ref.upper()
            alt = alt.upper()
            try:
                pos = int(pos_text)
            except ValueError:
                continue
            try:
                contig = "chr%02d" % int(chrom.lstrip("0") or "0")
            except ValueError:
                contig = chrom
            seq = seqs.get(contig)
            if seq is None:
                continue
            idx = pos - 1
            if idx < 0 or idx >= len(seq):
                continue
            fa = seq[idx]
            total += 1
            pair = "/".join(sorted([ref, alt]))
            if fa == ref:
                by_chrom[contig]["match_ref"] += 1
                by_pair[pair]["match_ref"] += 1
            elif fa == alt:
                by_chrom[contig]["match_alt"] += 1
                by_pair[pair]["match_alt"] += 1
            else:
                by_chrom[contig]["mismatch"] += 1
                by_pair[pair]["mismatch"] += 1

    print(f"total\t{total}")
    print("=== by chromosome ===")
    print("chrom\tmatch_ref\tmatch_alt\tmismatch\talt_frac")
    for chrom in sorted(by_chrom):
        c = by_chrom[chrom]
        n = c["match_ref"] + c["match_alt"] + c["mismatch"]
        print(f"{chrom}\t{c['match_ref']}\t{c['match_alt']}\t{c['mismatch']}\t{c['match_alt'] / max(n, 1):.4f}")

    print("=== by allele pair (sorted by match_alt desc) ===")
    print("pair\tmatch_ref\tmatch_alt\tmismatch\talt_frac")
    for pair in sorted(by_pair, key=lambda p: -by_pair[p]["match_alt"]):
        c = by_pair[pair]
        n = c["match_ref"] + c["match_alt"] + c["mismatch"]
        print(f"{pair}\t{c['match_ref']}\t{c['match_alt']}\t{c['mismatch']}\t{c['match_alt'] / max(n, 1):.4f}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("usage: python3 tally_720_refalt.py SNP FASTA", file=sys.stderr)
        sys.exit(2)
    main(sys.argv[1], sys.argv[2])
