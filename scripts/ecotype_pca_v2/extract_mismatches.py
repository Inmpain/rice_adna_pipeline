#!/usr/bin/env python3
"""Write mismatch sites (FASTA base is neither declared REF nor ALT) as TSV."""
import sys


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


def main(snp_path, fasta_path, out_path):
    seqs = load_fasta(fasta_path)
    n = 0
    with open(snp_path) as fh, open(out_path, "w") as out:
        out.write("snp_id\tchrom\tpos\tref\talt\tfasta_base\n")
        for line in fh:
            f = line.split()
            if len(f) != 6:
                continue
            sid, chrom, _gp, pos_text, ref, alt = f
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
            if fa != ref and fa != alt:
                out.write(f"{sid}\t{chrom}\t{pos}\t{ref}\t{alt}\t{fa}\n")
                n += 1
    print(f"wrote {n} mismatches to {out_path}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("usage: python3 extract_mismatches.py SNP FASTA OUT", file=sys.stderr)
        sys.exit(2)
    main(sys.argv[1], sys.argv[2], sys.argv[3])
