#!/usr/bin/env python3
"""Write the list of SNP IDs whose REF/ALT are flipped relative to the
reference FASTA (FASTA base == declared ALT, not REF).

For clean-swap panels (3K / Civán) this is "all sites" (so no need to write it);
for the mixed 720 panel it is the ~605k flipped subset, which Phase B needs.
Read-only, stdlib only.
"""
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
    n_ref = n_alt = n_mismatch = 0
    with open(snp_path) as fh, open(out_path, "w") as out:
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
            if fa == ref:
                n_ref += 1
            elif fa == alt:
                n_alt += 1
                out.write(sid + "\n")
            else:
                n_mismatch += 1
    print(f"match_ref={n_ref} match_alt(flip)={n_alt} mismatch={n_mismatch}")
    print(f"wrote {n_alt} flip IDs to {out_path}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("usage: python3 make_flip_list.py SNP FASTA OUT", file=sys.stderr)
        sys.exit(2)
    main(sys.argv[1], sys.argv[2], sys.argv[3])
