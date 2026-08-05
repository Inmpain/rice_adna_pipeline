#!/usr/bin/env python3
import sys, gzip, random, re
import pysam

def open_maybe_gz(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")

def parse_line(parts, fmt):
    if fmt == "bim":
        chrom, snpid, gd, pos, a1, a2 = parts[:6]
    else:
        snpid, chrom, gd, pos, a1, a2 = parts[:6]
    return chrom, pos, a1, a2

def normalize_chrom(chrom):
    c = re.sub(r'(?i)^chr', '', chrom)
    c = c.lstrip('0')
    return c if c else '0'

def main():
    path = sys.argv[1]
    fmt = sys.argv[2]
    N = int(sys.argv[3]) if len(sys.argv) > 3 else 200
    REF = "/home/scratch/yinmt202607/db/asian_rice_panel_index/irgsp.fa"

    fasta = pysam.FastaFile(REF)
    fasta_contigs = set(fasta.references)

    raw_chrom_forms = set()
    sites = []
    skipped_chrom = set()

    with open_maybe_gz(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 6:
                continue
            chrom, pos, a1, a2 = parse_line(parts, fmt)
            raw_chrom_forms.add(chrom)
            if a1 in ("0", "N") or a2 in ("0", "N"):
                continue
            norm = normalize_chrom(chrom)
            try:
                chrom_num = int(norm)
            except ValueError:
                skipped_chrom.add(chrom)
                continue
            chrom_fa = "chr%02d" % chrom_num
            if chrom_fa not in fasta_contigs:
                skipped_chrom.add(chrom)
                continue
            sites.append((chrom_fa, int(pos), a1, a2))

    print("raw chrom forms (sample):", sorted(raw_chrom_forms)[:20])
    if skipped_chrom:
        print("skipped chrom names:", skipped_chrom)

    random.seed(42)
    sample = random.sample(sites, min(N, len(sites)))

    n_a2_ref = 0
    n_a1_ref = 0
    n_neither = 0
    mismatches = []

    for chrom_fa, pos, a1, a2 in sample:
        ref_base = fasta.fetch(chrom_fa, pos - 1, pos).upper()
        if ref_base == a2:
            n_a2_ref = n_a2_ref + 1
        elif ref_base == a1:
            n_a1_ref = n_a1_ref + 1
            mismatches.append((chrom_fa, pos, a1, a2, ref_base))
        else:
            n_neither = n_neither + 1
            mismatches.append((chrom_fa, pos, a1, a2, ref_base))

    print("checked sites:", len(sample))
    print("A2_is_ref:", n_a2_ref)
    print("A1_is_ref:", n_a1_ref)
    print("neither:", n_neither)

    if mismatches:
        print("mismatches (up to 20):")
        for m in mismatches[:20]:
            print(m)

if __name__ == "__main__":
    main()
