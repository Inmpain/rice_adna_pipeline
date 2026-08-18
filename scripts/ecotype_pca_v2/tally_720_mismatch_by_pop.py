#!/usr/bin/env python3
"""For the 14 mismatch sites, tally 720-panel genotypes by population label.

Genotype codes (EIGENSTRAT): 0 = homozygous ALT, 1 = heterozygous,
2 = homozygous REF, 9 = missing. Population labels come from the .ind's
third column (OrA-OrF wild groups, OrADM, RAY, and cultivated anchors
IND/AUS/ARO/TRJ/TEJ/ADM).
"""
import gzip
import sys
from collections import Counter, defaultdict


def open_maybe_gz(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_labels(ind_path):
    labels = []
    with open_maybe_gz(ind_path) as fh:
        for line in fh:
            f = line.split()
            if len(f) >= 3:
                labels.append(f[2])
            elif len(f) == 2:
                labels.append(f[1])
            else:
                labels.append("?")
    return labels


def read_mismatch(mismatch_path):
    ids = []
    info = {}
    with open(mismatch_path) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if not f or not f[0]:
                continue
            sid = f[0]
            ids.append(sid)
            info[sid] = (f[1], f[2], f[3], f[4], f[5])
    return ids, info


def main(snp_path, geno_path, ind_path, mismatch_path, out_path):
    labels = read_labels(ind_path)
    n_ind = len(labels)
    target_ids, info = read_mismatch(mismatch_path)
    target_set = set(target_ids)

    # pass 1: find 1-based line number of each target SNP ID in the .snp
    line_for = {}
    with open_maybe_gz(snp_path) as fh:
        for line_no, line in enumerate(fh, 1):
            f = line.split()
            if len(f) >= 2 and f[0] in target_set:
                line_for[f[0]] = line_no

    target_lines = {ln: sid for sid, ln in line_for.items()}
    line_set = set(target_lines)

    # pass 2: stream the .geno, extract only the target rows
    rows = {}
    with open_maybe_gz(geno_path) as fh:
        for line_no, line in enumerate(fh, 1):
            if line_no in line_set:
                rows[target_lines[line_no]] = line.strip()

    with open(out_path, "w") as out:
        out.write("snp_id\tchrom\tpos\tref\talt\tfasta\tpop\tn_ind\tn_ALT0\tn_het1\tn_REF2\tn_miss9\n")
        for sid in target_ids:
            chrom, pos, ref, alt, fasta = info[sid]
            g = rows.get(sid, "")
            if len(g) != n_ind:
                print(f"WARNING: {sid}: genotype row length {len(g)} != n_ind {n_ind}", file=sys.stderr)
            counts = defaultdict(Counter)
            for i, ch in enumerate(g):
                lab = labels[i] if i < len(labels) else "?"
                counts[lab][ch] += 1
            for lab in sorted(counts):
                c = counts[lab]
                out.write(
                    f"{sid}\t{chrom}\t{pos}\t{ref}\t{alt}\t{fasta}\t{lab}\t{sum(c.values())}"
                    f"\t{c.get('0', 0)}\t{c.get('1', 0)}\t{c.get('2', 0)}\t{c.get('9', 0)}\n"
                )
    print(f"wrote {out_path}")


if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("usage: python3 tally_720_mismatch_by_pop.py SNP GENO IND MISMATCH OUT", file=sys.stderr)
        sys.exit(2)
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
