#!/usr/bin/env python3
"""Per-sample coverage funnel: raw-panel covered -> MAF-passed -> LD-passed.

Inputs:
  ancient_union_sites.tsv  (from 19_survey_ancient_coverage.py; has samples_covered)
  MAF snplist              (from geno_maf_filtered.bim, column 2)
  LD fixed snplist         (from 27_ancient_coverage_first_ld_prune.py)
Output:
  one row per ancient sample + a UNION row, columns raw / maf / ld.
"""
import sys
from collections import Counter


def load_sites(tsv_path):
    sites = {}
    with open(tsv_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(header):
                continue
            d = dict(zip(header, f))
            sc = d.get("samples_covered", "")
            sites[d["snp_id"]] = set(sc.split(",")) if sc else set()
    return sites


def load_ids(path):
    return set(line.strip() for line in open(path) if line.strip())


def main(union_tsv, maf_snplist, ld_snplist, out_path):
    sites = load_sites(union_tsv)
    maf = load_ids(maf_snplist)
    ld = load_ids(ld_snplist)

    samples = sorted(set().union(*sites.values()))
    per = {s: Counter() for s in samples}
    union = Counter()

    for sid, covered in sites.items():
        union["raw"] += 1
        for s in covered:
            per[s]["raw"] += 1
        if sid in maf:
            union["maf"] += 1
            for s in covered:
                per[s]["maf"] += 1
        if sid in ld:
            union["ld"] += 1
            for s in covered:
                per[s]["ld"] += 1

    with open(out_path, "w") as out:
        out.write("sample\traw_panel_covered\tmaf_passed_covered\tld_passed_covered\n")
        for s in samples:
            out.write(f"{s}\t{per[s]['raw']}\t{per[s]['maf']}\t{per[s]['ld']}\n")
        out.write(f"UNION\t{union['raw']}\t{union['maf']}\t{union['ld']}\n")
    print(f"wrote {out_path}")


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("usage: coverage_funnel.py UNION_TSV MAF_SNPLIST LD_SNPLIST OUT", file=sys.stderr)
        sys.exit(2)
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
