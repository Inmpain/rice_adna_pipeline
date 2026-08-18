#!/usr/bin/env python3
"""Convert one sample's pileupCaller PLINK output to a .calls.txt.

Output is one 0/2/9 per line, in the sample .bim order. pileupCaller --randomHaploid
produces haploid calls, so 0=hom ALT, 2=hom REF, 9=missing (any '1'/NA also -> 9).
"""
import argparse
import subprocess


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bfile", required=True, help="sample PLINK prefix (.bed/.bim/.fam)")
    ap.add_argument("--out", required=True, help="output prefix (writes OUT.calls.txt)")
    args = ap.parse_args()

    # Read the .bim marker order first, so we are robust to extra .raw columns
    # (plink2 --export A may include non-variant columns like PHENO/SEX).
    bim_ids = []
    with open(args.bfile + ".bim") as fh:
        for line in fh:
            bim_ids.append(line.split()[1])

    raw = args.out + ".raw"
    subprocess.run(["plink2", "--bfile", args.bfile, "--export", "A", "--out", args.out],
                   check=True)

    with open(raw) as fh:
        header = fh.readline().rstrip("\n").split()
        line = fh.readline().rstrip("\n").split()
    # Map variant column name -> value index (0-based within the data row).
    col = {name: i for i, name in enumerate(header)}
    calls = []
    for sid in bim_ids:
        if sid not in col:
            calls.append("9")
            continue
        val = line[col[sid]]
        if val == "0":
            calls.append("0")
        elif val == "2":
            calls.append("2")
        else:
            calls.append("9")  # '1' (unexpected for haploid) or NA -> missing

    out_calls = args.out + ".calls.txt"
    with open(out_calls, "w") as fh:
        fh.write("\n".join(calls) + "\n")
    print(f"wrote {out_calls}: {len(calls)} markers")


if __name__ == "__main__":
    main()
