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

    raw = args.out + ".raw"
    subprocess.run(["plink2", "--bfile", args.bfile, "--export", "A", "--out", args.out],
                   check=True)

    calls = []
    with open(raw) as fh:
        header = fh.readline().rstrip("\n").split()
        n = len(header) - 2  # skip FID IID
        line = fh.readline().rstrip("\n").split()
        for val in line[2:2 + n]:
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
