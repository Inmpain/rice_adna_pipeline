#!/usr/bin/env python3
"""
Merge one or more ancient samples' pseudo-haploid genotype calls
(pseudo_haploid_call.py output) into a panel's own EIGENSTRAT genotype
file, as extra individual columns.

WHY THIS SHAPE, NOT SEPARATE FILES: this is the standard smartpca
-lsqproject setup (see ECOTYPE_PCA_PANEL.md section 2 / 5.2). smartpca
takes ONE genotype/snp/ind triple. Eigenvectors are computed only from
individuals whose population label (last column of .ind) appears in the
poplistname file passed to smartpca; any individual NOT in that list is
projected onto the resulting fixed axes without influencing them. So the
ancient samples must live in the SAME genotype file as the modern panel,
just given a population label that poplistname deliberately excludes --
there is no other mechanism to keep the modern eigenvectors fixed while
projecting new samples.

EFFICIENCY: merges ALL ancient samples for a given panel in ONE pass
over the panel's .eigenstratgeno file, not one pass per sample -- for
29M_3k that file has ~29.6 million lines, so re-reading it once per
sample (16 times) would be wasteful compared to once total.

SAFETY (2026-08-12 hardening, see docs/ECOTYPE_PCA_EXECUTION_PLAN.md):
the original version of this script could silently write a corrupt
merged file (or crash with a raw IndexError past the point where a
warning would ever print) if the panel's .eigenstratgeno row count
didn't exactly match the call files' SNP count. All of the checks below
are hard failures (sys.exit, non-zero status, nothing written to the
final output path) rather than warnings:
  - call files must all be the same length as each other
  - every character in every call file must be one of 0/1/2/9
  - new ancient sample IDs must not collide with each other or with an
    existing ID already in the panel's .ind file
  - the panel .eigenstratgeno row count must exactly equal the call
    file length (checked as an immediate hard stop the moment a row
    index would run past the call length -- not after the fact)
  - output is written to a temp file in the same directory as the
    final path and atomically renamed into place (os.replace) only
    after a full, verified-length pass succeeds -- a failed run never
    leaves a partial-but-plausible-looking file at the real output path

Usage:
  python3 merge_ancient_into_panel.py \\
    --panel-geno panel.eigenstratgeno --panel-ind panel.ind \\
    --calls SAMPLE1=sample1.calls.txt SAMPLE2=sample2.calls.txt ... \\
    --ancient-poplabel Ancient \\
    --out-geno merged.eigenstratgeno --out-ind merged.ind

Each calls entry's file must be pseudo_haploid_call.py's output run
against the SAME panel .snp file as --panel-geno/--panel-ind, i.e. same
SNP count and order -- this script checks call-file lengths match each
other and match the panel's own row count exactly, but (being unable to
see the .snp file here) cannot independently verify the calls were
actually generated against THIS panel's SNP order rather than some
other panel that happens to have the same SNP count. That verification
is check_ref.py's/pseudo_haploid_call.py's job, upstream of this script.
"""
import argparse
import os
import sys
import tempfile

VALID_CALL_CHARS = frozenset("0129")


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--panel-geno", required=True)
    p.add_argument("--panel-ind", required=True)
    p.add_argument(
        "--calls", nargs="+", required=True,
        help="SAMPLE_ID=path/to/calls.txt, one per ancient sample (pseudo_haploid_call.py output)",
    )
    p.add_argument(
        "--ancient-poplabel", default="Ancient",
        help="population label written to .ind for every ancient sample added here -- "
             "this label must NOT be listed in the poplistname file given to smartpca, "
             "so these samples get lsqproject'd instead of used to build eigenvectors",
    )
    p.add_argument("--out-geno", required=True)
    p.add_argument("--out-ind", required=True)
    return p.parse_args()


def read_existing_ind_ids(panel_ind_path):
    ids = []
    with open(panel_ind_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            ids.append(line.split()[0])
    return ids


def atomic_writer(final_path):
    """Open a temp file next to final_path; caller must call .commit() or
    the temp file is left behind (never silently replaces final_path)."""
    d = os.path.dirname(os.path.abspath(final_path)) or "."
    fd, tmp_path = tempfile.mkstemp(prefix=os.path.basename(final_path) + ".", suffix=".tmp", dir=d)
    fh = os.fdopen(fd, "w")

    class _Handle:
        def write(self_inner, s):
            fh.write(s)

        def commit(self_inner):
            fh.close()
            os.replace(tmp_path, final_path)

        def abort(self_inner):
            fh.close()
            os.remove(tmp_path)

    return _Handle()


def main():
    args = parse_args()

    samples, call_strings = [], []
    for spec in args.calls:
        if "=" not in spec:
            sys.exit(f"--calls entries must be SAMPLE_ID=path, got: {spec}")
        sample_id, path = spec.split("=", 1)
        with open(path) as f:
            s = f.read().replace("\n", "")
        samples.append(sample_id)
        call_strings.append(s)

    if len(set(samples)) != len(samples):
        dupes = sorted({s for s in samples if samples.count(s) > 1})
        sys.exit(f"duplicate sample IDs within --calls: {dupes}")

    existing_ids = set(read_existing_ind_ids(args.panel_ind))
    collisions = sorted(set(samples) & existing_ids)
    if collisions:
        sys.exit(
            f"ancient sample ID(s) already present in panel .ind, refusing to "
            f"add a duplicate individual: {collisions}"
        )

    n_snps = len(call_strings[0])
    for sid, s in zip(samples, call_strings):
        if len(s) != n_snps:
            sys.exit(
                f"call file length mismatch: {sid} has {len(s)} calls, "
                f"expected {n_snps} (from {samples[0]}) -- these must all "
                f"be pseudo_haploid_call.py runs against the SAME panel "
                f".snp file"
            )
        bad_chars = set(s) - VALID_CALL_CHARS
        if bad_chars:
            sys.exit(
                f"call file for {sid} contains invalid character(s) {sorted(bad_chars)} "
                f"-- only '0','1','2','9' are valid EIGENSTRAT genotype calls"
            )

    sys.stderr.write(
        f"[merge] {len(samples)} ancient samples ({', '.join(samples)}), "
        f"{n_snps} SNP calls expected per sample, validated against {len(existing_ids)} "
        f"existing panel individuals (no ID collisions)\n"
    )

    geno_out = atomic_writer(args.out_geno)
    n_lines = 0
    try:
        with open(args.panel_geno) as fin:
            for i, line in enumerate(fin):
                if i >= n_snps:
                    raise SystemExit(
                        f"panel .eigenstratgeno has MORE rows than the call files' SNP "
                        f"count ({n_snps}) -- hit row {i + 1} while still reading. This "
                        f"means the call files were not generated against this panel's "
                        f".snp file. Aborting, no output written to {args.out_geno}."
                    )
                line = line.rstrip("\n")
                extra = "".join(s[i] for s in call_strings)
                geno_out.write(line + extra + "\n")
                n_lines += 1
        if n_lines != n_snps:
            raise SystemExit(
                f"panel .eigenstratgeno has FEWER rows ({n_lines}) than the call files' "
                f"SNP count ({n_snps}) -- these must match exactly (both are supposed to "
                f"follow the same panel .snp order). Aborting, no output written to "
                f"{args.out_geno}."
            )
    except SystemExit:
        geno_out.abort()
        raise
    except BaseException:
        geno_out.abort()
        raise
    geno_out.commit()
    sys.stderr.write(f"[merge] OK: {n_lines} SNP rows, counts matched exactly\n")

    ind_out = atomic_writer(args.out_ind)
    try:
        with open(args.panel_ind) as fin:
            for line in fin:
                ind_out.write(line if line.endswith("\n") else line + "\n")
        for sid in samples:
            ind_out.write(f"{sid}\tU\t{args.ancient_poplabel}\n")
    except BaseException:
        ind_out.abort()
        raise
    ind_out.commit()

    sys.stderr.write(f"[merge] wrote {args.out_geno} and {args.out_ind}\n")


if __name__ == "__main__":
    main()
