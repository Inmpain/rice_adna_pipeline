#!/usr/bin/env python3
"""
Pseudo-haploid genotype calling for one ancient sample against one PCA
panel's SNP list, for smartpca -lsqproject projection (see
docs/ECOTYPE_PCA_PANEL.md section 2 / 3.2 item 3 / 5.2).

DESIGN NOTES (read before trusting output):

1. Sparse-by-design: output has exactly one line per SNP in the panel's
   .snp file, IN THE SAME ORDER. Sites the sample's reads don't cover
   are "9" (missing). This is required so the output can be appended as
   one extra individual column onto the panel's own .eigenstratgeno file
   (same SNP set, same row order) -- see merge_ancient_into_panel.py
   (next script). Do NOT reorder or filter rows out; smartpca requires
   every individual's genotype file to have exactly the same SNP rows.

2. Pseudo-haploid = one read drawn at random per covered site, its base
   becomes the call. This is always homozygous by construction (0 or 2,
   never 1/het) -- standard practice for low-coverage ancient DNA, not a
   simplification specific to this project.

3. Transition SNPs (A/G, C/T pairs) are EXCLUDED from calling by default
   (emitted as missing, not dropped from the panel) via
   --transversions-only, which defaults to ON. This is deliberate: aDNA
   terminal damage produces C->T (5' end) and G->A (3' end) miscalls,
   and restricting genotype calls to transversion sites is a standard,
   widely-used way to sidestep this entirely without needing a
   damage-curve-specific read-trimming parameter (this project's own
   damage window calibration is still an open question -- see besthit
   branch ORYZA_BESTHIT_HANDOFF.md 7.5). Pass --no-transversions-only
   only if you have a specific reason to trust the damage profile at
   transition sites for a given sample/panel.

4. !!! REF/ALT COLUMN CONVENTION MUST BE VERIFIED PER PANEL BEFORE
    TRUSTING REAL OUTPUT !!! This script follows the EIGENSOFT
    CONVERTF/README's literal definition: .snp file column 5 = reference
    allele, column 6 = variant/alt allele, and genotype 0 = zero copies
    of the reference (i.e. homozygous alt), 2 = two copies of reference
    (homozygous ref). ECOTYPE_PCA_PANEL.md section 3.1 already found
    that this convention was NOT perfectly clean for one of our three
    panels (6.7M_720 matched the expected direction in only 183/200
    spot-checked sites, 91.5%, not 100%) -- so before running this
    script for real on a given panel, re-run check_ref.py in its
    existing "snp" mode against that panel's .snp file and confirm which
    column is actually acting as reference for that specific panel. If
    a panel's convention is flipped, use --swap-ref-alt to invert the
    0/2 assignment for that run rather than silently trusting the
    column order.

Usage:
  python3 pseudo_haploid_call.py \
    --bam <sample>.besthit_oryza.irgsp.bam \
    --panel-snp <panel>.snp \
    --out <sample>.<panel>.pseudohap.txt \
    --report <sample>.<panel>.pseudohap.report.tsv
"""
import argparse
import random
import sys

import pysam

TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--bam", required=True, help="sample BAM, mapped+markdup'd against irgsp.fa (see map_besthit_to_irgsp.sh)")
    p.add_argument("--panel-snp", required=True, help="panel .snp file (EIGENSTRAT format)")
    p.add_argument("--out", required=True, help="output: one 0/1/2/9 call per line, matching panel SNP order")
    p.add_argument("--report", help="optional: write call-rate summary here")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--min-mapq", type=int, default=20)
    p.add_argument("--min-baseq", type=int, default=20)
    p.add_argument("--no-transversions-only", action="store_true",
                    help="disable default transition-SNP exclusion (see docstring point 3 -- not recommended without a reason)")
    p.add_argument("--swap-ref-alt", action="store_true",
                    help="invert which .snp column (5 vs 6) is treated as the reference allele (see docstring point 4)")
    return p.parse_args()


def is_transition(ref, alt):
    return (ref, alt) in TRANSITIONS


def build_coverage_index(bam_path, min_mapq, min_baseq):
    """
    One pass over the whole BAM -> {(contig, 1-based pos): [(base, qual), ...]}.
    Ancient-sample coverage against a panel is sparse by design (see
    ECOTYPE_PCA_PANEL.md section 2), so this dict stays small -- bounded
    by total aligned bases in the sample's besthit-filtered read set, not
    by the panel size. This lets the per-SNP lookup below be O(1) instead
    of doing a separate pileup query for every one of a panel's (up to
    29 million) SNP positions.
    """
    idx = {}
    bam = pysam.AlignmentFile(bam_path, "rb")
    for contig in bam.references:
        for col in bam.pileup(
            contig,
            min_mapping_quality=min_mapq,
            stepper="samtools",
            ignore_overlaps=False,
            ignore_orphans=False,
        ):
            reads = []
            for pread in col.pileups:
                if pread.is_del or pread.is_refskip:
                    continue
                aln = pread.alignment
                if aln.is_duplicate or aln.is_secondary or aln.is_supplementary:
                    continue
                if aln.mapping_quality < min_mapq:
                    continue
                qpos = pread.query_position
                if qpos is None:
                    continue
                base = aln.query_sequence[qpos]
                qual = aln.query_qualities[qpos]
                if qual is None or qual < min_baseq:
                    continue
                reads.append(base.upper())
            if reads:
                idx[(contig, col.reference_pos + 1)] = reads
    bam.close()
    return idx


def main():
    args = parse_args()
    random.seed(args.seed)
    transversions_only = not args.no_transversions_only

    sys.stderr.write("[pseudo_haploid_call] indexing BAM coverage...\n")
    cov = build_coverage_index(args.bam, args.min_mapq, args.min_baseq)
    sys.stderr.write(f"[pseudo_haploid_call] {len(cov)} covered positions in BAM\n")

    n_total = n_transition_skipped = n_no_allele_info = n_uncovered = n_called = 0

    with open(args.panel_snp) as fin, open(args.out, "w") as fout:
        for line in fin:
            n_total += 1
            fields = line.split()
            chrom_field = fields[1]
            pos = int(fields[3])
            if len(fields) < 6:
                n_no_allele_info += 1
                fout.write("9\n")
                continue
            ref, alt = fields[4].upper(), fields[5].upper()
            if args.swap_ref_alt:
                ref, alt = alt, ref

            if transversions_only and is_transition(ref, alt):
                n_transition_skipped += 1
                fout.write("9\n")
                continue

            try:
                chrom_int = int(chrom_field.lstrip("0") or "0")
            except ValueError:
                fout.write("9\n")
                continue
            contig = "chr%02d" % chrom_int

            reads = cov.get((contig, pos))
            if not reads:
                n_uncovered += 1
                fout.write("9\n")
                continue

            base = random.choice(reads)
            if base == ref:
                fout.write("2\n")
                n_called += 1
            elif base == alt:
                fout.write("0\n")
                n_called += 1
            else:
                # neither known allele -- sequencing error / third allele /
                # residual damage the transversion filter didn't catch.
                # Treat as missing rather than guessing.
                fout.write("9\n")
                n_uncovered += 1

    sys.stderr.write(
        f"[pseudo_haploid_call] total={n_total} transition_skipped={n_transition_skipped} "
        f"no_allele_info={n_no_allele_info} called={n_called} "
        f"missing={n_total - n_called}\n"
    )
    if args.report:
        with open(args.report, "w") as r:
            r.write("metric\tvalue\n")
            r.write(f"total_panel_snps\t{n_total}\n")
            r.write(f"transition_skipped\t{n_transition_skipped}\n")
            r.write(f"no_allele_info\t{n_no_allele_info}\n")
            r.write(f"called\t{n_called}\n")
            r.write(f"missing\t{n_total - n_called}\n")


if __name__ == "__main__":
    main()
