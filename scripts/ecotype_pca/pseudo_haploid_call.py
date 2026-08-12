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
   branch ORYZA_BESTHIT_HANDOFF.md 7.5).
   OPERATIONAL NOTE (2026-08-12, see docs/ECOTYPE_PCA_EXECUTION_PLAN.md):
   the intended primary/sensitivity split is to run this script TWICE
   per sample per panel -- once with the default (TV = transversion-only,
   the primary result) and once with --no-transversions-only (ALL =
   every biallelic site called, a sensitivity check). "ALL" here still
   means "no damage-aware trimming is applied", not "damage-corrected" --
   this project has no read-trimming/damage-rescaling step yet, so ALL
   should be read as "less conservative", not as "safe because trimmed".
   If TV and ALL give the same population call for a sample, that
   agreement is the actual evidence of robustness; if a given site's
   functional interest depends on a transition (e.g. DROT1, a C/T SNP --
   see docs/ECOTYPE_PCA_EXECUTION_PLAN.md's ecotype-panel section), it
   can only ever appear in the ALL track, never in TV.

4. REF/ALT COLUMN CONVENTION MUST BE VERIFIED PER PANEL BEFORE TRUSTING
   REAL OUTPUT. This script follows the EIGENSOFT CONVERTF/README's
   literal definition: .snp file column 5 = reference allele, column 6 =
   variant/alt allele, and genotype 0 = zero copies of the reference
   (i.e. homozygous alt), 2 = two copies of reference (homozygous ref).

   !!! DO NOT decide --swap-ref-alt from check_ref.py's FASTA-match rate
   alone !!! (2026-08-12 correction, see docs/ECOTYPE_PCA_EXECUTION_PLAN.md
   P0-1 -- an earlier version of this note conflated two different
   questions). check_ref.py's "snp" mode answers "does .snp column 5 (or
   6) match the base actually present in irgsp.fa at that position" --
   that is a statement about how the panel's source data labeled REF/ALT
   relative to the genome. It is NOT the same question as "does this
   script's 0/2 encoding match how the panel's OWN .eigenstratgeno matrix
   already encodes its existing (modern) individuals at that site" --
   which is the thing that actually matters here, because
   merge_ancient_into_panel.py concatenates this script's output as new
   columns directly onto that same matrix. A panel whose source data used
   an unusual REF/ALT labeling relative to the genome can still be 100%
   internally self-consistent (this script's reading of column 5/6 can
   still agree with the modern genotype matrix) -- FASTA mismatch alone
   does not prove a real 0/2 encoding bug, and a high FASTA match rate
   does not prove the absence of one either.

   The authoritative check is a LEAVE-ONE-OUT SIMULATION using a modern
   sample already in the panel with a known population label (see
   docs/ECOTYPE_PCA_EXECUTION_PLAN.md Phase 0.5): mask that sample down to
   an ancient sample's covered-site pattern, run this script's same
   read-simulation/allele-matching logic against it, merge the result in
   as an extra column, and confirm it still lsqprojects back near its own
   known population. If the 0/2 encoding here disagreed with the panel's
   existing matrix, this simulated individual would systematically
   project AWAY from its true population (often toward whichever group is
   genotypically "opposite" at those sites), not merely with more noise --
   that is a stronger and more direct signal than any FASTA spot-check.
   Use --swap-ref-alt only once that simulation (or an equivalent direct
   comparison against the panel's own already-known genotypes for a
   sample it contains) has actually shown a flip is needed for a given
   panel -- not as a default response to check_ref.py's percentage.

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
                    help="disable default transition-SNP exclusion (see docstring point 3 -- "
                         "this is the 'ALL' sensitivity track, run in addition to the default "
                         "'TV' track, not instead of it)")
    p.add_argument("--swap-ref-alt", action="store_true",
                    help="invert which .snp column (5 vs 6) is treated as the reference allele "
                         "-- only pass this after a leave-one-out simulation has actually shown "
                         "it's needed for this panel, see docstring point 4")
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

    # NOTE (2026-08-12): n_no_coverage and n_allele_mismatch used to be
    # merged into a single n_uncovered counter, which made it impossible
    # to tell "this site simply had no reads" apart from "reads existed
    # but matched neither panel allele" (sequencing error / third allele /
    # residual damage the transversion filter didn't catch) -- these are
    # different diagnostic signals and are now reported separately. See
    # docs/ECOTYPE_PCA_EXECUTION_PLAN.md's reporting-field list.
    n_total = n_transition_skipped = n_no_allele_info = 0
    n_no_coverage = n_allele_mismatch = n_called = 0

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
                n_no_coverage += 1
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
                # Treat as missing rather than guessing. Distinct from
                # n_no_coverage: reads exist here, they just don't match
                # either panel allele.
                n_allele_mismatch += 1
                fout.write("9\n")

    n_uncovered = n_no_coverage + n_allele_mismatch
    sys.stderr.write(
        f"[pseudo_haploid_call] total={n_total} transition_skipped={n_transition_skipped} "
        f"no_allele_info={n_no_allele_info} no_coverage={n_no_coverage} "
        f"allele_mismatch={n_allele_mismatch} called={n_called} "
        f"missing={n_total - n_called}\n"
    )
    if args.report:
        with open(args.report, "w") as r:
            r.write("metric\tvalue\n")
            r.write(f"total_panel_snps\t{n_total}\n")
            r.write(f"transition_skipped\t{n_transition_skipped}\n")
            r.write(f"no_allele_info\t{n_no_allele_info}\n")
            r.write(f"no_coverage\t{n_no_coverage}\n")
            r.write(f"allele_mismatch\t{n_allele_mismatch}\n")
            r.write(f"uncovered_total\t{n_uncovered}\n")
            r.write(f"called\t{n_called}\n")
            r.write(f"missing\t{n_total - n_called}\n")


if __name__ == "__main__":
    main()
