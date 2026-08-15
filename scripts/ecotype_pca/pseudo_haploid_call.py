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

   CORRECTNESS REQUIREMENT THIS RELIES ON (2026-08-15 fix, see point 5
   below): TV and ALL must make the IDENTICAL pseudo-haploid call at any
   site both tracks call, so that a genuine TV/ALL disagreement can only
   mean "this site's population signal depends on whether transitions are
   trusted", not "these two runs happened to draw different reads by
   chance". Before this fix that guarantee did not hold -- see point 5.

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

5. PER-SITE INDEPENDENT RNG (2026-08-15 fix, GPT review of the original
   version, actually reproduced and verified -- see
   scripts/ecotype_pca/tests/ in this same commit). The original version
   called `random.seed(args.seed)` once globally and then drew from that
   single shared stream via `random.choice(reads)` only at sites that
   reached that line -- so in ALL mode, every transition site consumed
   one extra draw from the stream that TV mode never took, shifting the
   RNG state at every transversion site that came after it in the panel's
   .snp file order. Confirmed with a synthetic BAM: at a shared
   transversion site, seed=1, TV called genotype 0 and ALL called
   genotype 2 for THE SAME SITE, THE SAME READS, THE SAME SEED -- purely
   because ALL had consumed one extra random.choice() call at an earlier
   transition site. This meant TV vs ALL disagreement was contaminated by
   run-order-dependent sampling noise and could not be safely attributed
   to "the population signal depends on transitions" alone.

   Fixed by deriving each site's random state ONLY from
   (--seed, contig, 1-based position) via a stable hash (hashlib, not the
   builtin hash() which is process-salted and not reproducible across
   runs), independent of what other sites were visited before it, in
   either track, in either run. This also happens to let the per-site
   pick be done as a single-pass reservoir sample during BAM indexing
   (point 6 below) instead of storing every covered read.

6. MEMORY: the coverage index used to store the FULL LIST of read bases
   at every covered genomic position across the WHOLE BAM (not just panel
   SNP positions) before this fix -- for a real low-input ancient sample
   this is bounded by total aligned bases, not panel size, but is still a
   real risk for any sample with non-trivial coverage or duplication.
   Fixed via reservoir sampling (Algorithm R, one candidate kept in
   memory per site at a time) combined directly into the per-site
   deterministic RNG from point 5, so the index now stores one base per
   covered position, not a list.

Usage:
  python3 pseudo_haploid_call.py \
    --bam <sample>.besthit_oryza.irgsp.bam \
    --panel-snp <panel>.snp \
    --out <sample>.<panel>.pseudohap.txt \
    --report <sample>.<panel>.pseudohap.report.tsv
"""
import argparse
import hashlib
import random
import sys
from pathlib import Path

import pysam

TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
VALID_BASES = {"A", "C", "G", "T"}


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
    p.add_argument("--contig-format", default="chr%02d",
                    help="printf-style format used to turn the panel .snp file's numeric "
                         "chromosome column into a BAM contig name (default 'chr%%02d' -> "
                         "'chr01'). Validated against the BAM's actual contig names at "
                         "startup -- see docstring point about contig validation.")
    return p.parse_args()


def is_transition(ref, alt):
    return (ref, alt) in TRANSITIONS


def site_seed(seed, contig, pos):
    """Stable (process- and run-independent) per-site seed from
    (--seed, contig, 1-based pos). hashlib, not the builtin hash(), because
    builtin hash() of strings is randomized per-process (PYTHONHASHSEED)
    for security and is NOT reproducible across separate script runs --
    using it here would silently reintroduce a different version of the
    same non-reproducibility problem this fix exists to eliminate."""
    digest = hashlib.sha256(f"{seed}:{contig}:{pos}".encode()).digest()
    return int.from_bytes(digest[:8], "big")


def build_coverage_index(bam_path, min_mapq, min_baseq, seed):
    """
    One pass over the whole BAM -> {(contig, 1-based pos): base}.
    Reservoir-samples (Algorithm R) exactly one base per covered position
    using that site's own deterministic RNG (site_seed) as candidates are
    encountered, so the index never holds more than one base per site
    (fixes point 6) and the chosen base depends only on (seed, contig,
    pos, the set of qualifying reads at that site) -- never on iteration
    order or on what other sites were visited (fixes point 5).
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
            pos = col.reference_pos + 1
            rng = None
            n_seen = 0
            chosen = None
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
                base = base.upper()
                n_seen += 1
                if rng is None:
                    rng = random.Random(site_seed(seed, contig, pos))
                if rng.randrange(n_seen) == 0:
                    chosen = base
            if chosen is not None:
                idx[(contig, pos)] = chosen
    bam.close()
    return idx


def validate_contig_naming(bam_references, contig_format, sample_chrom_ints, out=sys.stderr):
    """Fail loudly at startup if --contig-format doesn't actually match any
    real contig in the BAM, instead of silently producing all-missing
    output for every SNP (the failure mode being fixed here)."""
    bam_ref_set = set(bam_references)
    matched = 0
    for c in sample_chrom_ints:
        if (contig_format % c) in bam_ref_set:
            matched += 1
    if matched == 0 and sample_chrom_ints:
        out.write(
            f"FATAL: none of the panel's chromosome numbers, formatted with "
            f"--contig-format '{contig_format}', match any contig name in the BAM.\n"
            f"  tried (sample): {[contig_format % c for c in sorted(sample_chrom_ints)[:5]]}\n"
            f"  BAM references (sample): {sorted(bam_ref_set)[:10]}\n"
            f"Fix --contig-format to match the BAM's actual naming, do not proceed --\n"
            f"every SNP would silently come out as coverage-9 otherwise.\n"
        )
        sys.exit(2)


def main():
    args = parse_args()
    transversions_only = not args.no_transversions_only

    panel_path = Path(args.panel_snp)
    if not panel_path.is_file():
        sys.exit(f"FATAL: --panel-snp not found: {panel_path}")
    bam_path = Path(args.bam)
    if not bam_path.is_file():
        sys.exit(f"FATAL: --bam not found: {bam_path}")
    if not (bam_path.with_suffix(bam_path.suffix + ".bai").is_file()
            or Path(str(bam_path) + ".bai").is_file()
            or bam_path.with_suffix(".bai").is_file()):
        sys.exit(f"FATAL: no .bai index found next to {bam_path} -- run `samtools index` first "
                  f"(pysam.pileup requires an indexed BAM; without one this would either crash "
                  f"or silently misbehave depending on pysam version).")

    # --- pre-scan panel .snp for structural problems + contig-naming sanity,
    # before opening the (possibly large) BAM ---
    seen_ids = set()
    seen_pos = set()
    dup_ids = dup_pos = 0
    sample_chrom_ints = set()
    n_lines = 0
    with open(panel_path) as fh:
        for lineno, line in enumerate(fh, 1):
            n_lines += 1
            fields = line.split()
            if len(fields) < 4:
                sys.exit(f"FATAL: {panel_path}:{lineno}: expected >=4 columns "
                          f"(id chrom genpos pos[, ref, alt]), got {len(fields)}: {line!r}")
            snp_id, chrom_field, _genpos, pos_field = fields[0], fields[1], fields[2], fields[3]
            try:
                pos = int(pos_field)
            except ValueError:
                sys.exit(f"FATAL: {panel_path}:{lineno}: position column is not an integer: "
                          f"{pos_field!r}")
            if snp_id in seen_ids:
                dup_ids += 1
            seen_ids.add(snp_id)
            key = (chrom_field, pos)
            if key in seen_pos:
                dup_pos += 1
            seen_pos.add(key)
            try:
                sample_chrom_ints.add(int(chrom_field.lstrip("0") or "0"))
            except ValueError:
                pass  # unsupported-chromosome lines are handled (and counted) in the main loop
            if len(fields) >= 6:
                ref, alt = fields[4].upper(), fields[5].upper()
                if len(ref) != 1 or len(alt) != 1 or ref not in VALID_BASES \
                        or alt not in VALID_BASES or ref == alt:
                    sys.exit(f"FATAL: {panel_path}:{lineno}: ref/alt must be two distinct single "
                              f"bases from {{A,C,G,T}}, got ref={fields[4]!r} alt={fields[5]!r}. "
                              f"Do not proceed with an invalid panel -- fix the source .snp.")
    if n_lines == 0:
        sys.exit(f"FATAL: {panel_path} is empty")
    if dup_ids:
        sys.exit(f"FATAL: {panel_path} has {dup_ids} duplicate SNP IDs -- output would have "
                  f"ambiguous row identity when merged into the panel's genotype matrix, fix "
                  f"the panel .snp before calling against it")
    if dup_pos:
        sys.exit(f"FATAL: {panel_path} has {dup_pos} duplicate (chrom,pos) positions -- fix "
                  f"the panel .snp before calling against it")

    bam_probe = pysam.AlignmentFile(str(bam_path), "rb")
    bam_references = list(bam_probe.references)
    n_mapped = bam_probe.mapped
    bam_probe.close()
    if n_mapped == 0:
        sys.exit(f"FATAL: {bam_path} has zero mapped reads (empty BAM) -- refusing to produce "
                  f"an all-missing output silently; confirm this is the intended input")
    validate_contig_naming(bam_references, args.contig_format, sample_chrom_ints)

    sys.stderr.write("[pseudo_haploid_call] indexing BAM coverage...\n")
    cov = build_coverage_index(str(bam_path), args.min_mapq, args.min_baseq, args.seed)
    sys.stderr.write(f"[pseudo_haploid_call] {len(cov)} covered positions in BAM\n")

    n_total = n_transition_skipped = n_no_allele_info = n_unsupported_chromosome = 0
    n_no_coverage = n_allele_mismatch = n_called = 0

    with open(panel_path) as fin, open(args.out, "w") as fout:
        for line in fin:
            n_total += 1
            fields = line.split()
            # structural validity of this line was already confirmed in the
            # pre-scan above (len>=4, integer pos) -- safe to index directly now.
            chrom_field, pos = fields[1], int(fields[3])
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
                n_unsupported_chromosome += 1
                fout.write("9\n")
                continue
            contig = args.contig_format % chrom_int

            base = cov.get((contig, pos))
            if base is None:
                n_no_coverage += 1
                fout.write("9\n")
                continue

            if base == ref:
                fout.write("2\n")
                n_called += 1
            elif base == alt:
                fout.write("0\n")
                n_called += 1
            else:
                n_allele_mismatch += 1
                fout.write("9\n")

    n_uncovered = n_no_coverage + n_allele_mismatch
    callable_sites = n_total - n_transition_skipped - n_no_allele_info - n_unsupported_chromosome
    # Two different denominators answer two different questions -- conflating
    # them into one "call_rate" (2026-08-15 correction, GPT review of a4fb1e6)
    # made the number look artificially tiny and unusable, because it's
    # dominated by no_coverage (a coverage-depth fact, not a calling-quality
    # fact):
    #   eligible_site_call_rate: of every site this run was even ALLOWED to
    #     attempt (excludes transition-skipped/no-allele-info/unsupported-
    #     chromosome, which are excluded by design, not by data), what
    #     fraction actually got called -- this is coverage-dominated, most
    #     of the denominator is typically no_coverage for a sparse ancient
    #     sample, and that's expected, not a quality signal.
    #   allele_match_rate_among_covered: of sites where a read WAS actually
    #     drawn (called + allele_mismatch), what fraction matched a known
    #     panel allele -- this is the actual data-quality signal (sequencing
    #     error / third allele / residual damage the TV filter didn't catch
    #     would show up here, not in the coverage-dominated number above).
    eligible_site_call_rate = (n_called / callable_sites) if callable_sites else 0.0
    n_drawn = n_called + n_allele_mismatch
    allele_match_rate_among_covered = (n_called / n_drawn) if n_drawn else 0.0
    sys.stderr.write(
        f"[pseudo_haploid_call] total={n_total} transition_skipped={n_transition_skipped} "
        f"no_allele_info={n_no_allele_info} unsupported_chromosome={n_unsupported_chromosome} "
        f"no_coverage={n_no_coverage} allele_mismatch={n_allele_mismatch} called={n_called} "
        f"missing={n_total - n_called} callable_sites={callable_sites} "
        f"eligible_site_call_rate={eligible_site_call_rate:.6f} "
        f"allele_match_rate_among_covered={allele_match_rate_among_covered:.4f}\n"
    )
    assert n_called + n_no_coverage + n_allele_mismatch + n_transition_skipped \
        + n_no_allele_info + n_unsupported_chromosome == n_total, \
        "category counts don't sum to total_panel_snps -- this is a bug, every site must land " \
        "in exactly one category"

    if args.report:
        with open(args.report, "w") as r:
            r.write("metric\tvalue\n")
            r.write(f"total_panel_snps\t{n_total}\n")
            r.write(f"transition_skipped\t{n_transition_skipped}\n")
            r.write(f"no_allele_info\t{n_no_allele_info}\n")
            r.write(f"unsupported_chromosome\t{n_unsupported_chromosome}\n")
            r.write(f"no_coverage\t{n_no_coverage}\n")
            r.write(f"allele_mismatch\t{n_allele_mismatch}\n")
            r.write(f"uncovered_total\t{n_uncovered}\n")
            r.write(f"called\t{n_called}\n")
            r.write(f"missing\t{n_total - n_called}\n")
            r.write(f"callable_sites\t{callable_sites}\n")
            r.write(f"eligible_site_call_rate\t{eligible_site_call_rate:.6f}\n")
            r.write(f"allele_match_rate_among_covered\t{allele_match_rate_among_covered:.6f}\n")


if __name__ == "__main__":
    main()
