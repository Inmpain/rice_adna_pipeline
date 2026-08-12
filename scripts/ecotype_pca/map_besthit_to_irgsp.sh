#!/bin/bash
set -euo pipefail

# Maps besthit-filtered ancient Oryza reads to irgsp.fa -- the shared
# reference coordinate system all three ecotype-pca panels (29M_3k,
# 6.7M_720, Civan) are confirmed to use (see ECOTYPE_PCA_PANEL.md 3.1).
# One BAM per sample feeds pseudo-haploid genotype calling against any
# of the three panels -- mapping is only done once per sample, not once
# per panel.
#
# Input: <sample>.besthit_oryza.fastq.gz (codex/oryza-competitive-mapping
#   branch's ORYZA_BESTHIT_HANDOFF.md section 5.3 output). This FASTQ has
#   NOT been deduplicated yet (explicitly deferred "further downstream"
#   per that doc's section 5.2) -- deduplication happens here, since this
#   is that downstream step.
#
# 2026-08-12 hardening (see docs/ECOTYPE_PCA_EXECUTION_PLAN.md P0-2): this
# revision aligns the mapping/filtering/dedup chain with the project's
# main pipeline (scripts/server_originals/mapping.sh) instead of a
# from-literature reconstruction that had drifted from it:
#   - bwa aln -l 1024 -n 0.01 -o 2 is UNCHANGED -- this already matched
#     mapping.sh exactly, confirmed by direct diff, so the previous
#     "not verified against main pipeline" caveat is removed.
#   - samtools view -bh -F 0x904 is now applied right after samse, same
#     position as mapping.sh -- drops unmapped/secondary/supplementary
#     alignments before sorting, not left for a later ad-hoc read-count
#     fix. (bwa samse without -a does not itself emit secondary/
#     supplementary records, so in practice this mainly drops unmapped
#     reads here, but the filter is applied at the same pipeline stage
#     as the main pipeline for consistency and future-proofing.)
#   - dedup now goes through the same collate -> fixmate -m -> sort ->
#     markdup chain as mapping.sh, instead of running markdup directly
#     on a coordinate-sorted BAM with fixmate left as an "add it if it
#     errors" fallback. fixmate needs mate-score info that a plain
#     coordinate sort does not guarantee is present/correct for markdup
#     to make optimal duplicate-of-pair decisions; for single-end data
#     this matters less than for paired-end, but there is no reason to
#     run the riskier path when the safe one costs one extra sort pass.
#   - DELIBERATE remaining difference from mapping.sh: markdup here does
#     NOT use -r (does not remove duplicates), only flags them. This is
#     intentional, not an oversight -- pseudo_haploid_call.py already
#     filters flagged duplicates at pileup time (aln.is_duplicate check),
#     and the Phase 0 IRGSP coverage census (see execution plan) wants
#     duplicate-rate visible in the BAM as a QC signal, not silently
#     removed upstream. If a future consumer needs duplicates physically
#     removed, add -r at the call site rather than changing this default.
#   - the per-sample summary line now reports mapped-primary reads via
#     samtools view -c -F 4 (excludes unmapped; matches what "mapped
#     reads" should mean) instead of a bare `samtools view -c` that
#     could include unmapped records depending on BAM state, plus
#     separate MAPQ>=30 (primary SNP-calling threshold, matches main
#     pipeline's q30 filter step) and MAPQ>=20 (sensitivity-analysis
#     threshold, see ECOTYPE_PCA_PANEL.md / execution plan) counts.
#
# Usage: bash map_besthit_to_irgsp.sh <besthit_fastq_dir> <irgsp_fa> <out_dir> <sample1> [sample2 ...]

BESTHIT_DIR="${1:?usage: map_besthit_to_irgsp.sh <besthit_fastq_dir> <irgsp_fa> <out_dir> <sample1> [sample2 ...]}"
IRGSP_FA="${2:?see usage above}"
OUT_DIR="${3:?see usage above}"
shift 3

mkdir -p "$OUT_DIR"

for sample in "$@"; do
  fq="$BESTHIT_DIR/${sample}.besthit_oryza.fastq.gz"
  if [ ! -f "$fq" ]; then
    echo "[skip] $sample: $fq not found"
    continue
  fi

  echo "[map] $sample: bwa aln + samse + filter (-F 0x904) + sort"
  bwa aln -l 1024 -n 0.01 -o 2 -t 4 "$IRGSP_FA" "$fq" > "$OUT_DIR/${sample}.sai"
  bwa samse "$IRGSP_FA" "$OUT_DIR/${sample}.sai" "$fq" \
    | samtools view -@ 4 -bh -F 0x904 - \
    | samtools sort -@ 4 -o "$OUT_DIR/${sample}.sorted.bam" -
  rm -f "$OUT_DIR/${sample}.sai"

  echo "[dedup] $sample: collate | fixmate -m | sort | markdup (flag only, not -r)"
  samtools collate -@ 4 -O "$OUT_DIR/${sample}.sorted.bam" \
    | samtools fixmate -@ 4 -m - - \
    | samtools sort -@ 4 -o "$OUT_DIR/${sample}.fixmate_sorted.bam" -
  samtools markdup -@ 4 "$OUT_DIR/${sample}.fixmate_sorted.bam" "$OUT_DIR/${sample}.besthit_oryza.irgsp.bam"
  samtools index "$OUT_DIR/${sample}.besthit_oryza.irgsp.bam"
  rm -f "$OUT_DIR/${sample}.sorted.bam" "$OUT_DIR/${sample}.fixmate_sorted.bam"

  n_mapped=$(samtools view -@ 4 -c -F 4 "$OUT_DIR/${sample}.besthit_oryza.irgsp.bam")
  n_dup=$(samtools view -@ 4 -c -f 1024 "$OUT_DIR/${sample}.besthit_oryza.irgsp.bam")
  n_q30=$(samtools view -@ 4 -c -F 1028 -q 30 "$OUT_DIR/${sample}.besthit_oryza.irgsp.bam")
  n_q20=$(samtools view -@ 4 -c -F 1028 -q 20 "$OUT_DIR/${sample}.besthit_oryza.irgsp.bam")
  echo "[done] $sample: mapped=$n_mapped duplicates_flagged=$n_dup mapq>=30_nondup=$n_q30 mapq>=20_nondup=$n_q20 (duplicates not removed -- pseudo_haploid_call.py filters them at pileup time)"
done
