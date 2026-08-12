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
# aDNA-appropriate bwa aln settings (-l 1024 disables seeding since it's
# longer than any read, -n 0.01 relaxes the edit-distance threshold to
# tolerate damage-derived mismatches -- standard literature settings for
# ancient DNA mapping, e.g. Schubert et al. 2012). NOTE: these are NOT
# verified against this project's own main-pipeline BWA parameters in
# scripts/server_originals/ -- if you want exact consistency with that
# pipeline instead of literature defaults, check those scripts first.
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

  echo "[map] $sample: bwa aln + samse"
  bwa aln -l 1024 -n 0.01 -o 2 -t 4 "$IRGSP_FA" "$fq" > "$OUT_DIR/${sample}.sai"
  bwa samse "$IRGSP_FA" "$OUT_DIR/${sample}.sai" "$fq" \
    | samtools sort -@ 4 -o "$OUT_DIR/${sample}.sorted.bam" -
  rm -f "$OUT_DIR/${sample}.sai"

  echo "[dedup] $sample: samtools markdup"
  echo "  NOTE: running markdup directly on coordinate-sorted single-end"
  echo "  data without a prior fixmate pass -- this is the common SE"
  echo "  shortcut, but if samtools errors/warns about missing mate"
  echo "  score info, insert 'samtools collate -O | samtools fixmate -m'"
  echo "  before the sort step above and re-run."
  samtools markdup -@ 4 "$OUT_DIR/${sample}.sorted.bam" "$OUT_DIR/${sample}.besthit_oryza.irgsp.bam"
  samtools index "$OUT_DIR/${sample}.besthit_oryza.irgsp.bam"
  rm -f "$OUT_DIR/${sample}.sorted.bam"

  n_reads=$(samtools view -c "$OUT_DIR/${sample}.besthit_oryza.irgsp.bam")
  n_dup=$(samtools view -c -f 1024 "$OUT_DIR/${sample}.besthit_oryza.irgsp.bam")
  echo "[done] $sample: $n_reads mapped reads, $n_dup flagged as duplicates (not removed, just flagged -- pseudo_haploid_call.py filters them at pileup time)"
done
