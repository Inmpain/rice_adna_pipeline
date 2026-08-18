#!/usr/bin/env bash
set -euo pipefail

# =====================================================================
# submit_map.sh
# sbatch 提交阶段②最终定量比对任务: 每样本 × 每个 mapper(bwa|bt2new) 一个任务。
# 注意: 需先跑完 merge(合并), 即 01.merge/ 里已有 combined.fastq.gz。
#
# 用法: bash submit_map.sh [threads]
# =====================================================================

THREADS="${1:-20}"
PARTITION="${PARTITION:-comp}"

BASE=/home/scratch/yinmt202607/nanzuo
MERGE_DIR="$BASE/01.merge"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SLURM_LOG="$BASE/_logs/slurm"
mkdir -p "$SLURM_LOG"

MAPPERS=(bwa bt2new)

for fq in "$MERGE_DIR"/*.combined.fastq.gz; do
  [[ -f "$fq" ]] || continue
  sample=$(basename "$fq" .combined.fastq.gz)
  for mapper in "${MAPPERS[@]}"; do
    out="$BASE/02.map_irgsp/$mapper/${sample}.dedup.bam"
    [[ -s "$out" ]] && { echo "[$sample][$mapper] 已存在, 跳过"; continue; }
    jobid=$(sbatch --parsable --job-name="nz_map_${mapper}_${sample}" \
      --cpus-per-task "$THREADS" --mem 24G --time 12:00:00 \
      --partition "$PARTITION" \
      --output "$SLURM_LOG/map_${mapper}_${sample}.%j.out" \
      "$SCRIPT_DIR/map_irgsp_single.sh" "$sample" "$mapper" "$THREADS")
    echo "[$sample][$mapper] 已提交: $jobid"
  done
done
