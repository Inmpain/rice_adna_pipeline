#!/usr/bin/env bash
set -euo pipefail

# =====================================================================
# submit_map.sh
# sbatch 提交阶段②最终定量比对任务: 每样本一个任务(Bowtie2 比对 irgsp.fa)。
# 注意: 需先跑完 merge(合并), 即 01.merge/ 里已有 combined.fastq.gz。
#
# 用法: bash submit_map.sh [threads]
# =====================================================================

THREADS="${1:-20}"
PARTITION="${PARTITION:-comp}"
EXCLUDE="${EXCLUDE:-node05,node06}"

BASE=/home/scratch/yinmt202607/nanzuo
MERGE_DIR="$BASE/01.merge"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SLURM_LOG="$BASE/_logs/slurm"
mkdir -p "$SLURM_LOG"

for fq in "$MERGE_DIR"/*.combined.fastq.gz; do
  [[ -f "$fq" ]] || continue
  sample=$(basename "$fq" .combined.fastq.gz)
  out="$BASE/02.map_irgsp/${sample}.dedup.bam"
  [[ -s "$out" ]] && { echo "[$sample] 已存在, 跳过"; continue; }
  jobid=$(sbatch --parsable --job-name="nz_map_${sample}" \
    --cpus-per-task "$THREADS" --mem 24G --time 12:00:00 \
    --partition "$PARTITION" --exclude "$EXCLUDE" \
    --output "$SLURM_LOG/map_${sample}.%j.out" \
    "$SCRIPT_DIR/map_irgsp_single.sh" "$sample" "$THREADS")
  echo "[$sample] 已提交: $jobid"
done
