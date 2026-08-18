#!/usr/bin/env bash
set -euo pipefail

# =====================================================================
# submit_extract.sh
# sbatch 提交阶段①全部提取任务: popgen + shotgun(bowtie2 提取) + function(bam→fq),
# 每文件一个任务。
#
# 用法: bash submit_extract.sh [threads]
# =====================================================================

THREADS="${1:-20}"
PARTITION="${PARTITION:-comp}"
EXCLUDE="${EXCLUDE:-node05,node06}"

BASE=/home/scratch/yinmt202607/nanzuo
POPGEN_DIR=/home/scratch/yinmt202607/2.nanzuo_popgen_yancheng
SHOTGUN_DIR=/home/scratch/yinmt202607/1.nanzuo_shotgun
FUNCTION_DIR=/home/scratch/yinmt202607/6.nanzuo_function_yancheng
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SLURM_LOG="$BASE/_logs/slurm"
mkdir -p "$SLURM_LOG"

extract_one() {
  local fq=$1 src=$2
  local sample out
  sample=$(basename "$fq" | grep -oE 'YWL1[-_]A[0-9]+' | head -1 | tr '_' '-')
  [[ -n "$sample" ]] || { echo "[$(basename "$fq")][$src] 跳过(非 YWL1 样本, 如阴性对照)"; return; }
  out="$BASE/00.extract/$src/${sample}.bt2.primary_mapped.fastq.gz"
  [[ -s "$out" ]] && { echo "[$sample][$src] 已存在, 跳过"; return; }
  local jobid
  jobid=$(sbatch --parsable --job-name="nz_ext_${sample}" \
    --cpus-per-task "$THREADS" --mem 32G --time 09:00:00 \
    --partition "$PARTITION" --exclude "$EXCLUDE" \
    --output "$SLURM_LOG/ext_${src}_${sample}.%j.out" \
    "$SCRIPT_DIR/extract_bt2_single.sh" "$fq" "$src" "$THREADS")
  echo "[$sample][$src] 已提交: $jobid"
}

bam2fq_one() {
  local bam=$1
  local sample out
  sample=$(basename "$bam" | grep -oE 'YWL1[-_]A[0-9]+' | head -1 | tr '_' '-')
  [[ -n "$sample" ]] || { echo "[$(basename "$bam")][function] 跳过(非 YWL1 样本)"; return; }
  out="$BASE/00.extract/function/${sample}.bam2fq.fastq.gz"
  [[ -s "$out" ]] && { echo "[$sample][function] 已存在, 跳过"; return; }
  local jobid
  jobid=$(sbatch --parsable --job-name="nz_bam2fq_${sample}" \
    --cpus-per-task 4 --mem 8G --time 02:00:00 \
    --partition "$PARTITION" --exclude "$EXCLUDE" \
    --output "$SLURM_LOG/bam2fq_${sample}.%j.out" \
    "$SCRIPT_DIR/bam_to_fastq_single.sh" "$bam" 4)
  echo "[$sample][function] 已提交: $jobid"
}

for fq in "$POPGEN_DIR"/*.bbduk.lowcomp_filtered.fq; do
  [[ -f "$fq" ]] || continue
  extract_one "$fq" popgen
done

for fq in "$SHOTGUN_DIR"/*.taxa_cleaned.fq.gz; do
  [[ -f "$fq" ]] || continue
  extract_one "$fq" shotgun
done

for bam in "$FUNCTION_DIR"/*.bam; do
  [[ -f "$bam" ]] || continue
  bam2fq_one "$bam"
done
