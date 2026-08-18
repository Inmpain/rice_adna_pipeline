#!/usr/bin/env bash
set -euo pipefail

# =====================================================================
# bam_to_fastq_single.sh
# function 来源: 把已比对好的 BAM 转回 FASTQ(不跑 bwa 提取), 只取 primary mapped。
#
# 用法: bash bam_to_fastq_single.sh <bam_path> [threads]
#   输出: nanzuo/00.extract/function/{YWL1-AXXXX}.bam2fq.fastq.gz
# =====================================================================

bam="$1"
THREADS="${2:-4}"

module load samtools/ 2>/dev/null || module load samtools 2>/dev/null || true

BASE=/home/scratch/yinmt202607/nanzuo
OUT_DIR="$BASE/00.extract/function"
LOG_DIR="$BASE/_logs/bam_to_fastq"

# 0x004 unmapped + 0x100 secondary + 0x800 supplementary
EXCLUDE_FLAGS=0x904

[[ -n "$bam" ]] || { echo "用法: $0 <bam_path> [threads]" >&2; exit 1; }

sample=$(basename "$bam" | grep -oE 'YWL1[-_]A[0-9]+' | head -1 | tr '_' '-')
[[ -n "$sample" ]] || { echo "ERROR: 无法从 $bam 解析样本名" >&2; exit 1; }

out_fq="$OUT_DIR/${sample}.bam2fq.fastq.gz"
if [[ -s "$out_fq" ]]; then
    echo "[$sample][function] 已存在, 跳过"
    exit 0
fi

mkdir -p "$OUT_DIR" "$LOG_DIR"
tmp_fq="$out_fq.tmp"
rm -f "$tmp_fq"
trap 'rm -f "$tmp_fq"' EXIT

samtools fastq -@ "$THREADS" -n -F "$EXCLUDE_FLAGS" "$bam" 2> "$LOG_DIR/${sample}.log" \
  | gzip -c > "$tmp_fq"
mv "$tmp_fq" "$out_fq"
echo "[$sample][function] 完成"
