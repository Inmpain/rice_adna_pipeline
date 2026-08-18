#!/usr/bin/env bash
set -euo pipefail

# =====================================================================
# extract_bwa_single.sh
# 阶段① 提取: 单个 FASTQ 用 BWA 比对 asian_rice_panel.fa, 只留 primary mapped
# reads, 转回 FASTQ。popgen / shotgun 两个来源共用本脚本(单文件粒度, 供 sbatch)。
#
# 用法: bash extract_bwa_single.sh <fq_path> <source> [threads]
#   source: popgen | shotgun
#   输出:  nanzuo/00.extract/{source}/{YWL1-AXXXX}.bwa.primary_mapped.fastq.gz
#
# 说明: 提取阶段的 BAM 只作中间产物(用完即删), 最终只保留 FASTQ;
#       最终定量比对在阶段②重新比对 irgsp.fa, 无需保留这里的 BAM。
# =====================================================================

fq="$1"
source="$2"
THREADS="${3:-20}"

# 工具加载: 优先 module, 失败则假定已在 PATH
module load bwa/ 2>/dev/null || module load bwa 2>/dev/null || true
module load samtools/ 2>/dev/null || module load samtools 2>/dev/null || true

PANEL_REF=/home/scratch/yinmt202607/db/asian_rice_panel_index/asian_rice_panel.fa
BASE=/home/scratch/yinmt202607/nanzuo
OUT_DIR="$BASE/00.extract/$source"
LOG_DIR="$BASE/_logs/extract_bwa"

# 与主线 / 9格矩阵测试一致
BWA_ALN_PARAMS=(-l 1024 -n 0.01 -o 2)
# 0x004 unmapped + 0x100 secondary + 0x800 supplementary
EXCLUDE_FLAGS=0x904

[[ -n "$fq" && -n "$source" ]] || { echo "用法: $0 <fq_path> <popgen|shotgun> [threads]" >&2; exit 1; }

sample=$(basename "$fq" | grep -oE 'YWL1[-_]A[0-9]+' | head -1 | tr '_' '-')
[[ -n "$sample" ]] || { echo "ERROR: 无法从 $fq 解析样本名" >&2; exit 1; }

out_fq="$OUT_DIR/${sample}.bwa.primary_mapped.fastq.gz"
if [[ -s "$out_fq" ]]; then
    echo "[$sample][$source] 已存在, 跳过"
    exit 0
fi

mkdir -p "$OUT_DIR" "$LOG_DIR"
tmp_bam="$(mktemp "$OUT_DIR/.${sample}.XXXXXX.bam")"
tmp_fq="$out_fq.tmp"
rm -f "$tmp_fq"
trap 'rm -f "$tmp_bam" "$tmp_fq"' EXIT

bwa aln "${BWA_ALN_PARAMS[@]}" -t "$THREADS" "$PANEL_REF" "$fq" 2> "$LOG_DIR/${sample}.aln.log" \
  | bwa samse "$PANEL_REF" - "$fq" 2>> "$LOG_DIR/${sample}.aln.log" \
  | samtools view -@ "$THREADS" -bh -F "$EXCLUDE_FLAGS" -o "$tmp_bam" - 2>> "$LOG_DIR/${sample}.aln.log"

samtools fastq -@ "$THREADS" -n "$tmp_bam" | gzip -c > "$tmp_fq"
mv "$tmp_fq" "$out_fq"

mapped=$(samtools view -c "$tmp_bam")
echo "[$sample][$source] 完成: $mapped primary mapped reads"
