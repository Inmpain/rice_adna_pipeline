#!/usr/bin/env bash
set -euo pipefail

# =====================================================================
# 02a_extract_bt2_new_single.sh
# 单个capture panel文件的Bowtie2新参数(-N1)提取, 供sbatch按文件粒度调用。
#
# 用法: bash 02a_extract_bt2_new_single.sh <fq_path> <threads>
# =====================================================================

fq="$1"
THREADS="${2:-4}"

BASE=/home/scratch/yinmt202607/tests/param_matrix_bt2_vs_bwa
BT2_IDX=/home/scratch/yinmt202607/db/asian_rice_panel_index/asian_rice_panel.fa
OUT_DIR="$BASE/00.extraction/bt2_new"

BT2_NEW_EXTRA=( -k 3 -L 22 -N 1 -i S,1,1.15 --mp 1,1 --rdg 0,1 --rfg 0,1 --score-min L,0,-0.1 --no-unal )

sample=$(basename "$fq" .bbduk.lowcomp_filtered.fq)
bam="$OUT_DIR/bam/${sample}.bt2new.bam"
fastq_out="$OUT_DIR/fastq/${sample}.bt2new.fastq.gz"

if [[ -s "$fastq_out" ]]; then
    echo "[$sample] 跳过(已存在)"
    exit 0
fi

echo "[$sample] 开始提取..."
bowtie2 -p "$THREADS" "${BT2_NEW_EXTRA[@]}" -x "$BT2_IDX" -U "$fq" \
    2> "$BASE/logs/${sample}.bt2new_extract.log" \
    | samtools view -@ "$THREADS" -bh -F 0x904 -o "$bam" -
samtools fastq -@ "$THREADS" "$bam" | gzip -c > "$fastq_out"

echo "[$sample] 完成"
