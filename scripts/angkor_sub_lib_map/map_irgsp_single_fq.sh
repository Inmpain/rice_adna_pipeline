#!/usr/bin/env bash
# =====================================================================
# map_irgsp_single_fq.sh — 单份 fastq → irgsp dedup BAM（bowtie2 -N1）
#
# 与 nanzuo_popgen/map_irgsp_single.sh 同参数（bowtie2 -N1 + markdup
# 只标记 0x400 不删除）。输入任意 fastq/.fq.gz（单端），输出
# OUT_DIR/<out_stem>.dedup.bam(.bai)。用于 angkor 16 个 no-besthit 样
# 的三子库（shotgun / capture1 / capture2）单独 mapping，供 4 点投影。
#
# 用法: bash map_irgsp_single_fq.sh <input_fq> <out_stem> [threads]
#   OUT_DIR 由环境变量 SUB_LIB_BAM_DIR 控制（默认 angkor_sub_lib_bams_16）
# =====================================================================
set -euo pipefail

fq="$1"; OUT_STEM="$2"; THREADS="${3:-8}"
OUT_DIR="${SUB_LIB_BAM_DIR:-/home/scratch/yinmt202607/angkor_sub_lib_bams_16}"
IRGSP_BT2_IDX=/home/scratch/yinmt202607/db/asian_rice_panel_index/irgsp_bt2idx

[[ -s "$fq" ]] || { echo "FATAL: no input fastq: $fq" >&2; exit 1; }
mkdir -p "$OUT_DIR" "$OUT_DIR/logs"
out_bam="$OUT_DIR/${OUT_STEM}.dedup.bam"
[[ -s "$out_bam" ]] && { echo "[$OUT_STEM] exists, skip"; exit 0; }

# 0x004 unmapped + 0x100 secondary + 0x800 supplementary
EXCLUDE_FLAGS=0x904
# 与 nanzuo_popgen 阶段②同参数（bt2_new / -N1）
BT2_PARAMS=(-k 3 -L 22 -N 1 -i S,1,1.15 --mp 1,1 --rdg 0,1 --rfg 0,1 --score-min L,0,-0.1 --no-unal)

work="$(mktemp -d "$OUT_DIR/.${OUT_STEM}.XXXXXX")"
trap 'rm -rf "$work"' EXIT
mapped="$work/mapped.bam"

bowtie2 -p "$THREADS" "${BT2_PARAMS[@]}" -x "$IRGSP_BT2_IDX" -U "$fq" \
  2> "$OUT_DIR/logs/${OUT_STEM}.aln.log" \
  | samtools view -@ "$THREADS" -bh -F "$EXCLUDE_FLAGS" - \
  | samtools sort -@ "$THREADS" -o "$mapped" - 2>> "$OUT_DIR/logs/${OUT_STEM}.aln.log"

# markdup 只标记不删除(无 -r)：collate + fixmate -m + sort + markdup
samtools collate -@ "$THREADS" -O "$mapped" \
  | samtools fixmate -@ "$THREADS" -m - - \
  | samtools sort -@ "$THREADS" -o "$work/sorted.bam" -
samtools markdup -@ "$THREADS" "$work/sorted.bam" "$out_bam"
samtools index -@ "$THREADS" "$out_bam"
echo "OK $OUT_STEM -> $out_bam"
