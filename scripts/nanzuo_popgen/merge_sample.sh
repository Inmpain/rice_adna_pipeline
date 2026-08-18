#!/usr/bin/env bash
set -euo pipefail

# =====================================================================
# merge_sample.sh
# 合并单样本三份候选 FASTQ(popgen 提取 + shotgun 提取 + function bam→fq)。
# 登录节点直接跑, 不 sbatch。三个来源都是 .gz, 直接 cat 拼接(gzip 流可串联)。
#
# 用法: bash merge_sample.sh <YWL1-AXXXX>
#   输出: nanzuo/01.merge/{sample}.combined.fastq.gz
# =====================================================================

sample="$1"
[[ -n "$sample" ]] || { echo "用法: $0 <YWL1-AXXXX>" >&2; exit 1; }

BASE=/home/scratch/yinmt202607/nanzuo
OUT_DIR="$BASE/01.merge"
LOG_DIR="$BASE/_logs/merge"

sources=(
  "$BASE/00.extract/popgen/${sample}.bt2.primary_mapped.fastq.gz"
  "$BASE/00.extract/shotgun/${sample}.bt2.primary_mapped.fastq.gz"
  "$BASE/00.extract/function/${sample}.bam2fq.fastq.gz"
)

out_fq="$OUT_DIR/${sample}.combined.fastq.gz"
mkdir -p "$OUT_DIR" "$LOG_DIR"

present=()
for s in "${sources[@]}"; do
  if [[ -s "$s" ]]; then
    present+=("$s")
  else
    echo "[$sample] WARNING: 缺失来源 $s" >&2
  fi
done

[[ ${#present[@]} -gt 0 ]] || { echo "[$sample] ERROR: 三份来源都缺失" >&2; exit 1; }

tmp_fq="$out_fq.tmp"
rm -f "$tmp_fq"
trap 'rm -f "$tmp_fq"' EXIT

cat "${present[@]}" > "$tmp_fq" 2> "$LOG_DIR/${sample}.log"
mv "$tmp_fq" "$out_fq"
echo "[$sample] 合并完成(来源数: ${#present[@]}/3)"
