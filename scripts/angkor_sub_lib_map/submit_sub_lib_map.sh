#!/usr/bin/env bash
# =====================================================================
# submit_sub_lib_map.sh — 提交 16 个 no-besthit 样的三子库单独 irgsp mapping
#
# 16 样 × 3 子库（shotgun / capture1 / capture2）= 48 个 BAM，
# 每任务调 map_irgsp_single_fq.sh。reads 输入：
#   shotgun  = 5.angkor_shotgun_finished/.../{shotgun_lib}.prefiltered.IRGSP1.mapped.fq
#   capture1 = tests/param_matrix_bt2_vs_bwa/00.extraction/bt2_new/fastq/{robot}_RicePanel1.bt2new.fastq.gz
#   capture2 = 同上 RicePanel2
# 输出命名 {robot}_SG / {robot}_C1 / {robot}_C2（与 376 无冲突，可作投影样品 ID）。
#
# 用法: bash submit_sub_lib_map.sh [threads]
#   PARTITION / EXCLUDE 环境变量可覆盖（默认 comp / node05,node06）
# =====================================================================
set -euo pipefail

THREADS="${1:-8}"
PARTITION="${PARTITION:-comp}"
EXCLUDE="${EXCLUDE:-node05,node06}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
META="$SCRIPT_DIR/16_nobesthit_robot_lib.tsv"

SHOTGUN_DIR=/home/scratch/yinmt202607/5.angkor_shotgun_finished/data/reads
BT2_FASTQ=/home/scratch/yinmt202607/tests/param_matrix_bt2_vs_bwa/00.extraction/bt2_new/fastq
SLURM_LOG=/home/scratch/yinmt202607/angkor_sub_lib_bams_16/logs/sbatch
mkdir -p "$SLURM_LOG"

[[ -s "$META" ]] || { echo "FATAL: missing $META" >&2; exit 1; }
n=0
while IFS=$'\t' read -r robot shotgun_lib archive core; do
  [[ -n "$robot" ]] || continue
  for spec in "SG:$SHOTGUN_DIR/${shotgun_lib}.prefiltered.IRGSP1.mapped.fq" \
              "C1:$BT2_FASTQ/${robot}_RicePanel1.bt2new.fastq.gz" \
              "C2:$BT2_FASTQ/${robot}_RicePanel2.bt2new.fastq.gz"; do
    tag="${spec%%:*}"; fq="${spec#*:}"
    out="${robot}_${tag}"
    [[ -f "$fq" ]] || { echo "WARN 缺输入 fastq: $out ($fq)" >&2; continue; }
    sbatch --parsable --job-name="sl_${out}" \
      --cpus-per-task "$THREADS" --mem 16G --time 06:00:00 \
      --partition "$PARTITION" --exclude "$EXCLUDE" \
      --output "$SLURM_LOG/${out}.%j.out" \
      "$SCRIPT_DIR/map_irgsp_single_fq.sh" "$fq" "$out" "$THREADS" >/dev/null
    n=$((n+1))
  done
done < "$META"
echo "submitted $n sub-library mapping jobs (16x3=48)"
