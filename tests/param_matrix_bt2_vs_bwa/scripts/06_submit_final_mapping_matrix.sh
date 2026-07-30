#!/usr/bin/env bash
set -euo pipefail

# =====================================================================
# 06_submit_final_mapping_matrix.sh
# 对7个新增组合(①②已有历史数据, 不重新跑), 每个组合x每个样本各自
# 提交一个独立sbatch作业, 调用05_run_final_mapping_single.sh处理。
# 7组合 x 16样本 = 112个作业。
# =====================================================================

BASE=/home/scratch/yinmt202607/tests/param_matrix_bt2_vs_bwa
SCRIPT_DIR="$BASE/scripts"
SBATCH_LOGDIR="$BASE/logs/sbatch"
THREADS=4

mkdir -p "$SBATCH_LOGDIR"

# 需要新跑的7个组合 (①bt2_old_extract__bwa_map 和 ②bwa_extract__bwa_map 已有历史数据, 跳过)
COMBOS=(
    "bt2_old:bt2old"
    "bt2_old:bt2new"
    "bt2_new:bwa"
    "bt2_new:bt2old"
    "bt2_new:bt2new"
    "bwa:bt2old"
    "bwa:bt2new"
)

# 从任意一个reads_combined目录里提取robotid列表
ROBOTS=$(ls "$BASE/01.reads_combined/bwa" 2>/dev/null | sed -E 's/^([^._]+).*/\1/' | sort -u)

if [[ -z "$ROBOTS" ]]; then
    echo "[ERROR] 找不到样本列表, 请先跑 04_prepare_reads_combined.sh"
    exit 1
fi

n_submitted=0
for combo in "${COMBOS[@]}"; do
    extract_method="${combo%%:*}"
    map_tool="${combo##*:}"

    for robot in $ROBOTS; do
        job_name="pm_${extract_method}_${map_tool}_${robot}"

        sbatch --job-name "$job_name" \
            --cpus-per-task "$THREADS" \
            --mem 8G \
            --time 02:00:00 \
            --output "${SBATCH_LOGDIR}/${job_name}.%j.out" \
            --wrap "bash ${SCRIPT_DIR}/05_run_final_mapping_single.sh ${extract_method} ${map_tool} ${robot} ${THREADS}" \
            > /dev/null

        n_submitted=$((n_submitted+1))
    done
    echo "已提交组合: ${extract_method}_extract__${map_tool}_map (16个样本作业)"
done

echo ""
echo "共提交 $n_submitted 个作业"
squeue -u "$USER" | wc -l
