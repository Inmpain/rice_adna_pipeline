#!/usr/bin/env bash
set -euo pipefail

# =====================================================================
# 02b_submit_extract_bt2_new.sh
# 对capture panel1/2全部原始文件, 逐个sbatch提交02a的单文件处理。
# =====================================================================

BASE=/home/scratch/yinmt202607/tests/param_matrix_bt2_vs_bwa
SCRIPT_DIR="$BASE/scripts"
SBATCH_LOGDIR="$BASE/logs/sbatch_extract"
THREADS=4

mkdir -p "$BASE/00.extraction/bt2_new"/{bam,fastq} "$SBATCH_LOGDIR" "$BASE/logs"

READ_GLOBS=(
    /home/scratch/yinmt202607/3.angkor_capture_panel1/data/reads/*.bbduk.lowcomp_filtered.fq
    /home/scratch/yinmt202607/7_angor_capture_panel2/data/reads/*.bbduk.lowcomp_filtered.fq
)

n_submitted=0
for pattern in "${READ_GLOBS[@]}"; do
    for fq in $pattern; do
        [[ -f "$fq" ]] || continue
        sample=$(basename "$fq" .bbduk.lowcomp_filtered.fq)

        sbatch --job-name "bt2new_${sample}" \
            --cpus-per-task "$THREADS" \
            --mem 24G \
            --time 01:00:00 \
            --output "${SBATCH_LOGDIR}/${sample}.%j.out" \
            --wrap "bash ${SCRIPT_DIR}/02a_extract_bt2_new_single.sh ${fq} ${THREADS}" \
            > /dev/null

        n_submitted=$((n_submitted+1))
    done
done

echo "共提交 $n_submitted 个作业"
squeue -u "$USER" | wc -l
