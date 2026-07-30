#!/usr/bin/env bash
set -euo pipefail

# =====================================================================
# 05_run_final_mapping_single.sh
# 处理单个 (提取方法 x 定量比对工具 x 样本) 组合:
#   比对该样本的shotgun+panel1+panel2三个来源 → merge → 去重(标记不删除)
#   → q30过滤
# 设计成单样本粒度, 方便配合sbatch按(combo, sample)逐个提交, 最大化并行度。
#
# 用法: bash 05_run_final_mapping_single.sh <extract_method> <map_tool> <robot> [threads]
#   extract_method: bt2_old | bt2_new | bwa
#   map_tool:       bwa | bt2old | bt2new
#   robot:          样本robotid, 如 LV7008416379
# =====================================================================

extract_method="$1"
map_tool="$2"
robot="$3"
THREADS="${4:-4}"

BASE=/home/scratch/yinmt202607/tests/param_matrix_bt2_vs_bwa
IRGSP_BWA_REF=/home/scratch/yinmt202607/db/asian_rice_panel_index/irgsp.fa
IRGSP_BT2_IDX=/home/scratch/yinmt202607/db/asian_rice_panel_index/irgsp_bt2idx

BT2_OLD_EXTRA=( -k 3 -L 22 -i S,1,1.15 --mp 1,1 --rdg 0,1 --rfg 0,1 --score-min L,0,-0.1 --no-unal )
BT2_NEW_EXTRA=( -k 3 -L 22 -N 1 -i S,1,1.15 --mp 1,1 --rdg 0,1 --rfg 0,1 --score-min L,0,-0.1 --no-unal )

combo_name="${extract_method}_extract__${map_tool}_map"
in_dir="$BASE/01.reads_combined/$extract_method"
out_dir="$BASE/02.final_mapping/$combo_name"

mkdir -p "$out_dir/bam" "$out_dir/final"

# 如果这个组合已经是软链接指向历史数据(①②两个已有组合), 直接退出不重新跑
if [[ -L "$BASE/02.final_mapping/$combo_name" ]]; then
    echo "[$combo_name][$robot] 该组合是软链接历史数据, 跳过"
    exit 0
fi

if [[ -s "$out_dir/final/${robot}.dedup.q30.bam" ]]; then
    echo "[$combo_name][$robot] 已完成, 跳过"
    exit 0
fi

echo "[$combo_name][$robot] 开始处理..."

bam_list=()
for fq in "$in_dir/${robot}"*.fq "$in_dir/${robot}"*.fastq.gz; do
    [[ -f "$fq" ]] || continue
    base=$(basename "$fq"); base=${base%.fq}; base=${base%.fastq.gz}
    bam="$out_dir/bam/${base}.bam"

    if [[ "$map_tool" == "bwa" ]]; then
        bwa aln -l 1024 -n 0.01 -o 2 -t "$THREADS" "$IRGSP_BWA_REF" "$fq" \
            2> "$out_dir/bam/${base}.aln.log" \
            | bwa samse "$IRGSP_BWA_REF" - "$fq" 2>> "$out_dir/bam/${base}.aln.log" \
            | samtools view -@ "$THREADS" -bh -F 0x904 - \
            | samtools sort -@ "$THREADS" -o "$bam" -
    else
        if [[ "$map_tool" == "bt2old" ]]; then
            extra_params=("${BT2_OLD_EXTRA[@]}")
        else
            extra_params=("${BT2_NEW_EXTRA[@]}")
        fi
        bowtie2 -p "$THREADS" "${extra_params[@]}" -x "$IRGSP_BT2_IDX" -U "$fq" \
            2> "$out_dir/bam/${base}.aln.log" \
            | samtools view -@ "$THREADS" -bh -F 0x904 - \
            | samtools sort -@ "$THREADS" -o "$bam" -
    fi
    samtools index "$bam"
    bam_list+=("$bam")
done

if [[ ${#bam_list[@]} -eq 0 ]]; then
    echo "[$combo_name][$robot] 警告: 没有找到任何输入reads文件"
    exit 1
fi

merged="$out_dir/bam/${robot}.merged.bam"
samtools merge -@ "$THREADS" -f "$merged" "${bam_list[@]}"

samtools collate -@ "$THREADS" -O "$merged" \
    | samtools fixmate -@ "$THREADS" -m - - \
    | samtools sort -@ "$THREADS" -o "$out_dir/bam/${robot}.sorted.bam" -

samtools markdup -@ "$THREADS" \
    "$out_dir/bam/${robot}.sorted.bam" "$out_dir/final/${robot}.dedup.bam"
samtools index "$out_dir/final/${robot}.dedup.bam"

samtools view -@ "$THREADS" -bh -q 30 -F 0x400 \
    "$out_dir/final/${robot}.dedup.bam" > "$out_dir/final/${robot}.dedup.q30.bam"
samtools index "$out_dir/final/${robot}.dedup.q30.bam"

total=$(samtools view -c "$out_dir/final/${robot}.dedup.bam")
q30=$(samtools view -c "$out_dir/final/${robot}.dedup.q30.bam")
echo "[$combo_name][$robot] 完成: dedup_total=${total}, q30=${q30}"
