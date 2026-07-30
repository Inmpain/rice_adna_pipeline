#!/usr/bin/env bash
set -euo pipefail

# =====================================================================
# 10_recompute_dup_for_historical_combos.sh
# 针对①bt2_old_extract__bwa_map 和 ②bwa_extract__bwa_map 这两个历史组合
# (原流程用markdup -r直接删除重复, 无法事后统计dup_count), 回到markdup
# 之前的中间产物(.sorted.bam), 重新跑一次不加-r的markdup, 补算出真实的
# dup_count/dup_rate, 输出到独立文件, 不覆盖原始历史final目录。
# =====================================================================

THREADS=8

declare -A HIST_COMBOS=(
    ["bt2_old_extract__bwa_map"]="/home/scratch/yinmt202607/results/02.irgsp/01.mapping/bam"
    ["bwa_extract__bwa_map"]="/home/scratch/yinmt202607/results/02.irgsp/01.mapping_bwa/bam"
)

OUT_ROOT=/home/scratch/yinmt202607/tests/param_matrix_bt2_vs_bwa/dup_recompute
mkdir -p "$OUT_ROOT"

OUT_TSV=/home/scratch/yinmt202607/tests/param_matrix_bt2_vs_bwa/summary/historical_dup_recompute.tsv
echo -e "combo\tsample\tdedup_total_recheck\tdup_count\tdup_rate_pct" > "$OUT_TSV"

for combo in "${!HIST_COMBOS[@]}"; do
    bam_root="${HIST_COMBOS[$combo]}"
    out_dir="$OUT_ROOT/$combo"
    mkdir -p "$out_dir"

    echo "=== $combo (源: $bam_root) ==="

    for sorted_bam in "$bam_root"/*/*.sorted.bam; do
        [[ -f "$sorted_bam" ]] || continue
        robot=$(basename "$(dirname "$sorted_bam")")
        marked_bam="$out_dir/${robot}.marked.bam"

        if [[ ! -s "$marked_bam" ]]; then
            samtools markdup -@ "$THREADS" "$sorted_bam" "$marked_bam"
            samtools index "$marked_bam"
        fi

        total=$(samtools view -c "$marked_bam")
        dup=$(samtools view -c -f 0x400 "$marked_bam")
        rate=$(awk -v d="$dup" -v t="$total" 'BEGIN{ if(t>0) printf "%.2f", d/t*100; else print 0 }')

        echo -e "${combo}\t${robot}\t${total}\t${dup}\t${rate}" >> "$OUT_TSV"
        echo "  $robot: total=${total} dup=${dup} rate=${rate}%"
    done
done

echo ""
echo "完成: $OUT_TSV"
