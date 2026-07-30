#!/usr/bin/env bash
set -euo pipefail

# =====================================================================
# 08_summarize_final_mapping.sh
# 对9个组合(2个历史软链接+7个新增)的每个样本, 统计:
#   dedup_total / dup_count / dup_rate_pct / dedup_q30 /
#   genome_meandepth / genome_pct_covered / gene_hit_reads
# =====================================================================

BASE=/home/scratch/yinmt202607/tests/param_matrix_bt2_vs_bwa
FLOWER_BED=/home/scratch/yinmt202607/db/gene/flower_gene.sorted.bed
OUT="$BASE/summary/final_mapping_summary.tsv"
TMP_COV=$(mktemp)

mkdir -p "$BASE/summary"
echo -e "combo\tsample\tdedup_total\tdup_count\tdup_rate_pct\tdedup_q30\tgenome_meandepth\tgenome_pct_covered\tgene_hit_reads" > "$OUT"

# 用find -L 显式follow软链接, 不依赖glob的默认行为
while IFS= read -r combo_dir; do
    combo=$(basename "$combo_dir")

    # 有的组合final结果直接在combo_dir下(软链接场景), 有的在combo_dir/final下(新生成场景)
    if [[ -d "$combo_dir/final" ]]; then
        final_dir="$combo_dir/final"
    else
        final_dir="$combo_dir"
    fi

    echo "=== $combo (数据源: $final_dir) ==="

    for dedup_bam in "$final_dir"/*.dedup.bam; do
        [[ -f "$dedup_bam" ]] || continue
        robot=$(basename "$dedup_bam" .dedup.bam)
        q30_bam="$final_dir/${robot}.dedup.q30.bam"
        [[ -f "$q30_bam" ]] || { echo "  跳过($robot): 缺q30 bam"; continue; }

        dedup_total=$(samtools view -c "$dedup_bam")
        dup_count=$(samtools view -c -f 0x400 "$dedup_bam" 2>/dev/null || echo 0)
        dup_rate=$(awk -v d="$dup_count" -v t="$dedup_total" \
            'BEGIN{ if(t>0) printf "%.2f", d/t*100; else print 0 }')
        q30_count=$(samtools view -c "$q30_bam")

        samtools coverage "$q30_bam" > "$TMP_COV" 2>/dev/null
        read meandepth pctcov < <(awk 'NR>1{
            len=$3-$2+1; sum_len+=len; sum_cov_bases+=$5; sum_depth_x_len+=$7*len
        }
        END{
            if(sum_len>0){printf "%.4f %.3f\n", sum_depth_x_len/sum_len, sum_cov_bases/sum_len*100}
            else {printf "0 0\n"}
        }' "$TMP_COV")

        gene_hit=$(samtools view -c -L "$FLOWER_BED" "$q30_bam" 2>/dev/null || echo 0)

        echo -e "${combo}\t${robot}\t${dedup_total}\t${dup_count}\t${dup_rate}\t${q30_count}\t${meandepth}\t${pctcov}\t${gene_hit}" >> "$OUT"
        echo "  $robot: dedup=${dedup_total} q30=${q30_count} genehit=${gene_hit}"
    done
done < <(find -L "$BASE/02.final_mapping" -mindepth 1 -maxdepth 1 -type d | sort)

rm -f "$TMP_COV"
echo ""
echo "完成: $OUT"
