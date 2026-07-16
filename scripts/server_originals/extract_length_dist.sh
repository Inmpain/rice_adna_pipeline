
#!/usr/bin/env bash
# ============================================================
# extract_length_dist.sh
# 从16个真实古DNA样本的BAM里提取read长度的CDF分布文件，
# 格式直接匹配 NGSNGS 的 -lf 参数要求(length \t cumulative_prob)。
# 同时记录每个样本的reads总数，供后面 NGSNGS 的 -r 参数使用。
# ============================================================
set -euo pipefail

BAM_DIR="/home/scratch/yinmt202607/results/02.irgsp/01.mapping/final"
OUT_DIR="readlen_dist"
mkdir -p "$OUT_DIR"

echo -e "sample\ttotal_reads\tmean_len\tmedian_len\tcdf_file" > readlen_summary.tsv

for bam in "${BAM_DIR}"/*.dedup.bam; do
    sample=$(basename "$bam" .dedup.bam)
    raw_counts="${OUT_DIR}/${sample}.raw_counts.txt"
    cdf_file="${OUT_DIR}/${sample}.length_cdf.txt"

    echo "[RUN] ${sample} ..."

    # 每个长度出现的次数(未排序前是length，之后按长度升序排)
    samtools view "$bam" | awk '{print length($10)}' | sort -n | uniq -c \
        | awk '{print $2"\t"$1}' > "$raw_counts"

    # 转成NGSNGS要的CDF格式: length \t cumulative_probability(累加到1)
    awk -F'\t' '
        {len[NR]=$1; cnt[NR]=$2; total+=$2}
        END{
            cum=0
            for(i=1;i<=NR;i++){
                cum+=cnt[i]
                printf "%d\t%.10g\n", len[i], cum/total
            }
        }
    ' "$raw_counts" > "$cdf_file"

    total_reads=$(awk -F'\t' '{s+=$2} END{print s}' "$raw_counts")
    mean_len=$(awk -F'\t' '{s+=$1*$2; n+=$2} END{printf "%.1f", s/n}' "$raw_counts")
    # 中位数: 找累积概率首次达到0.5的那个长度
    median_len=$(awk -F'\t' '$2>=0.5{print $1; exit}' "$cdf_file")

    echo -e "${sample}\t${total_reads}\t${mean_len}\t${median_len}\t${cdf_file}" >> readlen_summary.tsv
    echo "  reads=${total_reads} mean_len=${mean_len} median_len=${median_len}"
done

echo ""
echo "完成，汇总见 readlen_summary.tsv，每个样本的CDF文件在 ${OUT_DIR}/ 下"
echo "后面跑NGSNGS时，每个样本对应用:"
echo "  -r <total_reads>  -lf ${OUT_DIR}/<sample>.length_cdf.txt"

