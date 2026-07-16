#!/usr/bin/env bash
set -euo pipefail

module load samtools bedtools 2>/dev/null || true

BAM_DIR=/home/scratch/yinmt202607/results/02.irgsp/01.mapping/final
IN_DIR=/home/scratch/yinmt202607/results/02.irgsp/00.reads
BED=/home/scratch/yinmt202607/db/gene/flower_gene.sorted.bed
OUT_DIR=/home/scratch/yinmt202607/results/02.irgsp/02.gene_hits
THREADS=8

mkdir -p "$OUT_DIR"/{bam,fastq,coverage,per_gene}

SUMMARY="$OUT_DIR/sample_summary.tsv"
echo -e "sample\traw_reads\tdedup_mapped\tdedup_q30\tgenome_meandepth\tgenome_pct_covered\tgene_hit_reads" > "$SUMMARY"

# 从q30.bam文件名提取robotid列表 (排除.bai)
ROBOTIDS=$(find "$BAM_DIR" -maxdepth 1 -name "*.dedup.q30.bam" | xargs -n1 basename | sed 's/\.dedup\.q30\.bam$//' | sort -u)

for robot in $ROBOTIDS; do
    echo "=== Processing $robot ==="

    dedup_bam="$BAM_DIR/${robot}.dedup.bam"
    q30_bam="$BAM_DIR/${robot}.dedup.q30.bam"

    if [[ ! -f "$dedup_bam" || ! -f "$q30_bam" ]]; then
        echo "WARNING: missing BAM for $robot, skip"
        continue
    fi

    # --- 原始reads数 (00.reads里同robotid的所有fq/fq.gz求和) ---
    raw=0
    for fq in "$IN_DIR/${robot}."*.fq "$IN_DIR/${robot}_"*.fastq.gz; do
        [[ -f "$fq" ]] || continue
        if [[ "$fq" == *.gz ]]; then
            n=$(( $(zcat "$fq" | wc -l) / 4 ))
        else
            n=$(( $(wc -l < "$fq") / 4 ))
        fi
        raw=$((raw + n))
    done

    # --- mapped reads数 ---
    dedup_mapped=$(samtools view -c "$dedup_bam")
    q30_mapped=$(samtools view -c "$q30_bam")

    # --- 全基因组覆盖度 (基于q30 bam) ---
    samtools coverage "$q30_bam" > "$OUT_DIR/coverage/${robot}.percontig.tsv"
    read meandepth pctcov < <(awk 'NR>1{
        len=$3-$2+1
        sum_len+=len
        sum_cov_bases+=$5
        sum_depth_x_len+=$7*len
    }
    END{
        if(sum_len>0){
            printf "%.4f %.3f\n", sum_depth_x_len/sum_len, sum_cov_bases/sum_len*100
        } else {
            printf "0 0\n"
        }
    }' "$OUT_DIR/coverage/${robot}.percontig.tsv")

    # --- 基因区间命中reads (基于q30 bam) ---
    hit_bam="$OUT_DIR/bam/${robot}.flower_genes.bam"
    samtools view -bh -L "$BED" "$q30_bam" > "$hit_bam"
    samtools index "$hit_bam"
    gene_hit_reads=$(samtools view -c "$hit_bam")

    # 导出命中reads的fastq
    samtools fastq "$hit_bam" > "$OUT_DIR/fastq/${robot}.flower_genes.fq" 2>/dev/null

    # --- 每个基因单独计数 ---
    bedtools intersect -a "$BED" -b "$q30_bam" -c \
        > "$OUT_DIR/per_gene/${robot}.genecounts.tsv"

    echo -e "${robot}\t${raw}\t${dedup_mapped}\t${q30_mapped}\t${meandepth}\t${pctcov}\t${gene_hit_reads}" >> "$SUMMARY"
done

echo "=== 汇总统计完成: $SUMMARY ==="
column -t "$SUMMARY"
