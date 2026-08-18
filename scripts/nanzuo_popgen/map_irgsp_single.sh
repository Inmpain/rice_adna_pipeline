#!/usr/bin/env bash
set -euo pipefail

# =====================================================================
# map_irgsp_single.sh
# 阶段② 最终定量比对: 合并后的 FASTQ 比对到 irgsp.fa, markdup 只标记不删除,
# 产出 q20/25/30 计数 + 基因组覆盖度统计(全部只统计, 不删 reads)。
#
# 用法: bash map_irgsp_single.sh <sample> <mapper> [threads]
#   mapper: bwa | bt2new
#   输出:
#     nanzuo/02.map_irgsp/{mapper}/{sample}.dedup.bam(.bai)   # 重复只标记(0x400), 未删除
#     nanzuo/02.map_irgsp/{mapper}/stats/{sample}.tsv
# =====================================================================

sample="$1"
mapper="$2"
THREADS="${3:-20}"

module load bwa/ 2>/dev/null || module load bwa 2>/dev/null || true
module load bowtie2/ 2>/dev/null || module load bowtie2 2>/dev/null || true
module load samtools/ 2>/dev/null || module load samtools 2>/dev/null || true

BASE=/home/scratch/yinmt202607/nanzuo
IRGSP_REF=/home/scratch/yinmt202607/db/asian_rice_panel_index/irgsp.fa
IRGSP_BT2_IDX=/home/scratch/yinmt202607/db/asian_rice_panel_index/irgsp_bt2idx
MERGE_DIR="$BASE/01.merge"
OUT_DIR="$BASE/02.map_irgsp/$mapper"
STATS_DIR="$OUT_DIR/stats"
LOG_DIR="$BASE/_logs/map_irgsp"

# 0x004 unmapped + 0x100 secondary + 0x800 supplementary
EXCLUDE_FLAGS=0x904
BWA_ALN_PARAMS=(-l 1024 -n 0.01 -o 2)
# Bowtie2 新参数(-N1), 与 9 格矩阵测试的 bt2_new 一致
BT2_NEW_PARAMS=(-k 3 -L 22 -N 1 -i S,1,1.15 --mp 1,1 --rdg 0,1 --rfg 0,1 --score-min L,0,-0.1 --no-unal)

case "$mapper" in
  bwa|bt2new) ;;
  *) echo "ERROR: mapper 必须是 bwa 或 bt2new" >&2; exit 1 ;;
esac

fq="$MERGE_DIR/${sample}.combined.fastq.gz"
[[ -s "$fq" ]] || { echo "[$sample] ERROR: 缺 $fq" >&2; exit 1; }

out_bam="$OUT_DIR/${sample}.dedup.bam"
if [[ -s "$out_bam" && -s "$STATS_DIR/${sample}.tsv" ]]; then
  echo "[$sample][$mapper] 已存在, 跳过"
  exit 0
fi

mkdir -p "$OUT_DIR" "$STATS_DIR" "$LOG_DIR"
work="$(mktemp -d "$OUT_DIR/.${sample}.XXXXXX")"
trap 'rm -rf "$work"' EXIT

mapped_bam="$work/mapped.bam"

if [[ "$mapper" == "bwa" ]]; then
  bwa aln "${BWA_ALN_PARAMS[@]}" -t "$THREADS" "$IRGSP_REF" "$fq" 2> "$LOG_DIR/${sample}.bwa.aln.log" \
    | bwa samse "$IRGSP_REF" - "$fq" 2>> "$LOG_DIR/${sample}.bwa.aln.log" \
    | samtools view -@ "$THREADS" -bh -F "$EXCLUDE_FLAGS" - \
    | samtools sort -@ "$THREADS" -o "$mapped_bam" - 2>> "$LOG_DIR/${sample}.bwa.aln.log"
else
  bowtie2 -p "$THREADS" "${BT2_NEW_PARAMS[@]}" -x "$IRGSP_BT2_IDX" -U "$fq" 2> "$LOG_DIR/${sample}.bt2new.aln.log" \
    | samtools view -@ "$THREADS" -bh -F "$EXCLUDE_FLAGS" - \
    | samtools sort -@ "$THREADS" -o "$mapped_bam" - 2>> "$LOG_DIR/${sample}.bt2new.aln.log"
fi

# markdup 只标记不删除(无 -r): collate + fixmate -m + sort + markdup
samtools collate -@ "$THREADS" -O "$mapped_bam" \
  | samtools fixmate -@ "$THREADS" -m - - \
  | samtools sort -@ "$THREADS" -o "$work/sorted.bam" -
samtools markdup -@ "$THREADS" "$work/sorted.bam" "$out_bam"
samtools index -@ "$THREADS" "$out_bam"

# ---- 统计(只计数, 不删 reads) ----
merged_reads=$(( $(gzip -cd "$fq" | wc -l) / 4 ))
primary_mapped=$(samtools view -c -F "$EXCLUDE_FLAGS" "$mapped_bam")
dup_flagged=$(samtools view -c -f 0x400 "$out_bam")
# q20/25/30 口径同主线: 排除重复(0x400) + MAPQ 阈值, 不物理删除
q20=$(samtools view -c -q 20 -F 0x400 "$out_bam")
q25=$(samtools view -c -q 25 -F 0x400 "$out_bam")
q30=$(samtools view -c -q 30 -F 0x400 "$out_bam")

# 基因组覆盖度: 基于去重后的 reads
samtools view -@ "$THREADS" -b -F 0x400 "$out_bam" -o "$work/uniq.bam"
samtools index -@ "$THREADS" "$work/uniq.bam"
samtools coverage "$work/uniq.bam" > "$work/coverage.tsv" 2>/dev/null
read cov_bases cov_pct mean_depth < <(awk 'NR>1{
    len=$3-$2+1
    sum_len+=len
    sum_cov+=$5
    sum_depth_x_len+=$7*len
  }
  END{
    if(sum_len>0) printf "%.0f %.4f %.4f\n", sum_cov, sum_cov/sum_len*100, sum_depth_x_len/sum_len
    else print "0 0 0"
  }' "$work/coverage.tsv")

dup_rate=$(awk -v d="$dup_flagged" -v m="$primary_mapped" 'BEGIN{ if(m>0) printf "%.4f", d/m*100; else print "0" }')

printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
  "$sample" "$mapper" "$merged_reads" "$primary_mapped" "$dup_flagged" "$dup_rate" \
  "$q20" "$q25" "$q30" "$cov_bases" "$cov_pct" "$mean_depth" \
  > "$STATS_DIR/${sample}.tsv"

echo "[$sample][$mapper] 完成: mapped=$primary_mapped dup=$dup_flagged q30=$q30 cov_bases=$cov_bases"
