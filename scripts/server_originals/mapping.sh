set -euo pipefail

REF=/home/scratch/yinmt202607/db/asian_rice_panel_index/irgsp.fa
IN_DIR=/home/scratch/yinmt202607/results/02.irgsp/00.reads
OUT_DIR=/home/scratch/yinmt202607/results/02.irgsp/01.mapping
THREADS=20

mkdir -p "$OUT_DIR"/{bam,logs,final}

# 提取所有 robotid（取文件名第一个 "." 或 "_" 之前的部分）
ROBOTIDS=$(ls "$IN_DIR" | sed -E 's/^([^._]+).*/\1/' | sort -u)

for robot in $ROBOTIDS; do
    echo "=== $robot ==="
    mkdir -p "$OUT_DIR/bam/$robot"
    bam_list=()

    for fq in "$IN_DIR/${robot}"*.fq "$IN_DIR/${robot}"*.fastq.gz; do
        [[ -f "$fq" ]] || continue
        base=$(basename "$fq")
        base=${base%.fq}; base=${base%.fastq.gz}
        sai="$OUT_DIR/bam/$robot/${base}.sai"
        bam="$OUT_DIR/bam/$robot/${base}.bam"

        bwa aln -l 1024 -n 0.01 -o 2 -t "$THREADS" "$REF" "$fq" \
            > "$sai" 2> "$OUT_DIR/logs/${base}.aln.log"
        bwa samse "$REF" "$sai" "$fq" 2>> "$OUT_DIR/logs/${base}.aln.log" \
            | samtools view -@ "$THREADS" -bh -F 0x904 - \
            | samtools sort -@ "$THREADS" -o "$bam" -
        samtools index "$bam"
        bam_list+=("$bam")
    done

    # 合并同一robotid的所有bam
    merged="$OUT_DIR/bam/$robot/${robot}.merged.bam"
    samtools merge -@ "$THREADS" -f "$merged" "${bam_list[@]}"

    # 去重（单端数据，用samtools markdup: 需先collate+fixmate）
    samtools collate -@ "$THREADS" -O "$merged" \
        | samtools fixmate -@ "$THREADS" -m - - \
        | samtools sort -@ "$THREADS" -o "$OUT_DIR/bam/$robot/${robot}.sorted.bam" -
    samtools markdup -@ "$THREADS" -r \
        "$OUT_DIR/bam/$robot/${robot}.sorted.bam" \
        "$OUT_DIR/final/${robot}.dedup.bam"
    samtools index "$OUT_DIR/final/${robot}.dedup.bam"

    # MAPQ30 过滤版（用于SNP calling）
    samtools view -@ "$THREADS" -bh -q 30 \
        "$OUT_DIR/final/${robot}.dedup.bam" \
        > "$OUT_DIR/final/${robot}.dedup.q30.bam"
    samtools index "$OUT_DIR/final/${robot}.dedup.q30.bam"

    echo "$robot done: $(samtools view -c "$OUT_DIR/final/${robot}.dedup.bam") dedup reads, \
$(samtools view -c "$OUT_DIR/final/${robot}.dedup.q30.bam") after MAPQ30"
done
