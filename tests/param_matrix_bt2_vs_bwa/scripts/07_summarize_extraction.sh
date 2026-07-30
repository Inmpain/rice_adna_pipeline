#!/usr/bin/env bash
set -euo pipefail

BASE=/home/scratch/yinmt202607/tests/param_matrix_bt2_vs_bwa
SHOTGUN_DIR=/home/scratch/yinmt202607/results/02.irgsp/00.reads_bwa
OUT="$BASE/summary/extraction_summary.tsv"

mkdir -p "$BASE/summary"
echo -e "sample\tmethod\tpanel1_reads\tpanel2_reads\tpanel_total\tshotgun_reads\traw_total" > "$OUT"

count_fastq_reads() {
    local fq="$1"
    if [[ "$fq" == *.gz ]]; then
        echo $(( $(zcat "$fq" | wc -l) / 4 ))
    else
        echo $(( $(wc -l < "$fq") / 4 ))
    fi
}

for method in bt2_old bt2_new bwa; do
    dir="$BASE/00.extraction/$method"
    [[ -d "$dir" ]] || { echo "跳过($method): 目录不存在"; continue; }

    # 只从真正的angkor panel文件名里提取robotid, 排除MCP proxy样本
    # (兼容文件可能直接在$dir下, 也可能在$dir/fastq子目录下)
    robots=""
    for subdir in "$dir" "$dir/fastq"; do
        [[ -d "$subdir" ]] || continue
        for f in "$subdir"/*_RicePanel[12]*; do
            [[ -f "$f" ]] || continue
            r=$(basename "$f" | grep -oE '^LV[0-9]+')
            robots="$robots $r"
        done
    done
    robots=$(echo "$robots" | tr ' ' '\n' | sort -u | grep -v '^$')

    for robot in $robots; do
        p1_total=0
        for subdir in "$dir" "$dir/fastq"; do
            [[ -d "$subdir" ]] || continue
            for f in "$subdir/${robot}"*Panel1*; do
                [[ -f "$f" ]] || continue
                p1_total=$((p1_total + $(count_fastq_reads "$f")))
            done
        done

        p2_total=0
        for subdir in "$dir" "$dir/fastq"; do
            [[ -d "$subdir" ]] || continue
            for f in "$subdir/${robot}"*Panel2*; do
                [[ -f "$f" ]] || continue
                p2_total=$((p2_total + $(count_fastq_reads "$f")))
            done
        done

        panel_total=$((p1_total + p2_total))

        shotgun_fq="$SHOTGUN_DIR/${robot}.prefiltered.IRGSP1.mapped.fq"
        if [[ -f "$shotgun_fq" ]]; then
            shotgun=$(count_fastq_reads "$shotgun_fq")
        else
            shotgun=0
        fi

        raw_total=$((panel_total + shotgun))

        echo -e "${robot}\t${method}\t${p1_total}\t${p2_total}\t${panel_total}\t${shotgun}\t${raw_total}" >> "$OUT"
    done
done

echo "完成: $OUT"
column -t "$OUT"
