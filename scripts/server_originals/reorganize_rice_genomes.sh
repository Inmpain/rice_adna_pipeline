#!/usr/bin/env bash
# ============================================================
# reorganize_rice_genomes.sh
# 把 NCBI datasets 下载的深层嵌套目录，整理成 liftoff 脚本
# (run_liftoff_batch.sh) 期待的扁平结构:
#   asian_rice_panel/<friendly_name>/genome.fna
# 用软链接(不拷贝)，省磁盘。名字取自 16_3k.csv 的
# "Internal genome numbering" + "Acronyms" 两列，比如 genome15_CM。
# ============================================================
set -euo pipefail

DL_ROOT="/home/scratch/yinmt202607/db/16/asian_rice_panel_download"
CSV="${DL_ROOT}/16_3k.csv"
OUT_ROOT="/home/scratch/yinmt202607/db/16/asian_rice_panel"

[[ -s "$CSV" ]] || { echo "[ERROR] 找不到 $CSV"; exit 1; }
mkdir -p "$OUT_ROOT"

echo -e "accession\tfriendly_name\tstatus" > "${OUT_ROOT}/../reorganize_summary.tsv"

tail -n +2 "$CSV" | while IFS=',' read -r refname acronym genomenum accession rest; do
    acc_dir="${DL_ROOT}/unzipped/${accession}"
    friendly="${genomenum}_${acronym// /}"     # 去掉acronym里可能的空格

    if [[ ! -d "$acc_dir" ]]; then
        echo "[SKIP] ${accession} (${friendly}): 没成功下载，不在 unzipped/ 里"
        echo -e "${accession}\t${friendly}\tMISSING" >> "${OUT_ROOT}/../reorganize_summary.tsv"
        continue
    fi

    fna=$(find "$acc_dir" -name "*_genomic.fna" 2>/dev/null | head -1)
    if [[ -z "$fna" ]]; then
        echo "[WARN] ${accession} (${friendly}): 目录在但找不到 *_genomic.fna"
        echo -e "${accession}\t${friendly}\tNO_FNA" >> "${OUT_ROOT}/../reorganize_summary.tsv"
        continue
    fi

    dest_dir="${OUT_ROOT}/${friendly}"
    mkdir -p "$dest_dir"
    ln -sf "$fna" "${dest_dir}/genome.fna"

    # IRGSP这一份(GCF_001433935.1)在NCBI上自带RefSeq注释，顺手也链过去，
    # 这份gff和这份fna保证坐标严格一致，可以直接当liftoff的源注释用。
    gff=$(find "$acc_dir" -name "genomic.gff" 2>/dev/null | head -1)
    if [[ -n "$gff" ]]; then
        ln -sf "$gff" "${dest_dir}/genomic.gff"
        echo "[INFO] ${accession} (${friendly}): 自带RefSeq注释，已链接 genomic.gff"
    fi

    echo "[OK] ${accession} -> ${friendly}/genome.fna"
    echo -e "${accession}\t${friendly}\tOK" >> "${OUT_ROOT}/../reorganize_summary.tsv"
done

echo ""
echo "完成。结果见 ${OUT_ROOT}"
echo "汇总表: $(dirname "$OUT_ROOT")/reorganize_summary.tsv"
echo ""
ls "$OUT_ROOT"
