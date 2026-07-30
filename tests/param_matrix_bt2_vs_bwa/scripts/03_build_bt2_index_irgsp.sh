#!/usr/bin/env bash
set -euo pipefail

# =====================================================================
# 03_build_bt2_index_irgsp.sh
# 阶段②(定量比对)首次要用到Bowtie2, 之前irgsp.fa只建过BWA索引,
# 这里补建Bowtie2索引。只需要跑一次。
# =====================================================================

REF=/home/scratch/yinmt202607/db/asian_rice_panel_index/irgsp.fa
OUT_PREFIX=/home/scratch/yinmt202607/db/asian_rice_panel_index/irgsp_bt2idx

command -v bowtie2-build >/dev/null 2>&1 || { echo "[ERROR] bowtie2-build 未找到"; exit 1; }
[[ -f "$REF" ]] || { echo "[ERROR] 找不到参考基因组: $REF"; exit 1; }

if [[ -s "${OUT_PREFIX}.1.bt2" ]]; then
    echo "索引已存在，跳过: ${OUT_PREFIX}.1.bt2"
    exit 0
fi

echo "开始建索引..."
bowtie2-build "$REF" "$OUT_PREFIX"

echo "完成。索引文件:"
ls -la "${OUT_PREFIX}"*.bt2
