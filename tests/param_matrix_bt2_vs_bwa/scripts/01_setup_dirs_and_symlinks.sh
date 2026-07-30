#!/usr/bin/env bash
set -euo pipefail

# =====================================================================
# 01_setup_dirs_and_symlinks.sh
# 建立测试专用目录结构, 把已经跑过的历史数据(bt2_old提取/bwa提取/
# bt2_old→bwa_map/bwa→bwa_map这两个已有组合)软链接进来, 后续新生成的
# 数据直接写在这个新目录下, 新旧数据统一放在一起方便对比统计。
# =====================================================================

BASE=/home/scratch/yinmt202607/tests/param_matrix_bt2_vs_bwa

mkdir -p "$BASE"/{00.extraction,01.reads_combined,02.final_mapping,summary,scripts,logs}

echo "=== 建立目录结构 ==="
tree -L 2 "$BASE" 2>/dev/null || find "$BASE" -maxdepth 2

echo ""
echo "=== 软链接阶段①已有的两种提取结果 ==="

if [[ ! -e "$BASE/00.extraction/bt2_old" ]]; then
    ln -s /home/scratch/yinmt202607/results/asian_rice_compare/bowtie2/fastq \
          "$BASE/00.extraction/bt2_old"
    echo "已链接: bt2_old"
else
    echo "跳过(已存在): bt2_old"
fi

if [[ ! -e "$BASE/00.extraction/bwa" ]]; then
    ln -s /home/scratch/yinmt202607/results/asian_rice_compare/bwa/fastq \
          "$BASE/00.extraction/bwa"
    echo "已链接: bwa"
else
    echo "跳过(已存在): bwa"
fi

echo ""
echo "=== 软链接阶段②已有的两个组合 ==="

if [[ ! -e "$BASE/02.final_mapping/bt2_old_extract__bwa_map" ]]; then
    ln -s /home/scratch/yinmt202607/results/02.irgsp/01.mapping/final \
          "$BASE/02.final_mapping/bt2_old_extract__bwa_map"
    echo "已链接: bt2_old_extract__bwa_map"
else
    echo "跳过(已存在): bt2_old_extract__bwa_map"
fi

if [[ ! -e "$BASE/02.final_mapping/bwa_extract__bwa_map" ]]; then
    ln -s /home/scratch/yinmt202607/results/02.irgsp/01.mapping_bwa/final \
          "$BASE/02.final_mapping/bwa_extract__bwa_map"
    echo "已链接: bwa_extract__bwa_map"
else
    echo "跳过(已存在): bwa_extract__bwa_map"
fi

echo ""
echo "完成。目录结构:"
find "$BASE" -maxdepth 2 -exec ls -ld {} \;
