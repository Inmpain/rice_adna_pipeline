
#!/usr/bin/env bash

set -euo pipefail

# =====================================================================

# 04_prepare_reads_combined.sh

# 对3种提取方法(bt2_old/bt2_new/bwa), 各自把固定不变的shotgun reads +

# 对应方法提取出的panel1/2 reads, 软链接整理到统一目录, 供阶段②比对使用。

#

# 兼容两种目录结构: 提取产出可能直接在 00.extraction/{method}/ 这层,

# 也可能在 00.extraction/{method}/fastq/ 子目录下。

#

# ⚠️关键过滤: bt2_old/bwa的历史提取产出目录里混杂了MCP reshotgun proxy

# 样本(因为当初compare_rice_read_extractors.sh默认扫描范围包含了

# 4.mcp_reshotgun), 必须显式过滤掉, 只保留真正的angkor capture panel1/2

# 文件。真实panel文件名格式: LV{数字}_RicePanel{1|2}...

# MCP proxy文件名格式: LV{数字}-LV{数字}-proxy...(用短横线, 不是下划线)

# =====================================================================

BASE=/home/scratch/yinmt202607/tests/param_matrix_bt2_vs_bwa

SHOTGUN_SRC=/home/scratch/yinmt202607/results/02.irgsp/00.reads_bwa

for method in bt2_old bt2_new bwa; do

    src="$BASE/00.extraction/$method"

    dst="$BASE/01.reads_combined/$method"

    mkdir -p "$dst"

    if [[ ! -d "$src" && ! -L "$src" ]]; then

        echo "跳过($method): 提取结果目录不存在 $src"

        continue

    fi

    echo "=== $method ==="

    n_shotgun=0

    for f in "$SHOTGUN_SRC"/*.prefiltered.IRGSP1.mapped.fq; do

        [[ -f "$f" ]] || continue

        ln -sf "$f" "$dst/$(basename "$f")"

        n_shotgun=$((n_shotgun+1))

    done

    n_panel=0

    n_skipped_mcp=0

    for subdir in "$src" "$src/fastq"; do

        [[ -d "$subdir" ]] || continue

        for f in "$subdir"/*.gz "$subdir"/*.fq; do

            [[ -f "$f" ]] || continue

            base=$(basename "$f")

            # 只保留真正的angkor panel1/2文件, 排除MCP proxy样本

            if [[ "$base" =~ ^LV[0-9]+_RicePanel[0-9] ]]; then

                ln -sf "$f" "$dst/$base"

                n_panel=$((n_panel+1))

            else

                n_skipped_mcp=$((n_skipped_mcp+1))

            fi

        done

    done

    echo "  shotgun链接数: $n_shotgun, panel链接数: $n_panel, 跳过(非angkor如MCP): $n_skipped_mcp"

done

echo ""

echo "完成。检查各目录文件数:"

for method in bt2_old bt2_new bwa; do

    echo "$method: $(ls "$BASE/01.reads_combined/$method" 2>/dev/null | wc -l) 个文件"

done

