#!/usr/bin/env bash
set -euo pipefail

# =====================================================================
# collect_stats.sh
# 汇总每个样本 × mapper 的统计表为两张总表。登录节点跑, 不 sbatch。
#
# 用法: bash collect_stats.sh
#   输出:
#     nanzuo/03.stats/read_qc.tsv
#     nanzuo/03.stats/coverage_summary.tsv
# =====================================================================

BASE=/home/scratch/yinmt202607/nanzuo
STATS_ROOT="$BASE/02.map_irgsp"
OUT_DIR="$BASE/03.stats"
mkdir -p "$OUT_DIR"

READ_QC="$OUT_DIR/read_qc.tsv"
COV_SUM="$OUT_DIR/coverage_summary.tsv"

printf 'sample\tmapper\tmerged_reads\tprimary_mapped\tdup_flagged\tdup_rate_pct\tq20\tq25\tq30\n' > "$READ_QC"
printf 'sample\tmapper\tcov_bases\tcov_pct\tmean_depth\n' > "$COV_SUM"

for f in "$STATS_ROOT"/bwa/stats/*.tsv "$STATS_ROOT"/bt2new/stats/*.tsv; do
  [[ -f "$f" ]] || continue
  cut -f1-9 "$f" >> "$READ_QC"
  cut -f1,2,10,11,12 "$f" >> "$COV_SUM"
done

# 排序(保留表头在第一行)
{ head -1 "$READ_QC"; tail -n +2 "$READ_QC" | sort -k2,2 -k1,1; } > "$READ_QC.tmp" && mv "$READ_QC.tmp" "$READ_QC"
{ head -1 "$COV_SUM"; tail -n +2 "$COV_SUM" | sort -k2,2 -k1,1; } > "$COV_SUM.tmp" && mv "$COV_SUM.tmp" "$COV_SUM"

echo "=== read_qc.tsv ==="
column -t "$READ_QC"
echo
echo "=== coverage_summary.tsv ==="
column -t "$COV_SUM"
echo
echo "汇总完成:"
echo "  $READ_QC"
echo "  $COV_SUM"
