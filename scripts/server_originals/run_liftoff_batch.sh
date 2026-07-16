#!/usr/bin/env bash
# ============================================================
# run_liftoff_batch.sh
# 把 MSU7 的57基因(或全部)注释批量投影(liftover)到16个新下载的
# 亚洲栽培稻/野生稻基因组坐标系上，输出各自的 liftoff_from_msu7.gff3
# ============================================================
set -euo pipefail

# ---------- 1. 路径配置(已用真实目录确认过，染色体命名Chr1..Chr12两边一致) ----------
GENOME_ROOT="/home/scratch/yinmt202607/db/16/asian_rice_panel"   # 16个基因组所在目录(已用reorganize脚本整理成扁平结构)
SRC_GFF3="/home/scratch/yinmt202607/db/gene/msu7.gff3"                        # 待投影的MSU7注释
SRC_REF_FASTA="/home/scratch/yinmt202607/db/gene/msu7_pseudomolecule.fna"     # MSU7官方拟分子fasta(all.con)，
                                                                               # 已验证与msu7.gff3的Chr1..Chr12命名一致
THREADS=8
LOGDIR="liftoff_logs"
SUMMARY="liftoff_summary.tsv"

mkdir -p "$LOGDIR"
echo -e "genome\ttotal_features_input\tmapped\tunmapped\tunmapped_pct" > "$SUMMARY"

# ---------- 2. 前置检查 ----------
command -v liftoff >/dev/null 2>&1 || { echo "[ERROR] liftoff 未安装，先 conda install -c bioconda liftoff"; exit 1; }
[[ -s "$SRC_GFF3" ]] || { echo "[ERROR] 找不到源注释 $SRC_GFF3"; exit 1; }
[[ -s "$SRC_REF_FASTA" ]] || { echo "[ERROR] 找不到源参考基因组 $SRC_REF_FASTA"; exit 1; }

TOTAL_FEATURES=$(grep -vc '^#' "$SRC_GFF3" || true)

# ---------- 3. 批量跑 ----------
for genome_dir in "${GENOME_ROOT}"/*/; do
    genome_name=$(basename "$genome_dir")
    target_fasta="${genome_dir}/genome.fna"
    out_gff3="${genome_dir}/liftoff_from_msu7.gff3"
    unmapped_txt="${genome_dir}/unmapped_features.txt"
    log_file="${LOGDIR}/${genome_name}.log"

    if [[ ! -s "$target_fasta" ]]; then
        echo "[SKIP] ${genome_name}: 没找到 genome.fna，跳过"
        continue
    fi

    if [[ -s "$out_gff3" ]]; then
        echo "[SKIP] ${genome_name}: 已存在 liftoff_from_msu7.gff3，跳过(如需重跑先删除)"
    else
        echo "[RUN ] ${genome_name} ..."
        liftoff \
            -g "$SRC_GFF3" \
            -o "$out_gff3" \
            -u "$unmapped_txt" \
            -p "$THREADS" \
            -copies \
            -polish \
            "$target_fasta" \
            "$SRC_REF_FASTA" \
            > "$log_file" 2>&1 \
        || { echo "[FAIL] ${genome_name}: 看日志 $log_file"; continue; }
    fi

    # ---------- 4. 统计每个基因组的投影成功率 ----------
    mapped=$(grep -vc '^#' "$out_gff3" 2>/dev/null | awk '{print $1}' || echo 0)
    # unmapped_features.txt 每行一个feature ID
    unmapped=$(wc -l < "$unmapped_txt" 2>/dev/null || echo 0)
    if [[ "$TOTAL_FEATURES" -gt 0 ]]; then
        pct=$(awk -v u="$unmapped" -v t="$TOTAL_FEATURES" 'BEGIN{printf "%.2f", (u/t)*100}')
    else
        pct="NA"
    fi
    echo -e "${genome_name}\t${TOTAL_FEATURES}\t${mapped}\t${unmapped}\t${pct}%" >> "$SUMMARY"
done

echo ""
echo "全部完成，汇总见 $SUMMARY"
echo "重点看 unmapped_pct 明显偏高的基因组 —— 通常意味着组装质量差/结构变异多，"
echo "后续用57基因列表去定位时对这些基因组的结果要更谨慎。"

# ---------- 5. (可选) 只挑出57个开花基因，生成一份小的检索用GFF ----------
# 如果你已经有 flower_gene.txt (57个 LOC_Os... ID，一行一个)，可以用下面这段
# 从每个 liftoff_from_msu7.gff3 里只抽这57个基因的记录，方便下游脚本直接用：
#
# for genome_dir in "${GENOME_ROOT}"/*/; do
#     genome_name=$(basename "$genome_dir")
#     gff="${genome_dir}/liftoff_from_msu7.gff3"
#     [[ -s "$gff" ]] || continue
#     grep -Ff flower_gene.txt "$gff" > "${genome_dir}/flower57_liftoff.gff3"
# done

