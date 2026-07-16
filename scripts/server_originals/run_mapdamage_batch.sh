#!/usr/bin/env bash
# ============================================================
# run_mapdamage_batch.sh
# 对16个古DNA样本的BAM批量跑 mapDamage2，产出每个样本的损伤profile。
# 这些profile文件后面直接喂给 NGSNGS，让模拟数据的末端损伤模式
# 跟真实样本一致。
# ============================================================
set -euo pipefail

BAM_DIR="/home/scratch/yinmt202607/results/02.irgsp/01.mapping/final"
REF_FASTA="/home/scratch/yinmt202607/db/asian_rice_panel_index/irgsp.fa"   # 建库比对时用的那份参考(跟BAM header里的一致)
OUT_ROOT="mapdamage_out"
LOGDIR="mapdamage_logs"

mkdir -p "$OUT_ROOT" "$LOGDIR"
command -v mapDamage >/dev/null 2>&1 || { echo "[ERROR] mapDamage2 未安装: conda install -c bioconda mapdamage2"; exit 1; }

echo -e "sample\t5pC_to_T_pos1\t3pG_to_A_pos1\tstatus" > mapdamage_summary.tsv

for bam in "${BAM_DIR}"/*.dedup.bam; do
    sample=$(basename "$bam" .dedup.bam)
    outdir="${OUT_ROOT}/${sample}"
    log="${LOGDIR}/${sample}.log"

    echo "[RUN ] ${sample} ..."
    mapDamage \
        -i "$bam" \
        -r "$REF_FASTA" \
        -d "$outdir" \
        --merge-reference-sequences \
        > "$log" 2>&1 \
    || { echo "[FAIL] ${sample}: 看日志 $log"; \
         echo -e "${sample}\tNA\tNA\tFAILED" >> mapdamage_summary.tsv; continue; }

    # 从频率表里抓第一位(read最末端)的 C->T 和 G->A 频率，
    # 快速判断这个样本的古DNA损伤信号强不强
    c2t1=$(awk -F'\t' 'NR==2{print $2}' "${outdir}/5pC_to_T_freq.txt" 2>/dev/null || echo "NA")
    g2a1=$(awk -F'\t' 'NR==2{print $2}' "${outdir}/3pG_to_A_freq.txt" 2>/dev/null || echo "NA")
    echo -e "${sample}\t${c2t1}\t${g2a1}\tOK" >> mapdamage_summary.tsv
done

echo ""
echo "全部完成，汇总见 mapdamage_summary.tsv"
echo "经验判断: 5pC_to_T_pos1 明显 > 0.1 (即>10%)且随位置递减，说明这份数据有典型古DNA损伤特征。"
echo "如果普遍 < 0.03~0.05 且曲线很平，要么现代污染比例较高，要么建库做了UDG处理去除了损伤(需要跟湿实验那边确认protocol)。"
echo ""
echo "喂给 NGSNGS 时，用 \${OUT_ROOT}/<sample>/Stats_out_MCMC_correct_prob.csv"
echo "或直接用 \${OUT_ROOT}/<sample>/5pC_to_T_freq.txt + 3pG_to_A_freq.txt"
echo "(NGSNGS 的 -mf 参数接受 mapDamage 风格的错配频率表，无需再手工换算成Briggs的4个参数)"

