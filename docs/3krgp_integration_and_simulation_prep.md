# 3K RG数据整合 + 模拟数据准备工作记录

记录范围：从"OsGI精确坐标查询"开始，到"NGSNGS模拟数据环境搭建、mapDamage2
sbatch并行化"为止。承接 `flank1kb_msa_exploration.md`（该文档记录到DTH8/Ghd8
1116bp缺失簇的精确统计为止）。

最后更新: 2026-07-17

---

## 一、OsGI精确causal坐标查询——结论：没有单点causal突变

查文献未找到OsGI像Ghd7/DTH8/OsPRR37那样"一篇论文报道具体causal突变"的经典
故事。OsGI更多是作为**GWAS关联位点**出现（如181份aus种质研究里，OsGI是抽穗期
相关主成分PC1的显著关联基因，并在3K RG panel做过单倍型分析），但没有指出
具体是哪个SNP/碱基变化导致功能差异。

**结论**：即使reads覆盖到OsGI基因内部，也没有已知"精确causal坐标"可供对照
判断功能型/缺失型。退而求其次的方案：直接查这些坐标在3K RG里本来存在哪些
等位变异、aus与temperate japonica间频率是否有差异，不追求"causal"，只看
"是否存在自然分化"。

---

## 二、下游分析路径重新梳理——从"单基因深挖"转向"群体尺度"

回顾研究目标四层框架的实际进展：

- **第一层(可行性筛查)**：已完成。结论：现有深度不支持单基因层面可靠判断，
  另发现约35.7%命中数据本身有低复杂度序列比对歧义问题。
- **第二层(功能等位型判断)**：路径接近走到头。即使是证据"最干净"的OsGI，
  也因缺乏已知causal位点定义，无法给出确定性结论。
- **第三层(种植生态型推断)**：真正能推进的部分，且不需要精确到单基因单位点。
- **第四层(考古学意义)**：依赖前三层结果。

**建议**：重心从"逐个基因、逐条read找证据"转向"用全基因组尺度的群体遗传学
工具(PCA/f-statistics)直接回答生态型归属问题"——前者受限于深度和序列复杂度，
天花板已摸到；后者不依赖单点精度，是这批低深度古DNA数据真正能发挥价值之处。

---

## 三、3K RG SNP矩阵格式确认与提取

### 数据来源

```
/home/scratch/yinmt202607/db/3k/tmp/
├── 3kall_snpposition_map.tsv   # SNP_INDEX / CHROMOSOME / POSITION / REFCALL
├── 3kall_variety_map.tsv        # VARIETY_INDEX / NAME / IRIS_UNIQUE_ID / BOX_CODE / GS_ACCESSION
├── Universe_matrix_geno_NB      # 核心矩阵: 每行1个品种(共3024行), 每行是按SNP_INDEX
│                                  顺序拼接的完整基因型字符串(3200多万字符/行, "?"=缺失)
└── Universe_matrix_geno_NB.gz
```

README确认：这批SNP是**"called vs Nipponbare MSU7/IRGSP1.0 genome"**，与
`irgsp.fa`坐标系完全一致，无需转换。

### 提取57基因区间内的SNP

```bash
GENE_BED=/home/scratch/yinmt202607/db/gene/flower_gene.sorted.bed
SNP_MAP=/home/scratch/yinmt202607/db/3k/tmp/3kall_snpposition_map.tsv
OUT_DIR=/home/scratch/yinmt202607/results/04.3krgp_flowergenes
mkdir -p "$OUT_DIR"

tail -n+2 "$SNP_MAP" | awk -F'\t' 'BEGIN{OFS="\t"}{
    chr = sprintf("chr%02d", $2)
    print chr, $3-1, $3, $1
}' > "$OUT_DIR/3krgp_snp_positions.bed"
# 结果: 32,064,217行 (全部SNP)

bedtools intersect -a "$OUT_DIR/3krgp_snp_positions.bed" -b "$GENE_BED" -wa -wb \
    > "$OUT_DIR/flower_genes_3krgp_snps.tsv"
# 结果: 12,741个SNP落在57基因区间内
```

### 我们自己的aDNA reads是否精确匹配已知SNP位点

```bash
tail -n+2 /home/scratch/yinmt202607/results/02.irgsp/02.gene_hits/gene_hits_with_metadata.tsv \
    | awk -F'\t' 'BEGIN{OFS="\t"}{print $7, $8-1, $8, $1"|"$6}' \
    > /tmp/our_hit_positions.bed

bedtools intersect -a /tmp/our_hit_positions.bed \
    -b "$OUT_DIR/flower_genes_3krgp_snps.tsv" -wa -wb \
    > "$OUT_DIR/our_reads_matching_3krgp_snps.tsv"
```

**结果：42条命中reads里，6条精确匹配上3K RG已知SNP位点**：

```
SDG711/OsCLF (chr06:9354034)   PHYA (chr03:29172084)   GW2 (chr02:8121336)
OsSnRK1A (chr05:26349205)      Spl11 (chr12:23469623)  OsMYB8 (chr01:25593404)
```

**⚠️关键发现：DTH8/Ghd8、OsGI、DTH7/OsPRR37三个最初优先基因，一个都没有精确
匹配上**——这是低深度数据的必然结果，reads稀疏撒在基因组上，精确撞上已知
SNP的概率本身很低。这条路径目前只能提供次要基因(PHYA、GW2等)的信息。

### 提取矩阵基因型 —— 用 `cut -c` 而非awk内部getline(更可靠)

```bash
VARIETY_MAP=/home/scratch/yinmt202607/db/3k/tmp/3kall_variety_map.tsv
MATRIX=/home/scratch/yinmt202607/db/3k/tmp/Universe_matrix_geno_NB

mkdir -p "$OUT_DIR/per_snp_columns"

while read -r idx; do
    char_pos=$((idx + 1))   # SNP_INDEX是0-based, 字符位置是1-based
    cut -c${char_pos} "$MATRIX" > "$OUT_DIR/per_snp_columns/snp_${idx}.txt"
done < "$OUT_DIR/target_snp_indices.txt"

for idx in $(cat "$OUT_DIR/target_snp_indices.txt"); do
    paste "$VARIETY_MAP" "$OUT_DIR/per_snp_columns/snp_${idx}.txt" \
        > "$OUT_DIR/snp_${idx}_with_variety.tsv"
done
```

⚠️ **踩坑记录**：最初用一个awk内部`getline`一次性提取多个位点的写法失败了
(输出只有行号没有碱基字符)，原因未完全定位清楚，改用`cut -c`逐位点单独提取
后正常工作，更简单也更容易排查。

### 关于要不要做MAF过滤——不应该做

**MAF过滤**用于群体结构分析(PCA/admixture)，是为了排除低频变异的噪音干扰。
但这里做的是**单点诊断性核对**，目的是看古稻样本携带的碱基在3024份材料里
常见还是罕见——**罕见变异恰恰可能是最有价值的信号**，MAF过滤会直接把这类
信息踢掉。**结论：单点核对不做MAF过滤是对的，不是漏做**。但仍需做基础QC
(缺失率检查、多等位位点确认)。

---

## 四、3K RG结构变异(SV)数据 —— Ghd7/DTH8/DTH7三个基因区域对比

### 背景：为什么无法直接从我们自己的reads检测结构变异

- 单端测序，无双端不一致信号可用
- read长度仅30-70bp，物理上无法跨越几百至上千bp的大片段缺失断点
- 背景深度0.0002-0.001x，"看起来没有reads"本身就是低深度下的正常现象，
  无法区分"真实缺失"与"恰好没测到"

**结论**：结构变异检测的三大标准信号来源(双端、split-read、深度下降)在现有
数据条件下均不适用，这是数据物理特性决定的天花板，不是流程问题。

### 但3K RG的SV数据仍有价值——变成"已知变异坐标数据库"

```bash
mkdir -p /home/scratch/yinmt202607/db/3k/sv_extracted
cd /home/scratch/yinmt202607/db/3k/sv_extracted
tar xzf ../NB_DEL_mergesam_clustered.tar.gz
# 格式: 染色体 起点 终点 DEL;长度 品种名 clusterID
```

### 三个基因区域的DEL记录对比

**DTH8/Ghd8**(chr08:4332716-4336434)：发现一个**~1116-1359bp的大片段缺失簇**
(断点坐标因区域微同源性有抖动，判定为同一生物学事件)：

```bash
awk -F'\t' '$1=="chr08" && $2>=4332716 && $3<=4336434' NB_DEL_mergesam_clustered.txt \
    | awk -F'[;\t]' '{print $5}' | sort -n | uniq -c | sort -rn
#   107 12     (独立小indel, 背景多态性)
#    88 1117   ┐
#    62 1116   ├ 判定为同一缺失事件
#     2 1224   │ (断点微同源性导致坐标抖动)
#     1 1359   ┘ 等
#    15 19     (独立小indel, 位置在缺失簇外)
```

**精确统计**：163/3024份材料(约5.4%)携带这个大片段缺失型等位基因，与文献
报道的DTH8经典功能丧失型高度吻合。这是三个基因里唯一发现"大到有意义、
样本量充足"的结构变异候选。

**Ghd7**(chr07:9151401-9156185)：最常见是`DEL;24`(24bp)，100多个品种，属于
背景多态性，与文献报道的CCT结构域截断机制(移码/无义突变)量级不符，非causal候选。

**DTH7/OsPRR37**(chr07:29615704-29630223)：两个热点，`DEL;26`(chr07:29616463-
29616489)和`DEL;18`(chr07:29628375-29628393)，均为极大量品种共享的背景小
indel，与已知causal机制(CCT结构域错义点突变)不符，非causal候选。

### DTH8断点检验结果——16个古稻样本全部零覆盖

```bash
BAM_DIR=/home/scratch/yinmt202607/results/02.irgsp/01.mapping_bwa/final
for robot in $(ls "$BAM_DIR" | grep '\.dedup\.q30\.bam$' | sed 's/\.dedup\.q30\.bam$//'); do
    samtools view -c "$BAM_DIR/${robot}.dedup.q30.bam" chr08:4332834-4333950
done
# 全部16个样本 = 0
```

**结论：这是"零信息"，不是"阴性证据"**——现有DTH8覆盖坐标(4334006-4334634)
全部落在缺失区间(4332834-4333950)下游之外，现有数据对此断点完全没有观测力。
断点坐标已存档，供未来更深capture测序验证使用。

---

## 五、NGSNGS模拟数据环境搭建

### 目的

用5-16个已知生态型的现代/近现代品种，制造"数字孪生古样本"——用与真实古稻
样本一致的技术特征(深度、read长度分布、损伤模式)包装，作为PCA投影里带真值
标签的**校准点**，用于验证smartpca投影方法本身是否可信。

### 更优的锚点选择：16基因组liftoff-ready panel（优于原5基因组panel）

```
/home/scratch/yinmt202607/db/16/asian_rice_panel/
├── genome1_IRGSP  genome4_N22   genome5_AZ    genome6_IR64
├── genome7_ARC    genome8_LM    genome9_LX    genome10_KYG
├── genome11_LIMA  genome12_NABO genome13_PR106 genome14_KN
├── genome15_CM    genome16_GS   genome27_MH63  genome28_ZS97
```

每个目录下已有`genome.fna`(软链接到NCBI datasets原始文件)、`.fai`、`.mmi`索引，
以及**已经跑完的**`liftoff_from_msu7.gff3`(57基因在各自坐标系里的精确位置)。

优势：① 覆盖更全生态型谱系(aus/tropical japonica旱作/indica灌溉/多个地方
传统品种如KYG/LIMA/NABO/PR106/KN/CM/GS/ARC/LM/LX，比5个"教科书代表品种"更
贴近古稻可能来源)；② liftoff已完成，可顺便验证模拟reads比对回irgsp.fa后
坐标是否落回预期位置；③ 已有`.mmi`索引减少环境配置。

⚠️染色体命名各异(如N22是`CM007627.2`格式，非`chr01`)，不影响NGSNGS模拟
本身(只从fasta抽取序列)，只需在后续比对回irgsp.fa后看比对坐标即可。

### NGSNGS二进制定位与环境配置

**踩坑**：NGSNGS是从源码编译的(不是conda装的)，位于
`/home/usr/yinmt/software/NGSNGS/ngsngs`，未加入PATH导致直接敲命令找不到。

```bash
cd /home/usr/yinmt/software/NGSNGS
./ngsngs --version   # v0.9.2.2: 9422b47 (htslib: 1.21)

echo 'export PATH="/home/usr/yinmt/software/NGSNGS:$PATH"' >> ~/.bashrc
source ~/.bashrc
which ngsngs   # 之后可直接用
```

### 关键参数格式确认(来自`ngsngs -h`真实输出)

- `-r`：reads总数(与`-c`覆盖度二选一)
- `-lf`：read长度CDF文件(`length_cdf.txt`格式: 长度\t累积概率)
- `-mf`：**单个**核苷酸置换频率文件(不是分开传5'/3'两个文件)，格式示例：
  ```
  0.865434  0.888339  0.953086  1.000000
  0.882001  0.894563  0.979380  1.000000
  ...
  ```
  ⚠️ 需要把mapDamage2输出的`5pC_to_T_freq.txt`/`3pG_to_A_freq.txt`转换/
  合并成这个格式，具体转换逻辑待补(需要先理解NGSNGS这个格式每列的确切含义)
- `-seq SE`：单端测序(匹配你们的shotgun/capture数据都是单端)
- `-f fq.gz`：输出格式

### ⚠️关键纠正：模板数据必须来自BWA流程，不能沿用Bowtie2

发现`readlen_summary.tsv`和`mapdamage_summary.tsv`两个文件在预期路径下
**都不存在**(`extract_length_dist.sh`/`run_mapdamage_batch.sh`此前从未
真正跑完过，或输出到了未记录的相对路径)。且两个脚本原始`BAM_DIR`指向旧的
Bowtie2版`01.mapping/final`，即使跑过也是错误数据源。

**必须重新指向BWA版数据源跑一遍**：

```bash
# extract_length_dist改用BWA数据
BAM_DIR="/home/scratch/yinmt202607/results/02.irgsp/01.mapping_bwa/final"
OUT_DIR="/home/scratch/yinmt202607/results/02.irgsp/readlen_dist_bwa"
SUMMARY="/home/scratch/yinmt202607/results/02.irgsp/readlen_summary_bwa.tsv"
# (完整脚本见对话记录，逻辑与原脚本一致，仅替换BAM_DIR/OUT_DIR/SUMMARY路径)

# run_mapdamage_batch改用BWA数据
BAM_DIR="/home/scratch/yinmt202607/results/02.irgsp/01.mapping_bwa/final"
OUT_ROOT="/home/scratch/yinmt202607/results/02.irgsp/mapdamage_out_bwa"
SUMMARY="/home/scratch/yinmt202607/results/02.irgsp/mapdamage_summary_bwa.tsv"
```

**原因**：BWA提取的reads集合与Bowtie2不是同一批(BWA多捞回3倍以上reads)，
长度分布不同；损伤统计的样本基数也因reads集合变化而不同，理论上损伤模式
主要由DNA化学降解决定、与比对软件关系不大，但保持全链路口径一致仍需重跑。

### mapDamage2 性能与并行化

**现象**：mapDamage2运行缓慢，卡在"Performing Bayesian estimates"阶段——
这是调用R脚本做MCMC贝叶斯拟合估计损伤率的计算密集步骤，是算法固有特性，
与数据量大小关系不大，非卡死。

**解决：改用sbatch每样本独立提交，并行跑**

```bash
BAM_DIR="/home/scratch/yinmt202607/results/02.irgsp/01.mapping_bwa/final"
REF_FASTA="/home/scratch/yinmt202607/db/asian_rice_panel_index/irgsp.fa"
OUT_ROOT="/home/scratch/yinmt202607/results/02.irgsp/mapdamage_out_bwa"
LOGDIR="/home/scratch/yinmt202607/results/02.irgsp/mapdamage_logs_bwa"
SBATCH_LOGDIR="/home/scratch/yinmt202607/results/02.irgsp/mapdamage_sbatch_logs"
mkdir -p "$OUT_ROOT" "$LOGDIR" "$SBATCH_LOGDIR"

for bam in "${BAM_DIR}"/*.dedup.bam; do
    sample=$(basename "$bam" .dedup.bam)
    outdir="${OUT_ROOT}/${sample}"
    [[ -s "${outdir}/Stats_out_MCMC_correct_prob.csv" ]] && { echo "跳过(已完成): $sample"; continue; }

    sbatch --job-name "mapdmg_${sample}" --cpus-per-task 1 --mem 4G --time 04:00:00 \
        --output "${SBATCH_LOGDIR}/${sample}.%j.out" \
        --wrap "source activate /home/usr/yinmt/.local/mamba/snakemake && mapDamage -i ${bam} -r ${REF_FASTA} -d ${outdir} --merge-reference-sequences > ${LOGDIR}/${sample}.log 2>&1"
done
squeue -u "$USER"
```

参数说明：`--mem 4G`足够(基因组仅约380Mb, BAM文件本身很小)；`--cpus-per-task 1`
(mapDamage2核心统计单线程，多核用不上)；已加断点续跑检查避免重复提交。

---

## 待办清单(截至本次记录)

1. **`-mf`损伤文件格式转换**：需要把mapDamage2的`5pC_to_T_freq.txt`/
   `3pG_to_A_freq.txt`转换成NGSNGS认识的单文件格式，具体转换逻辑待写
2. **重新跑通`readlen_dist_bwa`和`mapdamage_out_bwa`**（sbatch并行进行中）
3. **确认16个基因组染色体命名后**，正式跑NGSNGS生成模拟reads
4. **passport分类表下载**（SNP-Seek）—— PCA路径最大卡点，仍未解决
5. **3K RG的163份DTH8缺失型品种**，需要拿到亚群标签后才能验证是否aus/
   tropical japonica富集
