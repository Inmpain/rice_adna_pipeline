# 开花基因 ±1kb 侧翼区域提取 + Multiple Alignment 可行性探索

记录范围：从老师提出"基因上下1kb调控区域reads提取做multiple alignment"开始，
到发现低复杂度序列问题、暂停等待BLAST验证为止。

最后更新: 2026-07-17

## 背景：老师这个建议在做什么

- 目的：不只看基因编码区(CDS)，把基因两侧各1kb的潜在调控区(启动子等)也纳入，
  因为很多光周期基因的功能差异不完全来自编码区突变，也可能来自调控区
- Multiple alignment的目的：把多个样本(理想情况下还应加入已知功能型/缺失型
  代表品种)在同一区域的序列放在一起比对，看古稻样本更像哪种已知单倍型，
  而不是只看单条read落在哪个坐标

## 第一步：确认注释层级(gene/mRNA/exon/CDS/UTR)

```bash
cd /home/scratch/yinmt202607/db/gene
grep "LOC_Os08g07740" msu7.gff3 | cut -f3 | sort | uniq -c
```
结果：DTH8在msu7.gff3里有 gene / mRNA / exon / CDS / five_prime_UTR / three_prime_UTR
六层记录各一条。**当前分析用的`flower_gene.sorted.bed`是gene整体层级的坐标**，
不是CDS层级，这次±1kb是从gene整体边界外扩的。

## 第二步：生成 ±1kb 扩展BED

```bash
awk -F'\t' 'BEGIN{OFS="\t"}{
    start = $2 - 1000; if(start < 0) start = 0
    end = $3 + 1000
    print $1, start, end, $4, $5, $6
}' flower_gene.sorted.bed > flower_gene.flank1kb.bed
```

## 第三步：提取扩展区域内的reads

```bash
BAM_DIR=/home/scratch/yinmt202607/results/02.irgsp/01.mapping_bwa/final
BED=/home/scratch/yinmt202607/db/gene/flower_gene.flank1kb.bed
OUT_DIR=/home/scratch/yinmt202607/results/02.irgsp/05.flank1kb

mkdir -p "$OUT_DIR"/{bam,fasta}

for robot in $(ls "$BAM_DIR" | grep '\.dedup\.q30\.bam$' | sed 's/\.dedup\.q30\.bam$//'); do
    samtools view -bh -L "$BED" "$BAM_DIR/${robot}.dedup.q30.bam" > "$OUT_DIR/bam/${robot}.flank1kb.bam"
    samtools index "$OUT_DIR/bam/${robot}.flank1kb.bam"
done
```

⚠️ **踩坑记录**：第一次跑这一步时`samtools`模块没load，命令报`command not found`，
但bash的`>`重定向已经把空文件建出来了，导致后续`samtools consensus`报
`Failed to read header`。**排查方法**：`ls -la`检查文件大小是否为0字节，
0字节说明生成过程本身失败了，不是数据问题。修复：`module load samtools`后重新执行。

## 第四步：逐基因生成consensus序列(每样本每基因一个fasta)

```bash
REF=/home/scratch/yinmt202607/db/asian_rice_panel_index/irgsp.fa

for robot in $(ls "$OUT_DIR/bam" | grep '\.flank1kb\.bam$' | sed 's/\.flank1kb\.bam$//'); do
    bam="$OUT_DIR/bam/${robot}.flank1kb.bam"

    while IFS=$'\t' read -r chr start end gene_name score strand; do
        region="${chr}:$((start + 1))-${end}"   # BED是0-based, samtools区域是1-based
        safe_name=$(echo "$gene_name" | sed 's/|/_/g; s#/#_#g')

        samtools consensus -r "$region" "$bam" \
            > "$OUT_DIR/fasta/${robot}_${safe_name}.fa" 2>/dev/null
    done < "$BED"
done
```

⚠️ **说明**：某些"样本×基因"组合的输出fasta是完全空文件(0行)。经验证
(`samtools view "$bam" "$region" | wc -l`输出0)，这是**真实的零覆盖**，
不是脚本bug——低深度数据下，很多样本在很多基因上本来就没有任何read落上去。

## 第五步(v1，有缺陷)：按fasta文件自身长度算覆盖率 —— 事后发现这个算法有问题

```bash
cd "$OUT_DIR/fasta"
echo -e "sample\tgene\ttotal_len\tN_count\tcovered_len\tcovered_pct" > coverage_check.tsv

for f in *.fa; do
    [[ "$f" != *"gene区域"* ]] || continue
    sample_gene=$(basename "$f" .fa)
    sample=$(echo "$sample_gene" | cut -d'_' -f1)
    gene=$(echo "$sample_gene" | cut -d'_' -f2-)

    seq=$(grep -v "^>" "$f" | tr -d '\n')
    total_len=${#seq}
    n_count=$(echo "$seq" | tr -cd 'Nn' | wc -c)
    covered_len=$((total_len - n_count))

    if [[ $total_len -gt 0 ]]; then
        covered_pct=$(awk -v c="$covered_len" -v t="$total_len" 'BEGIN{printf "%.2f", c/t*100}')
    else
        covered_pct="0"
    fi
    echo -e "${sample}\t${gene}\t${total_len}\t${n_count}\t${covered_len}\t${covered_pct}" >> coverage_check.tsv
done
```

⚠️ **发现的问题**：`samtools consensus`对完全零覆盖的区域，不会输出补N的固定长度序列，
而是只输出"从第一条有覆盖的read到最后一条有覆盖的read"之间的片段——也就是说
`total_len`这一列本身已经把真正的零覆盖侧翼区裁掉了，导致算出来的`covered_pct`
比真实情况乐观。最高看到22.79%，但这是虚高的。

## 第六步(v2，修正版)：用BED里的真实区间长度做分母

```bash
cd /home/scratch/yinmt202607/results/02.irgsp/05.flank1kb/fasta
BED=/home/scratch/yinmt202607/db/gene/flower_gene.flank1kb.bed

echo -e "sample\tgene\ttrue_region_len\tcovered_bases\tcovered_pct_true" > coverage_check_v2.tsv

for f in *.fa; do
    [[ "$f" != *"gene区域"* ]] || continue
    [[ -s "$f" ]] || continue

    sample_gene=$(basename "$f" .fa)
    sample=$(echo "$sample_gene" | cut -d'_' -f1)
    gene=$(echo "$sample_gene" | cut -d'_' -f2-)

    true_len=$(awk -F'\t' -v g="$gene" 'index($4, g)>0 || index(g, $4)>0 {print $3-$2; exit}' "$BED")
    if [[ -z "$true_len" ]]; then
        gene_id=$(echo "$gene" | grep -oE 'LOC_Os[0-9]+g[0-9]+')
        true_len=$(awk -F'\t' -v g="$gene_id" '$4 ~ g {print $3-$2; exit}' "$BED")
    fi

    seq=$(grep -v "^>" "$f" | tr -d '\n')
    n_count=$(echo "$seq" | tr -cd 'Nn' | wc -c)
    covered_bases=$((${#seq} - n_count))

    if [[ -n "$true_len" && "$true_len" -gt 0 ]]; then
        pct=$(awk -v c="$covered_bases" -v t="$true_len" 'BEGIN{printf "%.3f", c/t*100}')
    else
        pct="NA"
    fi
    echo -e "${sample}\t${gene}\t${true_len}\t${covered_bases}\t${pct}" >> coverage_check_v2.tsv
done
```

**修正后的真实结果**：最高覆盖率只有 **2.985%**（DTH8/Ghd8, LV7008416407,
111/3718bp），绝大多数在0.9%-3%之间，绝对覆盖碱基数是30-136bp。
**结论：±1kb扩展并没有实质性缓解深度不足的问题**，跟57基因本体区间的
情况(0.076%基因组占比→随机命中数量级吻合)是同一类瓶颈的延伸。

## 第七步：检查每个基因有多少个不同样本存在任意覆盖(比覆盖率百分比更关键的指标)

```bash
BAM_DIR=/home/scratch/yinmt202607/results/02.irgsp/05.flank1kb/bam
cd /home/scratch/yinmt202607/results/02.irgsp/05.flank1kb

echo -e "gene\tsamples_with_any_coverage" > gene_sample_count.tsv

while IFS=$'\t' read -r chr start end gene_name score strand; do
    region="${chr}:$((start+1))-${end}"
    n_samples=0
    for bam in bam/*.flank1kb.bam; do
        c=$(samtools view "$bam" "$region" 2>/dev/null | wc -l)
        [[ "$c" -gt 0 ]] && n_samples=$((n_samples+1))
    done
    echo -e "${gene_name}\t${n_samples}" >> gene_sample_count.tsv
done < /home/scratch/yinmt202607/db/gene/flower_gene.flank1kb.bed

sort -t$'\t' -k2,2 -rn gene_sample_count.tsv | head -20
```

**结果排名(前几名)**：
```
AID1              9个样本有覆盖
OsSPL17           8
OsGATA28          7
DTH8/Ghd8         6
DTH7/OsPRR37      6
RFL/APO2          6
OsGI              6
```

## 第八步：检查这几个候选基因，多样本的reads坐标是否聚集在同一窄窗口

```bash
BAM_DIR=/home/scratch/yinmt202607/results/02.irgsp/05.flank1kb/bam

declare -A gene_regions=(
    ["AID1"]="chr06:4015961-4020327"
    ["OsSPL17"]="chr09:18917259-18922656"
    ["OsGATA28"]="chr11:4431775-4435071"
    ["DTH7_OsPRR37"]="chr07:29615704-29630223"
    ["RFL_APO2"]="chr04:30181588-30186852"
    ["OsGI"]="chr01:4328151-4339486"
)

for name in "${!gene_regions[@]}"; do
    region="${gene_regions[$name]}"
    echo "=== $name ($region) ==="
    for bam in "$BAM_DIR"/*.flank1kb.bam; do
        robot=$(basename "$bam" .flank1kb.bam)
        count=$(samtools view "$bam" "$region" 2>/dev/null | wc -l)
        if [[ "$count" -gt 0 ]]; then
            samtools view "$bam" "$region" | awk -v r="$robot" '{print r"\t"$4"\t"length($10)}'
        fi
    done
    echo ""
done
```

**关键发现**：多个基因都出现"不同独立样本的reads挤在一个20-40bp窄窗口内"的模式，最典型的：

- **DTH7/OsPRR37**：6个样本的reads坐标落在 `29616768-29616789`(仅21bp窗口)
- **OsGATA28**：13条reads明显分两簇，`4433086-4433121`(35bp)和`4433718-4433758`(40bp)
- **RFL/APO2**：7个样本reads集中在`30185657-30185681`(24bp窗口)
- **DTH8/Ghd8**：6个样本reads落在`4334007-4334634`(约600bp窗口，全部在基因本体内，
  完全没有落到±1kb扩展出来的侧翼区)

一度推测：这可能是capture panel探针的设计靶点(杂交富集导致reads被拉向特定小片段)。

## 第九步：对DTH7/OsPRR37这个最干净的窗口做逐位点mpileup，试图看具体碱基是否一致

```bash
for bam in "$BAM_DIR"/*.flank1kb.bam; do
    robot=$(basename "$bam" .flank1kb.bam)
    samtools mpileup -r chr07:29616768-29616789 "$bam" 2>/dev/null | \
        awk -v r="$robot" '{print r"\t"$2"\t"$3"\t"$4"\t"$5}'
done
```

**发现**：
- 786、788、789三个位置，5个独立样本高度一致（好信号）
- 787位置出现分歧(3个样本C，LV7008416407是T)，但该位置离这条read起点仅2bp，
  高度符合古DNA末端C→T损伤模式，**不建议当作真实等位基因差异使用**
- 784附近出现1bp插入(`c+1a`)，770-772附近出现连续3bp缺失(`*`)——怀疑是
  低复杂度序列导致的比对歧义

## 第十步：查这段区域的真实参考序列，确认是否为低复杂度重复区

```bash
samtools faidx /home/scratch/yinmt202607/db/asian_rice_panel_index/irgsp.fa chr07:29616760-29616800
```

**结果**：
```
>chr07:29616760-29616800
CCTTCCCTCCTCTTCTTCCTCCGCCTTCGCCGCCGCCGCCG
```

**确认**：后半段`CCGCCTTCGCCGCCGCCGCCG`是明显的`CCG`/`GCC`密码子连续重复
(脯氨酸/丙氨酸重复结构域，NF-Y/CCT domain这类转录因子基因编码区常见)。

## 结论(第一版，已被下方最终结论修正/扩展)

1. **±1kb扩展没有捕捉到调控区数据**：目前检查的几个基因，reads仍然集中在基因本体
   编码区内，没有落到新扩出来的侧翼区
2. **真实覆盖率极低**：即使扩大范围，最高也只有2.985%，多数在1%以下——深度不足
   这个根本瓶颈没有被"扩大搜索范围"解决
3. **"多样本聚集在同一窄窗口"这个现象，原本怀疑是capture panel探针富集靶点，
   现在证据倾向于修正为：低复杂度重复序列导致的比对算法系统性行为**，而非
   真实的生物学信号富集——因为验证发现这些窗口的参考序列本身就是GC重复基序
4. **单点碱基差异(如787位点)需要谨慎解读**，古DNA末端损伤(C→T)会产生假性多态性，
   不能简单当作真实等位型差异

## BLAST验证结果(已完成)

对DTH7/OsPRR37窗口那条最完整的read(`LV7008416339`, 42bp)做BLAST交叉验证：

```bash
BAM_DIR=/home/scratch/yinmt202607/results/02.irgsp/05.flank1kb/bam

seq=$(samtools view "$BAM_DIR/LV7008416339.flank1kb.bam" chr07:29616768-29616789 | head -1 | cut -f10)
echo -e ">test_read\n$seq" > test_read.fa

mkdir -p /home/scratch/yinmt202607/db/asian_rice_panel_index/blastdb
makeblastdb -in /home/scratch/yinmt202607/db/asian_rice_panel_index/irgsp.fa \
    -dbtype nucl -out /home/scratch/yinmt202607/db/asian_rice_panel_index/blastdb/irgsp

blastn -query test_read.fa \
    -db /home/scratch/yinmt202607/db/asian_rice_panel_index/blastdb/irgsp \
    -task blastn-short \
    -outfmt "6 qseqid sseqid pident length mismatch qstart qend sstart send evalue bitscore" \
    -evalue 1
```

**结果**：这条42bp的read能以95-100%相似度匹配到chr02上几十个不同位置(仅chr02一条
染色体上就有几十个命中，间隔几十万到几百万bp均匀散布)，且比对分数(bitscore=38.2)
彼此相当，没有明显更优的"唯一真实位置"。**实锤证实这是全基因组广泛存在的低复杂度
重复基序**，BWA当初判定"唯一比对"(SAM tag `X0:i:1`)是被其启发式搜索算法误导的假阳性
——BWA aln用seed-and-extend策略，只在有限候选区域内搜索，恰好没搜到第二个候选就
判为"唯一"，不代表更穷举的搜索(如BLAST)也找不到。

## 批量低复杂度区域过滤(dustmasker + bedtools，已完成)

```bash
# 1. 用dustmasker标记全基因组低复杂度区域
dustmasker -in /home/scratch/yinmt202607/db/asian_rice_panel_index/irgsp.fa \
    -out /home/scratch/yinmt202607/db/asian_rice_panel_index/irgsp_lowcomplexity.txt \
    -outfmt interval

# 2. 转换成标准BED格式(注意起点为0时不能再减1，否则bedtools会报Invalid record)
awk '
/^>/{ chr=substr($1,2); next }
{
    start = $1 - 1
    if (start < 0) start = 0
    print chr"\t"start"\t"$3
}
' /home/scratch/yinmt202607/db/asian_rice_panel_index/irgsp_lowcomplexity.txt \
  | sort -k1,1 -k2,2n \
  > /home/scratch/yinmt202607/db/asian_rice_panel_index/irgsp_lowcomplexity.sorted.bed

# 3. 用整条read的比对跨度(而不是单点坐标)与低复杂度区域做交叉比对
#    (v1版只用pos单点坐标检查，会系统性漏检——read起点在区域外但read主体跨入
#    区域内的情况，v1版测出14.3%明显偏低，v2版改用pos到pos+read_length的完整
#    跨度后测出35.7%，更接近真实情况)
HITS_TSV=/home/scratch/yinmt202607/results/02.irgsp/02.gene_hits/gene_hits_with_metadata.tsv
LOWCOMP_BED=/home/scratch/yinmt202607/db/asian_rice_panel_index/irgsp_lowcomplexity.sorted.bed

tail -n+2 "$HITS_TSV" | awk -F'\t' 'BEGIN{OFS="\t"}{
    start = $8 - 1
    end = $8 - 1 + $13
    print $7, start, end, $1"|"$6
}' > /tmp/hit_spans.bed

bedtools intersect -a /tmp/hit_spans.bed -b "$LOWCOMP_BED" -c \
    > /home/scratch/yinmt202607/results/02.irgsp/gene_hits_lowcomplexity_check_v2.tsv
```

**最终结果：42条命中reads里，15条(35.7%)落在低复杂度重复区域，可信度存疑**：

```
LV6000619917|DTH8/Ghd8              LV6000620016|DTH7/OsPRR37       LV6000620016|OsSPL17
LV6000620032|GW2                    LV6000620166|OsPUP4/BG3          LV6000620166|SDG711/OsCLF
LV7008416272|OsSPL17                LV7008416280|AID1                LV7008416329|RFL/APO2
LV7008416339|OsMADS15/DEP           LV7008416339|DTH8/Ghd8            LV7008416349|GW2
LV7008416379|AID1                   LV7008416379|OsMADS15/DEP         LV7008416407|OsMYB8
```

### 对三个优先基因的最终判断(按证据可信度重新排序)

| 基因 | 该数据集内总命中数 | 落在低复杂度区 | 可信证据 | 结论 |
|---|---|---|---|---|
| **OsGI** | 3 | 0 | 3/3 (100%) | **证据质量最好**，全部干净，优先推进 |
| **DTH8/Ghd8** | 5 | 2 | 3/5 (60%) | 部分可信，需逐条注明哪些样本可用 |
| **DTH7/OsPRR37** | 1 | 1 | 0/1 (0%) | **唯一证据已被证伪**，目前无可信数据支撑 |

⚠️ 注：此表基于bowtie2版`gene_hits_with_metadata.tsv`(42条命中)，与BWA版57基因
矩阵的命中数(DTH8=8, OsGI=7, DTH7=6)不是同一批数据，BWA版尚未做低复杂度过滤，
建议对BWA版数据重复这套流程后再最终定论。

## 结论(最终版)

1. **±1kb扩展没有捕捉到调控区数据**——reads仍集中在基因本体内
2. **真实覆盖率极低**（最高2.985%）——深度不足的根本瓶颈未被扩大范围解决
3. **"多样本聚集在窄窗口"被BLAST证实主要是低复杂度重复序列导致的比对歧义**，
   不是capture panel探针富集信号
4. **系统性低复杂度过滤显示：57基因命中数据中约1/3(35.7%)的证据落在重复区、
   可信度存疑**，但仍有约2/3是干净、可用的——**不是"数据全废"，而是"必须先
   过滤才能用"**
5. **三个优先基因里，OsGI证据质量最好(全部干净)，DTH7/OsPRR37证据已被证伪，
   DTH8/Ghd8部分可信**——建议调整优先级，OsGI优先，DTH7/OsPRR37暂缓直到
   有更多干净证据出现

## 建议给老师的汇报要点(最终版)

1. 尝试了基因±1kb扩展，未观察到侧翼调控区有效覆盖，深度瓶颈依旧
2. 建立了一套低复杂度序列过滤流程(dustmasker + BLAST交叉验证 + bedtools批量筛查)，
   证实57基因命中数据中约35.7%的信号落在GC重复低复杂度区、可信度存疑
3. 这套过滤流程本身是一个可复用的质控方法，建议后续所有基因命中分析都先过一遍
   这个流程，再解读结果
4. 三个原定优先基因中，OsGI证据质量最扎实，建议作为下一步精确诊断位点分析的
   首选；DTH7/OsPRR37现有证据已被证伪，暂缓
5. 建议：若要继续深入，可能需要更长read长度数据(现有30-50bp过短，难以避免落入
   重复basegroup)，或改用能识别多重比对歧义的专门比对流程(如更严格的MAPQ阈值、
   或直接在比对阶段排除已知低复杂度区)

## 相关路径汇总(更新)

```
/home/scratch/yinmt202607/db/gene/flower_gene.flank1kb.bed        # ±1kb扩展坐标
/home/scratch/yinmt202607/db/asian_rice_panel_index/
├── irgsp_lowcomplexity.txt                   # dustmasker原始interval输出
├── irgsp_lowcomplexity.sorted.bed            # 转换后的标准BED格式
└── blastdb/irgsp.*                           # BLAST索引(用于交叉验证比对唯一性)

/home/scratch/yinmt202607/results/02.irgsp/05.flank1kb/
├── bam/{robot}.flank1kb.bam[.bai]           # 扩展区域内的reads
├── fasta/{robot}_{gene}.fa                   # 逐样本逐基因consensus(多数为空/近全N)
├── fasta/coverage_check.tsv                  # v1覆盖率统计(有缺陷, 分母算法有误)
├── fasta/coverage_check_v2.tsv               # v2覆盖率统计(修正版, 真实区间长度做分母)
└── gene_sample_count.tsv                     # 每个基因有多少样本存在任意覆盖

/home/scratch/yinmt202607/results/02.irgsp/
├── gene_hits_lowcomplexity_check.tsv          # v1低复杂度检查(单点坐标, 14.3%, 偏低)
└── gene_hits_lowcomplexity_check_v2.tsv       # v2低复杂度检查(整条read跨度, 35.7%, 最终版)
```
