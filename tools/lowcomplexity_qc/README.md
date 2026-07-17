# lowcomplexity_qc.sh

通用工具：检测"基因命中"类分析里，有多少条reads落在参考基因组的低复杂度/重复
序列区域内。这类区域(如GC重复的脯氨酸/丙氨酸编码密码子串联，常见于MADS-box、
NF-Y、CCT-domain等转录因子基因编码区)容易导致短read比对算法(BWA aln等)产生
"看似唯一、实则可比对到基因组多处"的假阳性唯一比对，不能直接当作可信的功能
位点证据使用。

## 来源

这套流程是在57个开花基因的BWA/Bowtie2命中数据分析中，逐步手动排查出来的
(详见 `../../docs/flank1kb_msa_exploration.md`)。关键发现：直接用BWA的
"唯一比对"判定(SAM tag `X0:i:1`)不可靠，一条被BWA判为"唯一比对"到DTH7/OsPRR37
的42bp read，BLAST交叉验证后发现能以95-100%相似度匹配到chr02上几十个不同位置。
批量筛查发现57基因命中数据中约35.7%的信号落在低复杂度区、可信度存疑。

现在把这套手动排查流程固化成可复用工具，任何新的"基因命中"数据集，只要格式
兼容，换个文件路径参数就能直接跑，不需要重新摸索。

## 依赖

- `dustmasker` (NCBI BLAST+ 套件自带)
- `makeblastdb` / `blastn` (NCBI BLAST+ 套件自带)
- `bedtools`
- `samtools` (仅 blast-spotcheck 子命令需要)

## 三个子命令

### 1. build-mask（每个参考基因组只需要跑一次）

```bash
bash lowcomplexity_qc.sh build-mask \
    --ref /home/scratch/yinmt202607/db/asian_rice_panel_index/irgsp.fa \
    --out-dir /home/scratch/yinmt202607/db/asian_rice_panel_index/lowcomplexity_qc
```

生成：
- `<out-dir>/lowcomplexity.sorted.bed` —— 标准BED格式的低复杂度区域清单
- `<out-dir>/blastdb/ref.*` —— BLAST索引，供第3步交叉验证用

### 2. check-hits（批量标注一份命中数据）

```bash
bash lowcomplexity_qc.sh check-hits \
    --hits /home/scratch/yinmt202607/results/02.irgsp/02.gene_hits/gene_hits_with_metadata.tsv \
    --mask /home/scratch/yinmt202607/db/asian_rice_panel_index/lowcomplexity_qc/lowcomplexity.sorted.bed \
    --out /home/scratch/yinmt202607/results/02.irgsp/gene_hits_lowcomplexity_check.tsv
```

输入TSV只要求：带表头，包含染色体列、比对起始位置列(1-based)、read长度列，
具体是第几列用 `--chr-col` / `--pos-col` / `--len-col` 指定（不需要固定列顺序，
换一份列结构不同的新数据集，改几个数字参数就行，不用改脚本本身）。

默认参数对应`gene_hits_with_metadata.tsv`的列结构（chr=第7列, pos=第8列,
read_length=第13列, 标签=第1列+第6列即sample|gene）。

输出会在终端打印汇总统计（总命中数/落在低复杂度区的数量/占比/具体清单），
同时生成完整的标注TSV文件。

### 3. blast-spotcheck（针对单条read的人工抽查/交叉验证）

```bash
bash lowcomplexity_qc.sh blast-spotcheck \
    --bam /home/scratch/yinmt202607/results/02.irgsp/05.flank1kb/bam/LV7008416339.flank1kb.bam \
    --region chr07:29616768-29616789 \
    --blastdb /home/scratch/yinmt202607/db/asian_rice_panel_index/lowcomplexity_qc/blastdb/ref
```

用于人工抽查某条具体read"BWA判定的唯一比对"是否真的可信——如果BLAST输出
只有1行结果，说明可信；如果输出多行且bitscore/pident相近，说明这条read
实际上能匹配基因组多个位置，BWA的判定是假阳性。

## 使用建议

任何新一批"基因命中"分析，建议在解读结果前，先跑一遍 `check-hits`，把落在
低复杂度区域的命中单独标记出来，报告时明确区分"干净证据"和"存疑证据"，
不要把两者混在一起下结论。
