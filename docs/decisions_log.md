## [结论] 9格提取方法 x 定量比对工具矩阵测试 —— 确认BWA提取是决定性因素

为了验证"坚持用BWA"这个决策是否稳健，做了完整的3x3矩阵测试：3种提取方法
(Bowtie2旧参数/Bowtie2新参数-N1/BWA aln) x 3种定量比对工具(BWA aln/
Bowtie2旧参数/Bowtie2新参数)，共9个组合，每个组合跑完整的16样本
比对→merge→去重→q30流程。详细方案、踩坑记录、完整数据见
`tests/param_matrix_bt2_vs_bwa/`及配套文档
`docs/09_extraction_mapping_matrix_final.md`。

**提取阶段(阶段①)**：3种方法在capture panel1+2上提取出的reads总量排名——
BWA > Bowtie2新参数(-N1) > Bowtie2旧参数。Bowtie2加`-N1`比旧参数有提升，
但离BWA仍有明显差距。

**定量比对阶段(阶段②)总量**：9个组合里`bwa_extract__bwa_map`(现有主线方案)
总q30 reads最高，比表现最差的`bt2_old_extract__bwa_map`(最早历史流程)
高出约1.67倍。

**关键细化(逐样本分析，避免总量堆积造成的误读)**：用16个样本各自的
"单样本最优组合"统计(而非只看总量)，发现**所有拿到"单样本最优"的组合，
提取阶段清一色都是BWA**——不管最后定量比对用BWA还是Bowtie2新参数，两者
互有胜负(q30维度9:7)；但提取阶段一旦用了Bowtie2(不管新旧参数)，基本
没有机会拿到单样本最优。

**结论：真正决定性的变量是提取方法(阶段①)，不是定量比对工具(阶段②)**——
这比笼统的"BWA+BWA最好"更精确：核心是提取阶段必须用BWA，定量比对阶段
用BWA或Bowtie2新参数都可以接受。这是"坚持用BWA"这一决策目前为止最完整、
最经得起检验的证据。

**附加发现**：q30最优与gene_hit最优只有6/16样本一致，说明"reads数量最多"
不完全等价于"基因命中数最多"，与之前"低复杂度序列使reads数量虚高但命中
质量存疑"的结论相互印证。

---

## [决策] 2026-08-16/17：Civán panel C 古代样本 fixed-marker calling —— TV vs ALL 覆盖度对比，瓶颈是古DNA信息量，不是选点bug

**背景**：Stage 50 prototype（LV7008416379, panel C, pooled_mixed, TV, 1000
markers）跑出 `callable_n=0/1000`。排查发现根因是 `50_civan_fixed_marker_
prototype.sh` 的marker选择取的是panel `.snp` 文件顺序的前1000个transversion
位点，全部挤在chr01最前面约0.2Mb，不是全基因组代表性抽样——不是深度问题，
是选点逻辑本身的bug。为此新增了 `19_survey_ancient_coverage.py`（对16个古代
BAM做覆盖度普查，宽松mapq/不卡baseq，输出 `ancient_union_sites.tsv`≥1样本
覆盖 / `ancient_core_sites.tsv`≥N样本覆盖 / `per_sample_coverage_summary.tsv`）
和 `20_filter_coverage_sites_to_transversions.py`（对已有普查结果按panel
REF/ALT做TV后过滤，不重扫BAM）。

**实测结果（16样本，Civán panel全量236万个位点）**：

- ALL track（union=3687个位点，未按TV过滤）：panel_sites_covered分布
  57–703，按`information_flags`(very_low=200/low=500/moderate=2000)分级，
  **14/16样本达到LOW或以上**（2个样本VERY_LOW：LV7008416294=57、
  LV7008416280=193）。
- TV track（union过滤后剩1211个位点，panel_tv_site_n=583812/2365188）：
  panel_sites_covered分布13–260，**15/16样本为VERY_LOW**，只有
  LV6000620166(260)达到LOW，无一样本达到MODERATE。

| sample | ALL covered | TV covered | ALL flag | TV flag |
|---|---|---|---|---|
| LV6000619499 | 222 | 61 | LOW | VERY_LOW |
| LV6000619917 | 251 | 93 | LOW | VERY_LOW |
| LV6000620016 | 286 | 87 | LOW | VERY_LOW |
| LV6000620032 | 513 | 174 | MODERATE | VERY_LOW |
| LV6000620166 | 703 | 260 | MODERATE | LOW |
| LV6000620172 | 456 | 152 | LOW | VERY_LOW |
| LV6000654686 | 210 | 73 | LOW | VERY_LOW |
| LV6000654698 | 385 | 127 | LOW | VERY_LOW |
| LV7008416272 | 202 | 74 | LOW | VERY_LOW |
| LV7008416280 | 193 | 74 | VERY_LOW | VERY_LOW |
| LV7008416294 | 57 | 13 | VERY_LOW | VERY_LOW |
| LV7008416329 | 316 | 105 | LOW | VERY_LOW |
| LV7008416339 | 326 | 114 | LOW | VERY_LOW |
| LV7008416349 | 230 | 83 | LOW | VERY_LOW |
| LV7008416379 | 627 | 182 | MODERATE | VERY_LOW |
| LV7008416407 | 438 | 148 | LOW | VERY_LOW |

（TV列的`panel_sites_covered`直接等于"若Stage 50用这1211个union-TV位点当
marker清单，该样本能call出多少个"——因为union定义就是"至少1个样本覆盖"，
某样本自己覆盖的TV位点必然落在union集合里，无需另算。）

**结论**：

1. Marker选点bug修好之后，瓶颈从"选点逻辑错误"变成"这批样本在
   damage-safe(TV-only)track下的真实信息量不足"——ALL/TV两组数字的巨大
   落差（14/16 LOW+ vs 15/16 VERY_LOW）证明这不是工程问题，是这批古DNA样本
   相对Civán panel密度的真实覆盖度问题。
2. **TV继续作为主分析轨道，不因为位点数量少而放宽**——不能靠转换成ALL
   track "凑数"来让样本看起来callable。
3. **ALL可以作为damage-sensitive的敏感性分析轨道保留，但不能冒充主结论**。
4. **`information_flags`阈值(very_low=200/low=500/moderate=2000)不修改**。
5. **Stage 60–80 不因为ALL track的数字看起来更好看而强行解锁**——解锁条件
   仍按原设计走。

No parameter change was made. No capture/shotgun separation was inferred.
capture_bait_bed remains null.

**实现更新（2026-08-17）**：coverage-aware marker universe已接入Stage 50，
ALL/TV分别输出，所有16个样本均要求完成技术calling；`technical_execution`
与基于`callable_n`的`scientific_projection`分开记录。分级固定为
`>=200 formal_validation_candidate`、`50–199 exploratory_projection`、
`<50 descriptive_only`。TV仍为主轨道，ALL仅为damage-sensitive敏感性轨道；
没有修改MAPQ/BaseQ、MAF或信息量阈值，也没有解锁Stage 60。

相关文件：`scripts/ecotype_pca_v2/19_survey_ancient_coverage.py`、
`scripts/ecotype_pca_v2/20_filter_coverage_sites_to_transversions.py`、
commit `29418bb`（普查脚本）、`66ec621`（TV后过滤脚本）。
