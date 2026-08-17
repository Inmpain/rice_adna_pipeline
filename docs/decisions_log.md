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

---

## [发现] 2026-08-17：Civán panel `civan_snp.snp` 的REF/ALT两列相对IRGSP参考基因组被系统性标反，但不影响已有PCA结果

**背景**：新增的 `23_validate_snp_ref_against_fasta.py`（搬自同实验室师兄
`Snakefile.pseudohaploid.from_panel`的REF校验逻辑）首次对`civan_snp.snp`全量
236万个位点跑，报`2365188/2365188`（100%）REF与IRGSP参考FASTA（
`/home/scratch/yinmt202607/db/asian_rice_panel_index/irgsp.fa`）不匹配。

**诊断**：逐条核对报告中的示例与`civan_snp.snp`自己的REF/ALT列，发现**每一条
例子的FASTA真实碱基都精确等于声明的ALT列，从不是其他碱基**（例：位点
1249，`civan_snp.snp`声明REF=C/ALT=A，FASTA chr01:1249实际为A）。236万个
位点如果是坐标系/基因组版本对不上，应该是雜乱无规律地两边都对不上（纯随机
大概仅25%-33%命中率），而不是100%精确地总是对上ALT这一列——这是**列递归
头被整体标反**的签名，不是数据损坏、也不是基因组版本/坐标不一致（位置本身是
对的，只是"哪个算REF、哪个算ALT"这个标签方向反了）。

**为什么不会污染已有结果**：

1. 现代595个axis builder的`.eigenstratgeno`矩阵是与`civan_snp.snp`同一次
   VCF→plink2→convertf转换的兄弟文件，0/1/2编码本来就是"相对`civan_snp.snp`
   自己声明的REF/ALT标签"来算的，与真实基因组无关，内部自洽（Stage 20那次
   PC1+2=27.79%精确复现论文数字，就是这个自洽性的证据）。
2. `10_call_ancient_fixed_markers.py`给古代样本call基因型时，拿BAM里的碱基
   去比对`civan_snp.snp`自己的`record["ref"]`/`record["alt"]`字符串——**不是**
   去比对真实FASTA。只要古代和现代两边用的是同一份`civan_snp.snp`的标签空间，
   0/2编码就是互相对得上的，PCA投影的数学本身没问题。
3. TV/transition判定用`frozenset((ref,alt))`，与哪个算ref无关，不受影响。
4. 唯一有实际影响的是**解释层面**：`decisions.tsv`里`CALLED_REF`/`CALLED_ALT`
   这两个状态名，如果拿去对应真实生物学意义（比如"这个古代样本在N个位点
   携带参考等位基因"），方向是反的——但这只影响措辞解读，不影响PCA坐标本身
   对不对。

**处理方案**：`23_validate_snp_ref_against_fasta.py`改成能区分"干净总体反
标"与"真的乱了"两种情况：每个位点分类match_ref/match_alt/true_mismatch/
no_such_contig/out_of_range，**只有**全部干净匹配REF（普通情况）或全部干净
匹配ALT且true_mismatch/no_such_contig/out_of_range均为0（确认反标）才会PASS；
**任何**match_ref与match_alt混合、或出现任一true_mismatch/no_such_contig/
out_of_range，仍硬阻FATAL——保证这道检查仍能拓住真正损坏或错基因组版本的panel
（那种情况会表现为散乱/不一致的不匹配，而不是这种干净的二分划分），不会因为
这次发现变成形同虚设。

**不做的事**：不反过来修改流程里实际的REF/ALT标签（比如交30写fixed_reference
.snp时交换REF/ALT列、同步翻转现代矩阵的0/2编码）——既然PCA数学不受影响，
这个更大面积的"修正"没必要，不在这次范围内做。不存在任何数值层面的risk。

No parameter change was made. No capture/shotgun separation was inferred.
capture_bait_bed remains null. MAPQ/BaseQ/MAF/LD参数未受影响。

相关文件：`scripts/ecotype_pca_v2/23_validate_snp_ref_against_fasta.py`。

---

## [决策] 2026-08-17：Civán panel C shared marker 集从"全景 LD 剪枝再交集"改成"ancient 覆盖度优先 LD 剪枝"

**背景**：51 号 runner（`51_civan_maf_ld_and_private_axis.sh`）第一次实跑（primary
敏感度，ALL track）得到的共享 marker 集只有 **47 个**（526,936 geno/MAF 过滤后
→ LD 剪枝(100kb/r²=0.2)剩 17,708 → 与 ancient union coverage(3,687 个位点)取
交集剩 47）。逐样本 callability 检查显示 16 个古样本全部落在 `VERY_LOW`（最高
4/47，多个样本 0/47），smartpca 的 `lsqproject` 因此把大量样本（包括不在 poplist
里、走投影路线的现代 `O._rufipogon`）判定为 `insufficient data` 直接从 evec 里
剔除，触发 `15_pca_qc.py` 的硬性 FATAL 门槛。

**排查过程**：先尝试把 sensitivity 从 `primary`(100kb) 换成 `S1`(50kb)，
但分析后判断瓶颈不在 LD 窗口大小——LD 剪枝是在全基因组 526,936 个位点上独立于
ancient 覆盖度进行的，选中的"每个 LD block 代表位点"和 ancient shotgun
reads 实际落点几乎是两个不相关的稀疏集合，交集天然很小，光调窗口大小无法
根本解决。

**结论与处理**：把执行顺序反过来——不再是"先在全景上做 LD 剪枝，再祈祷剩下
的位点恰好被 ancient 覆盖"，而是"先把 geno/MAF 过滤后的位点与 ancient union
coverage 取交集，得到 ancient 确实覆盖到的候选位点，再只在这个小范围内做 LD
剪枝"，保证每一个通过 LD 剪枝的位点都保证是 ancient 覆盖过的。新增
`scripts/ecotype_pca_v2/27_ancient_coverage_first_ld_prune.py`，输入是 07
`--stage geno_maf_only` 的产出加上一份候选 snplist（`25_intersect_snplists.py`
算出的 geno/MAF 位点 ∩ ancient coverage snplist），只在候选集合内部跑
`plink2 --indep-pairwise`（ld_window_kb/ld_r2 仍从 config 按 panel+sensitivity
解析，未改动任何冻结数值）。`51_civan_maf_ld_and_private_axis.sh` 的 Track
循环相应改成：07(geno_maf_only) → 21+25(候选交集) → 27(候选内部 LD 剪枝)，
产出文件名 `civan.<TRACK>.<SENSITIVITY>.ancient_first.fixed.snplist`。

**不做的事**：没有改动 `07_make_fixed_markers.sh` 本身——它是 A/B/C 三个面板
共用的冻结脚本，Panel A(3K)/Panel B(720) 目前都还没有实际跑过这套 marker
选择流程（`results/ecotype_pca_v2/` 在仓库里不存在，两者对应的 workflow 阶段
`60_panel_a_3k_prototype`/`70_panel_b_720_decision_and_audit` 都还是 locked），
所以这次改动只新增了一条可选路径，不影响它们将来沿用原有"全景 LD 剪枝→交集"
顺序（如果那样更合适的话）。MAF/geno/LD 数值本身未改动，只改了应用顺序。

相关文件：`scripts/ecotype_pca_v2/27_ancient_coverage_first_ld_prune.py`、
`scripts/ecotype_pca_v2/workflow/runners/51_civan_maf_ld_and_private_axis.sh`。
