# 旱稻/水稻生态型 PCA 分析——数据库与设计草案

> 本文档记录一条与 `codex/oryza-competitive-mapping`(besthit)平行、下游依赖
> 其输出的新分析线：把 besthit 过滤后确认为 Oryza 的古代 reads，投影到现代
> 群体 PCA 空间，判断每个古代样本更接近旱稻还是水稻生态型。
> 首次写下: 2026-08-05

---

## 📍📍 2026-08-15 新增：`scripts/ecotype_pca_v2/` ——冻结统计设计的独立新流水线

**跟本文档下面所有内容（`scripts/ecotype_pca/`，下称v1）是并行的两套代码，
不是替换关系**：

- **v1（`scripts/ecotype_pca/`）**：本文档主体描述的流水线。唯一在服务器
  真实数据上跑成功过的是Civáň panel、LV7008416379的leave-one-out冒烟测试
  （见📍原有section第2/3条），以及poplist bug修复后2026-08-15的重新验证
  （见📍原有section第3条，PC1/PC2从0.0283/-0.0074变到0.0393/-0.0103，排序
  结论aromatic仍最近未变——**这条验证记录不是本次Claude Code会话做的**，
  应该来自并行会话或用户直接操作，本次会话开始时被明确告知"除了最初那次
  bug版本的冒烟测试，其余全部没在服务器跑过"，与本节内容有出入，接手时
  以`squeue`/`ls`实测结果为准，不要假设任一版本的说法绝对准确）。
  `pseudo_haploid_call.py`（v1的伪单倍型调用脚本）本次会话修了两轮真实
  bug（见下）。
- **v2（`scripts/ecotype_pca_v2/`）**：本次会话新开的独立目录，用户给了一份
  "STATISTICAL DESIGN IS FROZEN"的完整冻结规范（`docs/ECOTYPE_PCA_V2_SPEC.md`），
  核心是修正v1的一个设计缺陷——v1的`build_sample_panel_subset.py`每个古样本
  各自把panel缩到自己覆盖的位点再建轴，导致不同古样本的PC1/PC2其实是不同
  marker set上算出来的，不能跨样本比较；v2改成reference-first：**先用现代
  参考样本冻结MAF/LD-pruned的marker set和PCA轴，古样本只做投影，永远不参与
  建轴**。v2目前**完全没有在服务器真实数据上跑过**——只在本机构造的合成
  EIGENSTRAT/BAM数据上跑通过33+18个单元/集成测试（pure-Python部分实际跑过，
  `test_integration_synthetic.sh`需要plink2，本机没有，没跑）。

**v2当前状态（3次commit，`codex/ecotype-pca-panel`分支）**：
1. `10878d7`——Batch 1初版：`00`-`08`共9个脚本+`lib_ecotype_v2.py`+config。
2. `03032ed`——Batch 1 correction：修复GPT review发现的15项实现问题（enum
   硬校验、04的LD计算改成分块流式避免OOM、07/08接口修复、manifest字段拆分、
   overwrite保护等），新增33个实测通过的单元测试。
3. 两个data question（Civáň japonica标签带不带括号；Panel B用raw 720还是
   filtered 718、2个UNK的技术状态）**没有替用户决定，等确认**。
4. **capture bait BED全仓库4个分支都搜不到**，capture track（CAPTURE.TV/
   CAPTURE.ALL）完全被卡住，shotgun track不受影响。

**v2脚本清单**（全部`scripts/ecotype_pca_v2/`下）：
`config/ecotype_pca_v2.yaml`（唯一参数来源，冻结）、`lib_ecotype_v2.py`
（共享库）、`00_validate_inputs.py`（输入/工具/标签校验，SHOTGUN_READY与
BATCH_1_FULL两级门禁）、`01_make_panel_manifest.py`、
`02_convert_eigenstrat_for_plink.sh`（convertf转换）、`03_audit_panel.py`
（MAF/missingness/spacing审计）、`04_audit_720_ld.py`（720面板LD衰减，
分块+halo算法）、`05_intersect_panel_baits.py`（capture bait交集，**当前
被阻塞**）、`06_build_reference_sample_set.py`（每个panel的axis-builder
keep-list，含595/718-720/五标签的硬校验）、`07_make_fixed_markers.sh`
（geno→MAF→LD剪枝，冻结`*.fixed.snplist`）、`08_make_5kb_thinned_markers.py`
（Panel B的paperlike_5kb路线）。`tests/`下4个测试文件。

**v1本次会话修的两处真实bug**（`scripts/ecotype_pca/pseudo_haploid_call.py`，
commit `a4fb1e6`+`c058189`）：
1. TV/ALL两次运行共享同一条全局随机数流，导致共同的transversion位点在
   两次运行里抽到不同read——用合成BAM实测复现（同一位点同一seed，TV调用
   出0，ALL调用出2），改成按`(seed, contig, position)`哈希算每个位点独立
   随机数后，两次运行在共享位点上完全一致（已实测验证）。这个bug如果不修，
   **今后任何"TV和ALL结果不一致，是否因为transition信号"的判断都不可信**。
   v1目前唯一跑过的Civáň LOO测试只用了TV单轨，不受这个bug影响，不需要
   重跑；但今后一旦要正式对比TV/ALL两条轨迹，必须用修复后的版本。
2. `call_rate`拆成`eligible_site_call_rate`（原公式，覆盖深度主导）和
   `allele_match_rate_among_covered`（只看真正抽到read的位点，是数据质量
   信号），避免一个数字混淆两个问题。

**两个仍未解决、需要你确认的数据问题**（本次会话已报告，未替你决定）：
1. `pseudo_haploid_call.py`的`ignore_overlaps=False`要不要改——取决于
   besthit过滤后的BAM里是不是已经是合并(collapsed)的古DNA read，还是原始
   未合并的paired mate。需要在服务器跑：
   ```bash
   samtools view -c sample.bam
   samtools view -c -f 1 sample.bam
   samtools view -c -f 2 sample.bam
   ```
2. v2的两个data question（见上）。

---

## 📍 给接手人的启动指令（新会话/新窗口先读这里，2026-08-14更新）

**⚠️2026-08-14更新：上一版(2026-08-13深夜)说"批量铺开已开始"是不准确的
——用户已明确说明，除了下面第①项那一次Civáň smoke test，其余全部PCA
（包括poplist bug修复后的重跑、16样本×3panel批量铺开、任何MAF/LD
pruning）**实际上都还没有在服务器上真正跑过**，之前给的都只是命令，
不要假设它们已经执行，接手后第一件事是用`squeue -u $USER`/`ls`核实
服务器上到底有什么，不要相信文档里任何"已提交"的措辞。真实状态：

1. **群体标签（3.2节待办1/2/2c）+ UNK剔除：三个panel全部完成，且已
   验证过**。`build_29m3k_population_labels.py`/`build_720_population_labels.py`/
   `build_civan_population_labels.py`，匹配率 99.2%/99.7%/99.9%；
   `filter_panel_by_label.py`产出`.filtered.*`文件——**新脚本一律用
   这批文件，不要用没过滤过的原始panel文件**。服务器路径见
   `docs/ECOTYPE_PCA_PHASE1_COMMANDS.md`。
2. **唯一真正在服务器上跑成功过的PCA**：Civáň panel、LV7008416379、
   leave-one-out正对照——`docs/ECOTYPE_PCA_PHASE1_COMMANDS.md`第6节，
   结果真实（aromatic最近，dist=0.0083；LOO模拟坐标跟已知样本真实
   投影几乎完全重合，PC1/PC2差距在小数点后4位，REF/ALT编码方向由此
   直接证实无需`--swap-ref-alt`）。**但这次跑的时候，axis-building
   的poplistname里混进了456个O._rufipogon野生稻样本**——见下一条，
   这个结果的"aromatic最近"结论建立在一个已知bug之上，还没有用修复
   后的版本重新验证过。
3. **已修复且已重新验证的bug**（commit `cd5de6a`修复，2026-08-15验证）：
   `run_sample_panel_pca.sh`原来把merged.ind里除Ancient外的所有标签
   都塞进poplistname，导致Civáň panel的野生稻样本参与了axis-building
   ——这跟Civáň论文自己的方法相反(论文只用595份栽培稻定义axes，野生
   稻只做projection)。修复：新增可选第9个参数`REFERENCE_LABELS_FILE`，
   配套`scripts/ecotype_pca/civan_domesticated_reference_labels.txt`
   列出6个栽培稻标签。**验证结果**：用修复后的脚本重跑
   LV7008416379×Civáň（输出在
   `gene/results/ecotype_pca/civan_refonly_check/`）。**注意坐标不能
   跨版本直接比较数值**——旧结果的axes由栽培+野生共同定义，新结果的
   axes只由595份栽培稻定义，两次是不同的PCA坐标系(各自独立的旋转/
   尺度)，PC1/PC2从(0.0283,-0.0074)变到(0.0393,-0.0103)、distance从
   0.0083变到0.0086，这些数值层面的接近是巧合，不代表"结果几乎没变"。
   真正可比、也确实成立的是**排序**：两次都是aromatic第一，其次
   japonica系列，indica/aus/野生稻明显更远。结论：野生稻混入建轴**
   不是**"aromatic最近"这个结果的成因(已排除)，但这还不等于验证了
   "这个样本就是aromatic"——147个位点本身有没有区分aromatic/japonica
   的能力还没测过，见`docs/ECOTYPE_PCA_PANEL_QC_DESIGN.md`第0节新增的
   confusion-matrix待办(用`simulate_leaveoneout_projection.py
   --mask-from`批量遮蔽已知标签样本、比较真实标签vs投影最近标签，
   现有脚本循环用法即可，不需要新脚本)。**仍然只是1个样本、
   147个位点、单次伪单倍型抽样的冒烟结果**，bootstrap不确定性量化
   (执行计划第6节)还没做，也不代表可以现在就把16样本×3panel的48
   组合批量跑开——见QC设计文档第6节第4项，reference-first架构重构
   还没落地，48组合批量目前跑了以后大概率要扔掉重跑。
4. **2026-08-13当天问了GPT关于三个panel MAF/LD pruning的问题**，得到
   一整套架构性反馈（reference-first、冻结marker set、单次smartpca
   覆盖所有古样本而不是每个样本单独跑一次子集），**全部记录在新文档
   `docs/ECOTYPE_PCA_PANEL_QC_DESIGN.md`，这是仅次于本节的第二重要
   入口文档**——里面有完整的待办顺序(第6节)、每个panel的具体参数
   建议、哪些是GPT转述还没独立核实过的，都写清楚了，不要在新会话里
   重新问一遍GPT或重新设计一遍。**720 panel被标记为最紧急**：这份
   6.7M SNP矩阵是论文第一作者私下发的加密版本，来源和处理流程完全
   不透明，跟论文本身发表的60,722-marker分析集不是一回事。
5. **三个panel目前都没有做过MAF/LD pruning**，只做过UNK剔除——这是
   `ECOTYPE_PCA_PANEL_QC_DESIGN.md`记录的核心待办，还没开始。**2026-
   08-15更新：这条已经在`scripts/ecotype_pca_v2/`里用冻结统计设计
   重新实现（见本文档最上方新增小节），但v2还没在服务器真实数据上
   跑过，v1这边仍然维持"没做过MAF/LD pruning"的原状，不要混淆两套
   代码的进度。**
6. **样本专属子集PCA(`build_sample_panel_subset.py`)的设计缺陷**：
   每个古样本单独跑一次、每次都把panel缩到该样本自己覆盖的SNP——这
   意味着不同古样本的"PC1"其实是不同marker set上算出来的，不能直接
   拿来比较（比如画时间轨迹）。`ECOTYPE_PCA_PANEL_QC_DESIGN.md`第1节
   给出了正确设计(reference-first、每个panel冻结一套marker set、单
   次smartpca覆盖全部古样本)，**这是当前最大的一块未完成重构，还没
   开始动手**。**2026-08-15更新：这条重构就是`scripts/ecotype_pca_v2/`
   要解决的问题，设计已经落地成脚本，但还没在服务器上跑，见本文档
   最上方新增小节。**

**集群基础设施问题（跟本分支设计无关，但会直接导致SLURM作业失败，
已写入长期记忆）**：`/itp`挂载点不覆盖`node06`(曾经也包括临时down掉
的`node05`)，任何作业路径经过itp软链接又调度到node06就会瞬间失败——
`sbatch`/`srun`一律加`--exclude=node05,node06`。

**读文档顺序**：先读本节 → `docs/ECOTYPE_PCA_PANEL_QC_DESIGN.md`(GPT
review + MAF/LD设计 + 待办顺序) → `docs/ECOTYPE_PCA_PHASE1_COMMANDS.md`
(标签匹配→矩阵瘦身→样本专属子集PCA→leave-one-out的服务器命令原样
记录，跟`ECOTYPE_PCA_PHASE0_COMMANDS.md`一个风格)。本文档回答"为
什么"，PHASE1_COMMANDS"照抄命令用"，QC_DESIGN"下一步做什么、参数
怎么定"。**v2的对应入口是`docs/ECOTYPE_PCA_V2_SPEC.md`（冻结规范）
+ 本文档最上方新增小节（现状）。**

---

## 0. 现状一句话总结

**⚠️2026-08-12更新：pseudo-haplotype调用+smartpca投影这条链路的四个脚本
已经全部写完并推送，但一个都还没在服务器真实数据上跑过**——①BWA映射
(`map_besthit_to_irgsp.sh`)、②pseudo-haplotype调用
(`pseudo_haploid_call.py`)、③合并进panel(`merge_ancient_into_panel.py`)、
④smartpca参数模板(`par.PROJECT.template`，**故意留空**，等群体标签
待办1/2/2c完成才能真正填)。①的运行命令已经给了用户，**结果未知**。
完整清单、路径、待运行/待核对命令，见📍节+5.8节，**这是当前接手最应该
先看的两处**。

**⚠️2026-08-11晚更新：PCA-C(Civáň)的convertf转换已经在服务器上跑完，
成功**——`civan_snp.eigenstratgeno/.snp/.ind`已经产出，2,365,188个SNP、
1056个样本，跟论文声称的数字精确吻合(这也顺带证实了VCF文件本身没有被
截断，见5.6节)。**PCA-A(29M_3k)、PCA-C(Civáň)现在都有现成的EIGENSTRAT
输入了，PCA-B(6.7M_720)本来就是EIGENSTRAT**——三个panel的格式转换全部
完成。

**⚠️2026-08-11更新：三级祖源框架定案，Civáň 2019面板作为PCA-C桥接层
插入到PCA-A/PCA-B之前，见第5节**（这是当前最新的设计层，第2节两个独立
PCA的设计仍然成立，只是不再是"最先一步"，现在的顺序是：先用Civáň统一面板
判断栽培/野生/混合，再分别用3K和720面板细分亚群/谱系）。

**⚠️2026-08-08深夜重大方向调整：放弃合并两个panel，改为两个独立PCA**。
`29M_3k`与`6.7M_720`求交集这条路线走到底了——先是mergeit按SNP ID字符串
匹配、两边ID命名方式不同(`1026` vs `1np1409`)导致0匹配（见下方历史记录），
用ID统一重命名的方式修复后重新跑，交集数量**依然太少**（chr1单条染色体
两边分别305万/78万个位点，真实重叠才132个，全基因组交集大概率只有一两千
的量级）——这个密度对"每个古代样本还要再对着这份交集去找自己reads覆盖到
的位点"这一步来说完全不够用，古代样本reads本来就稀疏，一两千个位点里能
抽到的更是所剩无几。**决定：不再合并，`29M_3k`和`6.7M_720`各自独立做一次
smartpca**，古代样本分别投影到"栽培稻多样性"和"野生稻多样性"两个坐标空间
里，分别回答"更接近哪个栽培亚群"和"更接近哪个野生谱系"，不强求画在同一张
图上。详见第2节。

好处：convertf那步(29M_3k PLINK→EIGENSTRAT，4.4小时CPU)完全没有浪费——
转换好的`NB_final_snp.eigenstratgeno/.snp/.ind`正好就是独立PCA-A需要的
完整输入；`6.7M_720`本来就是EIGENSTRAT格式，也是现成的。不用再纠结
mergeit的ID匹配、strand-ambiguous剔除这些坑。标签来源也顺带简化：
`29M_3k`这边不用再等`OrA-OrF`，3K RGP官方品种元数据
(`docs/references/3k_rice_genomes_project/rice_line_metadata_20141029.xlsx`，
main分支`docs/LITERATURE.md`第3节)自带IND/AUS/ARO/TRJ/TEJ/ADM标准亚群
标签，直接能用；`OrA-OrF`这条线只需要给`6.7M_720`野生稻panel那边用。

**mergeit求交集的调试历史（存档，不再是行动项，仅供以后类似问题参考）**：
坐标系核对完成(见3.1)：染色体命名一致、REF/ALT方向虽相反但mergeit本可
自动处理。convertf这步(29M_3k PLINK→EIGENSTRAT)在服务器上验证跑通
(2026-08-07)：29,635,224个SNP、3024个样本转换前后完全一致，无丢失。
`asn720data`一度被误判"密度太低、直接弃用"，2026-08-07推翻——它的`.fam`
FID列(`OrA-OrF`)是野生稻这边真实的群体标签来源。2026-08-08晚间跑
mergeit，`##end of mergeit`正常结束但`Histogram of checkmatch return codes`
显示`total: 0`——诊断出两个panel的SNP ID命名方式不同(29M_3k是纯数字物理
位置如`1026`，720是`{chrom}np{pos}`如`1np1409`)，mergeit按ID字符串匹配、
两者永远对不上；这个诊断在另一个并行会话窗口(疑似Codex CLI，工作副本在
本机`/Users/inmpain/Documents/angkor/rice_adna_pipeline_publish`)和本会话
里被独立验证过两次，用`comm`工具在chr1上实测确认真实物理位置重叠132个
(而非0)，证实ID格式不同确实是0匹配的原因。用ID统一重命名(`{chrom}_{pos}`)
的方式修复后重新跑，交集数量出来了但太少，促成了这次放弃合并的决定。
`scripts/ecotype_pca/par.MERGE`和`run_convert_merge.sh`仍保留在仓库里
作为这段调试过程的记录，**不再是当前分析路径的一部分**，不要在新设计里
继续用它们。

**⚠️2026-08-11更新：besthit的硬阻塞可能已经解除，但有一个新的前置问题
待评估**——besthit分支`docs/ORYZA_BESTHIT_HANDOFF.md`第0/6.3节显示v2脚本
已经在服务器真实数据上对**全部16个古代样本**跑完正式批(`submit all`+
`merge`)，16/16一致性自检通过，整体留存率13.72%，早就不是"只有4个样本"
的旧状态了(那是2026-08-08晚的快照，本文档之前一直没同步，接手时不要信
这句话之前的版本)。**但**besthit v2的KEEP范围是**整个Oryza属**(不只
sativa/rufipogon/nivara三个目标种)，而本文档PCA-A/PCA-B/PCA-C三个panel
都只覆盖目标AA基因组复合群——这意味着直接拿`<sample>.besthit_oryza.
fastq.gz`喂给PCA，可能混入其他15个Oryza种的reads，产生虚假投影信号。
已经在该分支记录了一条"按taxonomic tier分级输出KEEP集合"的建议(见该分支
第7.6节)，**这条建议是否已经实施、`target_aa_complex.fastq.gz`这类分级
输出是否已经存在，接手第一件事需要去besthit分支核实**，核实之前①A/①B/
①C不要直接用全量`besthit_oryza.fastq.gz`当作"确认是目标AA复合群"的输入。
**（2026-08-12现状：这一条尚未核实，目前①BWA映射脚本直接用的是全量
`besthit_oryza.fastq.gz`，如果besthit那边后续真的加了taxonomic tier
分级输出，①的输入路径需要改成`target_aa_complex.fastq.gz`重新跑）**

## 1. 数据库

### 1.1 29M_3k(驯化稻，已下载)

- 路径：`/home/scratch/yinmt202607/db/29M_3k/`
  ```
  NB_final_snp.bed.gz
  NB_final_snp.bim.gz
  NB_final_snp.fam.gz
  ```
- 内容：对应IRRI页面的`3K RG 29mio biallelic SNPs Dataset`，3024份3K RG材料，
  约2900万个biallelic SNP位点，PLINK bed/bim/fam格式(gzip压缩)，
  **called vs Nipponbare MSU7/IRGSP1.0 genome**——与本项目主参考
  `db/asian_rice_panel_index/irgsp.fa` 坐标系一致，理论上不需要转换。
- **为什么选29mio而不是页面上其他几个3K数据集**：
  - 不选 `32mio "Full 3K RG SNPs"(含multi-allelic)`——只有tabular格式，没有
    PLINK bed/bim/fam，且含多等位位点，跟 pseudo-haplotype/smartpca 这套
    流程(默认biallelic)不兼容，需要额外写解析器，收益(比29mio多约9%位点)
    不值得这个成本。
  - 不选 `18mio Base SNP` / `4.8mio filtered SNP` / `1M GWAS SNP` /
    `404k CoreSNP` / `160k HDRA共享位点`——密度依次更低。本分析的 SNP 密度
    诉求不是"现代群体PCA需要多少独立位点"(这个404k甚至160k就够)，而是
    "古代样本稀疏reads随机覆盖到panel位点的概率"——panel 越密，这个概率
    越高，所以反而应该选最密的可用集合，29mio 是密度和格式兼容性的最优点。
- **现在直接可用**：`scripts/ecotype_pca/`下convertf已经把这份数据转成
  EIGENSTRAT格式(`NB_final_snp.eigenstratgeno/.snp/.ind`)，可以直接作为
  独立PCA-A的输入，不需要再等交集。

### 1.2 6.7M_720 panel(野生稻为主，已上传)

- 路径：`/home/scratch/yinmt202607/db/6.7M_720/`
  ```
  asn720.6m.geno
  asn720.6m.geno.gz    (压缩版，跟未压缩版并存)
  asn720.6m.ind
  asn720.6m.snp
  ```
- **格式是EIGENSTRAT**(`.geno/.ind/.snp`三件套)，本来就是smartpca的原生
  输入格式，可以直接作为独立PCA-B的输入。
- 内容：720份样本，约670万个SNP位点，以野生稻为主
- 来源/参考坐标系：**尚未最终确认是否也是 called vs Nipponbare
  MSU7/IRGSP1.0**——但check_ref.py抽查200个位点，91.5%(183/200)的REF/ALT
  与IRGSP参考序列在对应位置的碱基能对上，间接支持坐标系基本一致（不是
  完全不同的组装版本），只是不如29M_3k那么干净。既然现在不再需要跟
  29M_3k求交集，这个问题的紧迫性降低——只影响PCA-B自身smartpca跑出来的
  可信度，不再是"求交集会不会静默出错"这个问题了。
- **⚠️`asn720data`是这条线群体标签的关键来源**：
  `asn720data/asn720.pop.{bed,bim,fam}`(720份、94,974个SNP位点)本身的
  **基因型数据**密度太低，不适合拿来做PCA底层数据；但`asn720.pop.fam`的
  **FID列(第一列)记录着`OrA`/`OrB`/`OrC`/`OrD`/`OrE`/`OrF`这样的群体标签**，
  例如：
  ```
  OrD    ERR068594    0    0    0    1
  OrD    ERR068597    0    0    0    1
  OrA    ERR068598    0    0    0    1
  OrF    ERR068600    0    0    0    1
  OrC    ERR068604    0    0    0    1
  ```
  这批`ERR0685xx`风格的样本ID，跟`asn720.6m.ind`里那批被记成笼统`control`
  标签的样本**是同一套ID体系**(见下面"样本ID观察")。**群体标签应该从
  `asn720data/asn720.pop.fam`按样本ID(IID)匹配过去**，而不是用
  `asn720.6m.ind`自己那个占位符式的`control`。
- **样本ID观察**：`asn720.6m.ind`里样本ID是两种风格混合——`ERR068594`一类
  (ENA测序run编号，群体标签目前记成笼统的`control`，但真实标签应该从
  `asn720data/asn720.pop.fam`按ID匹配，见上一条)和`B011_merged`一类(跟
  3K RG自己的品种编号风格接近，见besthit分支`ORYZA_BESTHIT_HANDOFF.md`
  第6节SV数据里出现过的`B071`类样本名)。说明这720份很可能是从不同批次/
  不同命名体系拼合而成，不是单一来源；`asn720data/asn720.pop.fam`目前只
  见到`ERR`风格的条目，`B0xx_merged`风格的样本有没有对应标签、标签在哪
  还未核实。
- **新发现，内容未读**：`db/wild_rice_pangenome_README.txt`——在`db/`
  顶层目录，跟`db/16/`、`db/3k/`、`db/29M_3k/`、`db/6.7M_720/`、
  `db/asn720data/`平级，此前从未被提及过，很可能是解释`OrA-OrF`具体定义、
  或`db/3k/wild/`140+野生稻组装身份(见besthit分支`ORYZA_BESTHIT_HANDOFF.md`
  第7.2节的老问题)的权威说明文件。已请用户`cat`出内容。

### 1.3 Civáň 2019 统一栽培-野生面板（PCA-C，2026-08-11 convertf转换
已完成，见5.6/5.7节完整记录）

- 论文引用与数据集DOI：main分支`docs/LITERATURE.md`第2.2节（已用WebFetch/
  WebSearch核实，非转述）
- **确切服务器路径**：`/home/scratch/yinmt202607/db/paper1/`
  ```
  sativa-rufipogon_SNPs.vcf.gz    # 核SNP矩阵，VCF格式，bgzip压缩
  1825_Oryza_cpDNA.fastq          # 叶绿体基因组集，⚠️真FASTQ格式(见下)
  1825_Oryza_cpDNA.fastq.tar      # 上面那份的tar打包版，内容应该相同
  Table_S1.csv                    # 样本元数据(物种/群体/来源等)
  Table_S2.csv                    # 叶绿体组装质控表(测序深度/覆盖度等)
  Table_S3.csv                    # qSH1/qSD1-2/SPS1三个驯化位点单倍型频率
  civan_snp.bed/.bim/.fam         # plink2从VCF转出的中间PLINK格式
  civan_snp.eigenstratgeno/.snp/.ind  # 【PCA-C直接输入】convertf最终产出
  ```
- 核SNP矩阵：**1,056份样本(595栽培[283 indica/154 japonica/124 aus/34
  aromatic] + 461野生)，2,365,188个双等位SNP**，坐标系IRGSP-1.0，格式VCF。
  **坐标系已核实一致**：`check_ref.py`的VCF模式对200个位点抽查，REF列与
  `irgsp.fa`对应位置碱基**200/200完全匹配**——VCF的REF列本来就该是参考
  基因组本身，跟29M_3k/720那种A1/A2方向模糊的情况不同，这个100%匹配是
  "坐标系统真的对得上"的强证据。**VCF→EIGENSTRAT转换已在服务器上跑完
  （2026-08-11晚）**：convertf不支持VCF输入(README确认，见5.6节)，走的是
  `VCF→(plink2 --vcf --make-bed)→PLINK bed/bim/fam→(convertf)→EIGENSTRAT`
  两步路径，脚本是`scripts/ecotype_pca/convert_civan_vcf.sh`+
  `par.CIVAN.PLINK.EIGENSTRAT`。**最终产出2,365,188个SNP、1056个样本，
  跟论文声称的数字完全精确吻合**——这也**间接证实了VCF文件本身没有被
  截断**（之前"总行数因Ctrl-C中断没确认"这个悬而未决的问题，现在可以
  认为已经解决：如果文件被截断，转换出的SNP数不可能刚好精确等于论文
  声称的数字）。转换过程中出现过一条`all individuals set ignore...
  resetting all individuals...`的警告，最终`indivs: 1056`(全量，没有
  样本被误删)，见5.6节对这条警告的解释——是PLINK`.fam`第6列phenotype
  默认值(-9)触发的convertf已知解析行为，convertf自己检测到异常并恢复，
  不是数据丢失，下次转换类似来源的PLINK文件遇到同样警告不用紧张，但
  最终产出的`indivs`数字仍然要核对一遍。**⚠️REF/ALT方向没有单独核对过**
  ——1.3节这里的"200/200完全匹配"核对的是VCF阶段，convertf转成
  `civan_snp.snp`之后列顺序有没有变、跟VCF阶段是否还一致，从未专门验证，
  见📍节第4条和5.8节。
- 叶绿体基因组集：1,825个，独立于核基因组的母系谱系证据层。**⚠️格式是
  真FASTQ，不是FASTA**（`grep -c "^>"`返回0是符合预期的，因为FASTQ用`@`
  不用`>`；`head`看到`@ERR605276 chloroplast, complete genome`这样的
  header、序列行、`+`分隔行、质量值行的标准四行一条记录结构）——质量值
  行看起来几乎全是`~`（Phred+33下`~`=Q93，接近满分），这种"整条质量值
  近乎恒定拉满"的模式是**拼接好的基因组序列被工具批量转成FASTQ格式时的
  典型占位符**（真实测序质量值不会这样几乎不变），不是这批叶绿体基因组
  本身测序质量异常好——**下游用这份数据时应该只取序列、丢弃/忽略质量值**
  （比如`seqtk seq -A`转回FASTA，或者直接用biopython的SeqIO按fastq解析
  拿序列部分）。**记录数已核实**：`wc -l`返回7300行，按FASTQ每条记录
  4行计算，7300/4=1825，跟论文声称的1825个叶绿体基因组**精确吻合**。
- **Table_S1/S2/S3.csv 完整行数（2026-08-11确认，含表头）**：
  ```
  Table_S1.csv   1064行 → 1063条数据(表头本身有大量尾随空列，Excel导出
                 常见现象，不是数据损坏)。跟论文摘要里"1,056份样本"这个
                 数字对不上(多7行)，原因未知——可能是额外的表头/脚注行，
                 也可能真的多出7个条目，不阻塞当前工作但记一笔待查
  Table_S2.csv   1646行 → 1645条数据。这张表是"叶绿体组装质控"表(测序
                 深度/覆盖度等)，1645跟叶绿体基因组总数1825对不上——一个
                 合理猜测是Table_S2只列了"本研究新测序组装"的部分，1825
                 的总数里还包含了论文引用的、之前已发表的叶绿体基因组，
                 这个猜测未经证实，不阻塞工作
  Table_S3.csv   86行 → 85条数据，是qSH1/qSD1-2/SPS1单倍型频率表(按
                 accession/群体分行)，量级看起来合理，不是这次待查重点
  ```
  这两处行数对不上的地方都不阻塞PCA-C convertf脚本的开发（convertf只
  需要VCF），留作以后需要精确核对样本身份/群体归属时再回头查。

### 1.4 两个panel的关系（历史设计，已废弃，存档）

~~分析设计里两者是互补而非替代关系：3K(29mio) 覆盖驯化稻的谱系多样性，
6.7M_720 补充野生稻/近缘种一侧——古稻样本理论上落在"驯化-野生"谱系的某个
位置，两个panel合起来才能给出有意义的PCA背景。~~ **2026-08-08已放弃**：
两个panel密度差23倍(chr1上305万 vs 78万)，真实交集只有一两千位点量级，
撑不起"合并成一个统一底盘"这个设计，改为两个独立PCA，见第2节。
**2026-08-11更新**：这个"统一坐标系"的诉求本身没有错，只是不该由
29M_3k+6.7M_720自己硬合并来满足——Civáň 2019面板天生就是栽培+野生同一
坐标系(见1.3/第5节)，是更合适的桥接工具。

## 2. 分析设计：两个独立PCA（2026-08-08版，仍然成立，现在是第5节三级
框架里的第2/3层）

```
besthit 过滤后的 Oryza reads (来自 codex/oryza-competitive-mapping 分支产出:
  <sample>.besthit_oryza.fastq.gz / 或对应BAM)
        │
        ├───────────────────────────────┬───────────────────────────────┐
        ▼                                ▼
   PCA-A：29M_3k(驯化稻)独立PCA      PCA-B：6.7M_720(野生稻为主)独立PCA
        │                                │
   ①A 该样本reads ∩ 29M_3k全部位点    ①B 该样本reads ∩ 6.7M_720全部位点
      (不再先跟另一个panel求交集，     (同左)
       直接用各自完整密度，古代样本
       覆盖到位点的概率更高)
        │                                │
   ②A pseudo-haplotype调用             ②B pseudo-haplotype调用
      (每个位点随机抽一条覆盖read       (同左)
       取碱基，标准古DNA pseudo-
       haploid做法)
        │                                │
   ③A smartpca -lsqproject投影         ③B smartpca -lsqproject投影
      现代样本(3K全部3024份)建PCA       现代样本(720份)建PCA参考空间，
      参考空间，古代样本投影上去        古代样本投影上去
      标签来源：3K RGP官方亚群标签      标签来源：asn720data的`OrA-OrF`
      (IND/AUS/ARO/TRJ/TEJ/ADM，       (第3.2节第1条，待核实覆盖率)
       docs/references/3k_rice_
       genomes_project/)
        │                                │
        ▼                                ▼
   古代样本在栽培稻多样性里的位置    古代样本在野生稻多样性里的位置
        │                                │
        └───────────两个独立坐标空间，不能直接叠加对比───────────┘
                     分别解读："更接近哪个栽培亚群" /
                     "更接近哪个野生谱系"，合起来综合判断生态型
```

**与旧设计(单一合并panel)的区别**：旧设计想要一张"驯化-野生"统一谱系图，
这个方案牺牲了这一点，换来两条独立、密度更高(古代样本覆盖概率更高)的
路径。这个取舍是2026-08-08晚上根据mergeit实测交集量级（chr1上只有132个
真实重叠位点）做出的决定。**2026-08-11更新**：这个牺牲现在被Civáň面板
(第5节)部分补回来了——PCA-A/PCA-B仍然各自独立、密度优先，但前面多了一层
Civáň统一面板先给出"栽培/野生/混合"这个粗判断，两边不再是完全脱节的两张
图，而是同一棵决策树的两个下游分支。

**①A/①B"先subset到该样本实际覆盖的位点"这一步的意义，不只是逻辑上必须
(古代样本本来就不可能在panel全部位点上都有read覆盖)，还直接决定②/③两步
的计算量**——不管是29M_3k还是6.7M_720，②pseudo-haplotype调用和③smartpca
投影都只需要在"该样本reads实际覆盖到的那一小撮位点"上跑，不需要在panel
全部2900万/670万个位点上跑。16个古代样本、每个样本都远比整个panel稀疏，
**每个样本应该各自生成一份自己专属的、小得多的位点子集**，而不是对所有
样本都套用同一份大面板重复计算——这是2026-08-08用户明确要求写清楚的设计
原则。**这一点2026-08-11被GPT质疑过，见第5.2节的详细技术核对**——结论是
这条设计原则本身没有问题，GPT的担忧建立在对`-lsqproject`工作方式的误解上，
但GPT指出的"缺失数据导致投影收缩"是真实存在的独立问题，已经采纳进第5.2节
的补充设计。**2026-08-12：这个设计已经落地成`pseudo_haploid_call.py`
（见5.8节），真实机制细节在那里。**

## 3. 待确认/待办(按优先级)

### 3.1 坐标系核对结果(已用 `scripts/ecotype_pca/check_ref.py` 完成)

**染色体命名**：`29M_3k`原始bim(裸数字1-12) 与 `6.7M_720`(`asn720.6m.snp`，
同样是裸数字1-12) **命名方式一致**——之前把`NB_final_snp.bim`手动改成过
`chr01`风格(`NB_final_snp.bim.chrfix`一类)，**这一步是多余的、不需要**，
后续统一使用原始的 `NB_final_snp.bim.orig`(裸数字版)即可，不要用改过
染色体名的那份。

**REF/ALT等位基因方向**：
- `29M_3k`(`NB_final_snp.bim`，200/200抽查)：**A2=REF, A1=ALT**，非常干净，
  完全符合PLINK从VCF导入的标准默认约定。
- `6.7M_720`(`asn720.6m.snp`，200个位点抽查)：**A1=REF, A2=ALT 占183/200
  (91.5%)**，跟29M_3k方向相反。⚠️ 剩下17/200(8.5%)不符合这个反向规律，
  不是单纯的"顺序颠倒"能解释的，可能是720面板内部个别位点的数据质量问题，
  暂不处理。**现在两个panel各自独立跑PCA，这个反向问题不再需要处理**——
  smartpca只看自己panel内部的REF/ALT一致性，不涉及跨panel比较。**但
  `pseudo_haploid_call.py`需要知道这个方向才能正确输出0/2，这个方向的
  权威判断依据是执行计划文档第6节的leave-one-out模拟，不是这里的FASTA
  匹配率本身——匹配率低只说明原始数据标注习惯，不直接等同于0/2编码方向
  错误，见`docs/ECOTYPE_PCA_EXECUTION_PLAN.md`2.1节。**
- **Civáň VCF核SNP矩阵(200/200抽查)**：**REF列与irgsp.fa 200/200完全
  匹配**，见1.3/5.6节。这是三个panel里坐标系验证最干净的一个。**但这是
  VCF阶段的核对，convertf转成`civan_snp.snp`之后有没有变、要专门重新
  核对，见5.8节。**

### 3.2 待办(按优先级，2026-08-12更新)

0. **【已完成】确认Civáň 2019数据在服务器上的确切路径，并核对下载
   完整性**——详见1.3/5.6节。
0c. **【已完成，2026-08-11晚】写并跑通Civáň VCF→EIGENSTRAT转换**——
    `scripts/ecotype_pca/convert_civan_vcf.sh`+`par.CIVAN.PLINK.EIGENSTRAT`，
    走`VCF→plink2→PLINK→convertf→EIGENSTRAT`两步路径(convertf本身不
    支持VCF，见5.6节)，`civan_snp.eigenstratgeno/.snp/.ind`已产出，
    2,365,188 SNP/1056样本精确匹配论文数字，详见1.3/5.7节。
0e. **【已完成，2026-08-12：四个脚本写完但都还没跑过】** pseudo-haplotype
    调用+smartpca投影四步流水线全部脚本已写完并推送，详见5.8节完整
    清单。**下一步不是继续写代码，是先跑通①(BWA映射)拿到结果**，见
    📍节和5.8节。
0d. **【已完成，2026-08-13】群体标签(待办1/2/2c)**——三个panel全部
    做完，④(smartpca)已经用真实`poplistname`跑通并且leave-one-out
    验证通过，见📍节。
0b. **【同等优先级，尚未核实】去besthit分支核实"按taxonomic tier分级
    输出KEEP集合"这条建议(第0节末尾/第4节末尾)是否已经实施**——如果
    还没有，需要决定是在besthit分支加这个分级逻辑，还是本分支自己对
    `besthit_oryza.fastq.gz`做一次后处理(用`decisions.tsv`里现成的
    `best_oryza_taxid`列筛出sativa/rufipogon/nivara子集)，两种做法都
    可行，但要先选一种，不要两边各写一份。**2026-08-12：①的BWA映射
    脚本目前直接用的全量`besthit_oryza.fastq.gz`作为输入，如果这条
    建议后续被besthit分支采纳，①要改成用`target_aa_complex.fastq.gz`
    重新跑一遍。**
0f. **【新增，2026-08-12同日稍晚，当前最高优先级的新开发工作】** 执行
    `docs/ECOTYPE_PCA_EXECUTION_PLAN.md`——先做Phase 0 IRGSP覆盖普查
    （不建PCA），再开发样本专属panel子集脚本和leave-one-out模拟脚本。
    该文档第1节有完整的10步顺序，第2节记录了本次已经落地的三处代码
    修复(REF/ALT判断逻辑纠正、比对去重对齐主流程、合并脚本硬检查)。
1. **【已完成，2026-08-13】`6.7M_720`独立PCA的群体标签**——
   `build_720_population_labels.py`，718/720 (99.7%) 匹配。真实情况
   比预想的复杂：`asn720.pop.fam`不是纯野生稻标签，混了栽培稻锚点
   (IND/AUS/ARO/TRJ/TEJ/ADM)；`_merged`后缀样本(310个)在`.pop.fam`里
   完全没有对应条目(只有ERR/SRR风格ID)，改用剥掉后缀后复用29M_3k同一份
   3K RGP元数据解决。`db/wild_rice_pangenome_README.txt`确认是**空
   文件**，`OrA-OrF`跟Nature 2025论文`Or-Ia/Or-Ib/Or-II/Or-IIIa/
   Or-IIIb`是否同一套体系**仍未确认**，标签原样保留`OrA-OrF`，不做
   等价假设。`OrADM`（推测是野生稻侧的admixed，未证实）、`RAY`（9个
   样本，含义完全未知）都原样保留，未强行解释。
2. **【已完成，2026-08-13】`29M_3k`独立PCA的群体标签**——
   `build_29m3k_population_labels.py`，3000/3024 (99.2%) 匹配。IRIS_313
   风格ID(2466个)按Table S1A精确匹配；B0xx/CX风格ID(总558个)都在
   Table S1B里(分属"MC"/"IRMBN"两个内部来源，之前以为CX完全没有元数据
   是错的，是没查全导致的误判)。8种原始Variety Group取值里6种干净映射
   到IND/AUS/ARO/TRJ/TEJ，另外`Intermediate type`(135)、`Japonica`
   未细分(132)**刻意不**强行归类，保留原样标签(`INTERMEDIATE_TYPE`/
   `JAPONICA_UNSPEC`)，不参与建轴但不删除。24个UNK已从矩阵物理删除
   (`filter_panel_by_label.py`)。
2c. **【已完成，2026-08-13】Civáň PCA-C的群体标签**——
    `build_civan_population_labels.py`，1055/1056 (99.9%) 匹配。样本ID
    是plink2自动生成的`ACCESSION_ACCESSION`自重复格式(部分样本本身带
    FID、没被重复)，`recover_accession()`按字符串中点位置切分而不是
    按下划线数量(因为accession本身可能带下划线，如`IRIS_313-9986`)。
    野生样本(461个)在panel里是ERR号，但`Table_S1.csv`按`W####`编号，
    靠`Table_S2.csv`(名义上是"叶绿体组装质控表"，实际带`Accession↔SRA
    dataset used`精确桥接)解决，460/461桥接成功。**发现一个真实
    数据问题并修正**：Civáň面板的野生样本不是清一色O. rufipogon，
    有4个是完全不同的种(O. meridionalis/O. glaberrima/O. barthii/
    O. longistaminata)，`Table_S1.csv`的Group列对所有野生样本都是
    占位符"-"，脚本改为Group为"-"时退回用Species列，避免这4个种被
    误并入O. rufipogon这一个标签。1个最终UNK(461-460桥接失败的那个)。
3. **【已完成并验证，2026-08-13】四步流水线在真实数据上跑通**——
   `build_sample_panel_subset.py`+`simulate_leaveoneout_projection.py`
   补齐了样本专属子集化和leave-one-out模拟这两块，在Civáň panel用
   LV7008416379完整跑通①②③④全链路，leave-one-out正对照通过(见📍节)。
   第5.2节的bootstrap点云版本仍然是下一层增量工作，不是当前优先级。
   **当前优先级：16样本×3panel批量铺开**（`run_sample_panel_pca.sh`，
   进行中，见`docs/ECOTYPE_PCA_PHASE1_COMMANDS.md`）。
4. **核查`asn720.6m.ind`与`NB_final_snp.ind`之间的样本ID重叠**(即
   `CX382`疑似重复案例)——去掉`_merged`后缀后系统性比对两份样本列表，
   确认重叠规模。优先级较低，不阻塞主线，但正式产出结果前必须处理。
5. **mergeit相关内容已废弃，不再是待办**——`scripts/ecotype_pca/par.MERGE`、
   `run_convert_merge.sh`保留在仓库里仅作调试历史记录，不要再执行。

## 4. 与 besthit 分支的关系

本分析线的输入直接依赖 `codex/oryza-competitive-mapping` 分支产出的
`<sample>.besthit_oryza.fastq.gz`(见该分支 `docs/ORYZA_BESTHIT_HANDOFF.md`
第5.3节)。besthit 那边的分析(第7/8节:acc2taxid 覆盖度诊断、野生稻组装
物种身份诊断)结果，间接也会影响本分析线——如果 besthit 阶段发现需要扩充
competitive mapping 的参考库(把 `db/3k/wild/` 那140+野生稻组装加进去)，
"确认是Oryza"的read集合本身会变化，本分析线①A/①B步要用的输入也要跟着
重跑。两条线目前互相独立推进，但下游会汇合，需要留意上游变动。

**2026-08-11新增，已转发并落地为besthit分支第7.6节开放问题**：GPT建议
besthit的KEEP集合按taxonomic tier分级输出(`oryza_genus.fastq.gz`全属 /
`target_aa_complex.fastq.gz`目标AA基因组复合群 / `other_oryza.fastq.gz`
其他Oryza种)，而不是现在v2版本的单一"全Oryza属都算KEEP"。理由：
Civáň/3K/720三个PCA panel主要覆盖sativa/rufipogon/nivara这个目标AA复合群，
如果besthit的KEEP集合混入了其他15个Oryza种的reads，这些reads投到PCA
里会被强行投影到某个不相关的现代栽培/野生群附近，产生虚假信号。**这条
建议是否采纳、什么时候做，由besthit分支接手人根据其自身进度决定**（见
该分支`docs/ORYZA_BESTHIT_HANDOFF.md`第7.6节），本文档3.2节待办0b是
本分支这边需要做的对应核实工作。

---

## 5. 2026-08-11：三级祖源框架 + Civáň 2019桥接面板（GPT策略讨论整理）

### 5.0 背景

用户把本文档、`docs/RESEARCH_ROADMAP.md`、besthit分支
`docs/ORYZA_BESTHIT_HANDOFF.md`三份文档的raw链接给GPT看，请GPT结合读、
给下一步研发策略建议。GPT给出了一套三级祖源框架的方案，同时下载了
Civáň et al. 2019(见1.3节、main分支`docs/LITERATURE.md`第2.2节)的核SNP
矩阵、叶绿体基因组集、论文三张补充表到服务器。以下是这次讨论整理进仓库
的内容，包含GPT的原始建议、以及Claude Code(我)核对后的技术修正——**不是
GPT说什么就直接照抄**，凡是我核实过跟标准做法不符的地方都明确标注。

### 5.1 三级祖源框架（采纳）

```
古DNA真实性QC + besthit(确认Oryza属)
        │
        ▼
PCA-C：Civáň 2019统一面板 —— 先判断栽培/野生/混合
        │
        ├─栽培侧──▶ PCA-A：3K面板 —— 细分IND/AUS/ARO/TEJ/TRJ
        │
        └─野生侧──▶ PCA-B：720面板 —— 细分OrA-OrF
        
（叶绿体：独立母系证据，与核基因组结果并行报告，不作替代）
```

**为什么这个框架比"两个完全独立、互不通气的PCA"更好，理由是站得住的**：
Civáň面板的2,365,188个SNP同时覆盖栽培(indica/japonica/aus/aromatic)和
野生稻，是**天然的统一坐标系**——不需要像29M_3k+6.7M_720那样费力去对齐
两个各自独立构建的panel(那次mergeit失败正是因为两边从来没设计成能对齐)。
用它先做一次粗分类(栽培/野生/混合)，比直接用两个互不通气的PCA各自给一个
说不清楚谁更可信的结论要扎实。这一条**采纳**，作为当前设计的最新层，
写进第2节顶部。

**需要做的技术工作**：Civáň面板是VCF格式，需要convertf转EIGENSTRAT——
**已完成，见1.3/5.6/5.7节**。

### 5.2 对"每个古样本各自建一套PCA"这条批评的技术核对（部分修正）

GPT原话：

> 当前文档提出：每个古样本按自己覆盖的位点生成一份小面板，再各自跑
> smartpca。这个做法会导致每个样本使用不同SNP、不同特征向量和不同PC
> 坐标系……不同古样本的PC1/PC2不能直接比较。

**核对结果：这个批评不完全准确，建立在对`-lsqproject`工作方式的误解上**。
重新读第2节①③两步的设计：③步明确写的是"现代样本(3K全部3024份/720全部
720份)建PCA参考空间，古代样本投影上去"——也就是说，**PCA的特征向量
(eigenvector)从始至终只用现代样本算一次、固定住**，不会因为古代样本
不同就重新计算。①步"每个样本各自生成一份小面板"，指的是**该样本在哪些
位点上有pseudo-haploid基因型**(genotype-calling范围)，不是"该样本用来
建PCA轴的SNP集合不同"——古代样本本来就只能在自己实际covered的位点上有
非缺失基因型，其余位点对`smartpca`来说就是缺失值，这正是`-lsqproject`
模式设计出来要处理的场景(EIGENSOFT官方就是给古DNA/低覆盖度样本准备的
投影模式)：不同古代样本各自不同的缺失模式，被投影到同一组**已经固定**
的现代PC轴上，PC1/PC2坐标之间是可比的，可以画在同一张图上。这是Reich
lab等古DNA群体遗传学论文的标准做法，不是本文档的设计问题。**2026-08-12
更新：这个机制现在已经落地成`par.PROJECT.template`里的
`poplistname`+`lsqproject: YES`组合，见5.8节，不再只是理论描述。**

**GPT批评里真正站得住的部分，已经采纳**：低覆盖度样本投影到固定PCA轴上
时，`-lsqproject`确实有已知的**收缩偏差**(shrinkage/attenuation
bias)——covered位点越少，投影出来的坐标越容易被拉向原点/现代均值，这是
真实存在的统计伪影，不是"能否比较"的问题，而是"比较时要不要校正"的问题。
**采纳的补充设计（⚠️2026-08-12：这部分还没有实现，见3.2节待办3）**：
1. 对每个现代参考样本按对应古代样本的实际覆盖位点数做**下采样(downsample)
   匹配对照**，量化"覆盖度这么低、就算是一个真实的现代样本，投影坐标会
   收缩到什么程度"，作为该古代样本投影结果的校正基准/置信区间参考
2. 对每个古代样本做**100-200次pseudo-haplotype重复抽样**(每次独立在
   每个覆盖位点随机抽一条覆盖read取碱基)，每次重复各自投影一次，画出**投影
   点云**而不是单一个点，点云的分布形态本身就是一种不确定性量化——点云
   集中在一个现代群体附近，可以报告"主要接近该群体"；点云横跨多个群体，
   应该报告"混合来源或数据不足"，而不是强行归到某一个群体

这两条是第3.2节待办3的**增量**部分，等①②③④四个脚本先跑通单次调用版本
之后再实现——`pseudo_haploid_call.py`目前的`--seed`参数已经为将来支持
多次重复抽样留了接口，但脚本本身还是单次确定性调用，没有循环生成N个
副本的逻辑。**2026-08-12更新：这套不确定性量化的组合方案(多次抽样+区块
bootstrap+现代样本降采样对照+现代reads走同一流程)已经在
`docs/ECOTYPE_PCA_EXECUTION_PLAN.md`第6节写成具体待办，与leave-one-out
模拟脚本共用同一套实现。**

### 5.3 沉积物DNA不能默认是单一二倍体个体（采纳）

GPT指出：样本来自沉积物，可能混合多株水稻/多个栽培群/栽培与野生/不同
时期沉积/现代污染，单次pseudo-haploid调用出来的一个PCA点，不应该直接
解释成"一粒古稻的基因型"。**这个提醒是对的，直接采纳**，与5.2节的重复
抽样点云设计是同一件事的两个角度——点云不只是量化"低覆盖度导致的统计
噪声"，也同时是在诚实地面对"沉积物样本本身可能不是单一个体"这个更根本
的不确定性来源，两者都用同一套重复抽样+点云可视化的方法处理，不需要
分开设计两套流程。

### 5.4 新的延伸方向：功能位点桥接生态型证据缺口（新方向，记录，未启动）

Civáň论文Table S3给出`qSH1`(落粒)/`qSD1-2`(休眠/株高)/`SPS1`(穗粒结构)
三个基因在IRGSP-1.0上的坐标，以及japonica单倍型 vs 印度本地(非japonica)
单倍型在aromatic群体里的频率。这是`docs/RESEARCH_ROADMAP.md`工作线4
("生态型与耕作方式，目前完全空白")第一批具体的、有明确坐标可查的候选
功能位点，值得做两件事(都还没做)：
1. 和main分支`db/gene/flower_gene.txt`的57个开花基因清单交叉核对坐标，
   看有没有重合——如果重合，说明这57基因清单本身可能已经间接覆盖了部分
   驯化/生态型信号，之前只是没有从这个角度去解读过
2. 参考main分支`docs/RESEARCH_ROADMAP.md`工作线5(8.3节，57基因SV判生态型
   的方法论)里已经验证过的教训——**覆盖度QC必须先做**，DTH8/Ghd8那次
   "明星基因古DNA照样零覆盖"的先例同样适用于`qSH1`/`qSD1-2`/`SPS1`，
   不要一上来就冲进去分析，先用现有BAM查这三个位点古代样本的覆盖度。
   **2026-08-11：besthit分支8.3节独立提出了同样的思路，两边结论一致
   （见该分支`docs/ORYZA_BESTHIT_HANDOFF.md`第8.3节末尾）**

**这条是否要开新的git分支**：不需要。这本质上是"生态型功能位点证据层"，
逻辑上属于`docs/RESEARCH_ROADMAP.md`工作线4，技术上又跟本分支的PCA
古代样本pileup/覆盖度分析工具高度重叠(同样是"看古代样本在哪些位点有没有
覆盖")，**建议先在本分支(`codex/ecotype-pca-panel`)下新开一个子目录
(比如`scripts/domestication_loci/`)原型验证覆盖度够不够用，如果验证后
发现这是一条独立、值得长期投入的分析线(比如真的能稳定call出SV有无)，
再考虑要不要拆成独立分支——不要在验证覆盖度是否可行之前就先拆分支，
参考main分支8.3节"先做半天工作量的覆盖度QC，再决定要不要投入"这个已经
验证过的性价比原则。**

**2026-08-12更新**：又新增了一条同类工作——旱稻/水稻(upland/irrigated)
生态型panel的整理(Lyu et al. 2014等)，同样面临"是否新开分支"的问题，
结论与这里完全一致(不需要)，完整推理见`docs/ECOTYPE_PCA_EXECUTION_PLAN.md`
第10节。

### 5.5 besthit读集合分级（已转发并落地，见第4节末尾/besthit分支第7.6节）

见第4节末尾，这条建议本质上是besthit分支的工作，**已经贴过去了**（该分支
`docs/ORYZA_BESTHIT_HANDOFF.md`第7.6节，2026-08-11），本分支这边对应的
核实待办是3.2节第0b条。

### 5.6 Civáň数据服务器路径与完整性核对结果（2026-08-11晚，已完成）

用户在服务器上跑完了第5.6节(旧版)列出的核对命令，结果如下：

**1. 确切路径**：`/home/scratch/yinmt202607/db/paper1/`

**2. 核SNP矩阵**：`sativa-rufipogon_SNPs.vcf.gz`，标准VCF格式，`#CHROM`行
之后是`1/2/.../12`裸数字染色体名(跟29M_3k/720两个panel命名方式一致，
不需要额外转换)，样本列包含`B0xx`(3K RGP风格)/`CX0xx`/`ERR0685xx`(ENA
测序run号)/`IRIS_313-xxxx`(IRIS ID)/`W0xxx`(野生稻)等多种ID风格混合，
跟720 panel"多来源拼合"的样本ID观察(1.2节)是同一种模式。

**3. 坐标系核对**：用升级后的`check_ref.py`(新增`vcf`模式，见commit
`ac8a8eb`，直接读VCF的REF列跟irgsp.fa比对，不需要像PLINK/EIGENSTRAT那样
猜A1/A2方向)，200个位点抽查，**REF列与irgsp.fa 200/200完全匹配，0个
mismatch**——这是三个panel里坐标系验证结果最干净的一次，不像720那样有
91.5%的模糊匹配空间。

**4. VCF总行数**：**间接确认，见5.7节**。用户直接跑`zcat|grep -vc "^#"`
两次都因文件太大被Ctrl-C中断，没拿到确切数字，但5.7节记录的convertf转换
最终产出2,365,188个SNP，跟论文声称的数字精确吻合，等同于间接确认了VCF
文件完整、没有被截断。

**5. 叶绿体基因组集格式**：`1825_Oryza_cpDNA.fastq`，**真FASTQ格式**——
`head`看到`@ERR605276 chloroplast, complete genome`风格的header、序列
行、`+ERR605276...`分隔行、质量值行的标准四行一条记录结构，`grep -c
"^>"`返回0(符合预期，FASTQ本来就不用`>`)。质量值行几乎全是`~`
(Phred+33下`~`=Q93)，这是拼接好的基因组序列被批量转成FASTQ时常见的
**占位符质量值**，不代表真实测序质量——下游使用时应该只取序列部分，
忽略/丢弃质量值(比如用`seqtk seq -A`转回FASTA)。`wc -l`返回7300行，
按FASTQ每条记录4行算，7300/4=1825，跟论文声称的1825个叶绿体基因组
**精确吻合**，记录数这块完整性确认无误。

**6. 三张表完整行数**：`Table_S1.csv`=1064行(1063条数据，跟论文摘要
"1,056份样本"差7行，原因未知，不阻塞)；`Table_S2.csv`=1646行(1645条
数据，叶绿体组装质控表，跟叶绿体总数1825对不上，猜测是只列了本研究
新测序的部分、总数里还含之前发表的组装，未证实)；`Table_S3.csv`=86行
(85条数据，qSH1/qSD1-2/SPS1单倍型频率表，量级合理)。这些差异都记在
1.3节，不阻塞convertf脚本开发。

**7. convertf不支持VCF输入**：`cat CONVERTF/README`确认convertf只支持
5种格式(ANCESTRYMAP/EIGENSTRAT/PED/PACKEDPED/PACKEDANCESTRYMAP)，没有
VCF，也**没有`inputformat:`这个参数**(这一点跟29M_3k那次的PACKEDPED
坑是同一类问题的再次确认——convertf的输入格式判断完全靠文件后缀，不靠
任何显式参数)。服务器上只有plink2(无plink1.9)，`which plink2`确认
`~/.local/mamba/snakemake/bin/plink2`，版本`v2.0.0-a.6.9LM`。

### 5.7 Civáň VCF→EIGENSTRAT转换：跑通结果（2026-08-11晚，已完成）

**转换路径**：`VCF → (plink2 --vcf --max-alleles 2 --vcf-half-call
missing --chr-set 12 no-xy no-mt --make-bed) → PLINK bed/bim/fam →
(convertf, symlink .fam→.pedind，.bim原生支持不用改名) → EIGENSTRAT`。
脚本：`scripts/ecotype_pca/convert_civan_vcf.sh` +
`par.CIVAN.PLINK.EIGENSTRAT`（commit `cef87d5`/`f59a966`）。

**实测结果**：
```
before compress: snps: 2365188 indivs: 1056
after compress:  snps: 2365188 indivs: 1056
numvalidind: 1056  maxmiss: 1056001
numsnps output: 2365188
##end of convertf: 442.305 seconds cpu  2100.630 Mbytes in use
```
**2,365,188个SNP、1056个样本，跟论文声称的数字完全精确吻合**——没有任何
SNP或样本在转换过程中丢失。运行时间442秒CPU(~7.4分钟)，远比29M_3k那次
(4.4小时)快，符合"这份VCF量级小得多"的预期。

**转换过程中出现过的警告(已确认无害)**：
```
all individuals set ignore.  Likely input problem (col 6)
resetting all individuals...
```
这是PLINK`.fam`文件第6列(phenotype)默认值是`-9`(未指定表型)触发的
convertf已知解析行为——convertf检测到"全部样本都被判定为ignore"这种
明显异常的情况，判断是输入问题而不是用户真的想丢弃全部样本，自动重置
恢复。最终`indivs: 1056`是全量，没有样本被误删，`numvalidind: 1056`
也确认了这一点。**这条警告本身不需要处理**，但以后遇到同样警告时，
务必像这次一样核对最终输出的`indivs`数字是不是符合预期全量，不能只看
警告文字就恐慌，也不能完全无视——核对数字才是关键。

**结论**：PCA-C的EIGENSTRAT输入(`civan_snp.eigenstratgeno/.snp/.ind`)
已经就绪。三个panel(29M_3k/6.7M_720/Civáň)现在**全部**有可直接喂给
`smartpca`的EIGENSTRAT格式数据了。

### 5.8 pseudo-haplotype调用+smartpca投影：四步流水线完整清单
（2026-08-12交接快照，⚠️2026-08-13深夜已整体跑通，见下方状态列更新
和📍节——本节结构保留作历史参照，状态描述已更新为当前真实情况）

这一节是专门为"新开一个Claude Code窗口接手"准备的，尽量做到不需要
回头翻整个对话历史。

**四步流水线（含2026-08-13新增的样本专属子集化+leave-one-out两步）
全部已在Civáň panel真实数据上跑通并验证通过**，16样本×3panel批量
铺开进行中：

| # | 脚本 | 仓库路径 | 状态 |
|---|---|---|---|
| ① | `map_besthit_to_irgsp.sh` | `scripts/ecotype_pca/` | **已完成**，16样本全部跑完，改用全基因组besthit读集（放弃了ORSC窄化，见`docs/ECOTYPE_PCA_PHASE1_COMMANDS.md`开头决策记录） |
| ② | `pseudo_haploid_call.py` | `scripts/ecotype_pca/` | **已完成并验证**，LV7008416379在Civáň panel上的调用结果（called=147）跟独立的`summarize_panel_overlap.py`普查预测完全一致。**2026-08-15：本脚本修了两轮真实bug（TV/ALL随机数流错位、call_rate指标拆分），见本文档最上方新增小节** |
| ②.5 | `build_sample_panel_subset.py`（**新增**） | `scripts/ecotype_pca/` | **已完成并验证**，把panel从百万级SNP瘦身到样本实际覆盖的几百个位点，`--mask-from`支持leave-one-out场景下"用真古样本的覆盖模式、取模拟样本的值"这种row-对齐需求 |
| ②.6 | `simulate_leaveoneout_projection.py`（**新增**） | `scripts/ecotype_pca/` | **已完成并验证**，leave-one-out正对照通过，REF/ALT方向确认无误 |
| ③ | `merge_ancient_into_panel.py` | `scripts/ecotype_pca/` | **已完成并验证** |
| ④ | `par.PROJECT.template` | `scripts/ecotype_pca/` | 模板仍保留占位符供参考；`run_sample_panel_pca.sh`（**新增**）会为每个样本×panel运行动态生成实际par文件和poplistname，不需要手工套模板 |
| 批量 | `run_sample_panel_pca.sh` + `summarize_projection_distances.py`（**新增**） | `scripts/ecotype_pca/` | 单条命令跑完①②③④全链路+产出"最近/次近现代群体"报告；16×3=48组合批量铺开**进行中** |

**①做什么**：把besthit过滤后的FASTQ(还没有位置信息、也还没去重)比对到
`irgsp.fa`，产出带位置信息、已用`samtools markdup`标记重复(不删除，②
在pileup阶段按需过滤)的BAM。三个panel共用同一份BAM，只映射一次。
`bwa aln`参数用的是文献里标准的古DNA设置(关闭seeding、放宽错配容忍度)，
**已核对与本项目主流程`scripts/server_originals/mapping.sh`完全一致**
(2026-08-12同日稍晚核实，`-l 1024 -n 0.01 -o 2`两边逐字相同，之前"未核对"
的说法已过时)；去重链路(`samtools collate|fixmate -m|sort|markdup`)和
过滤(`-F 0x904`)也已改为与主流程一致，唯一刻意保留的差异是不加`-r`
(不物理删除重复，只打标记)，理由见执行计划2.2节。

**①的运行命令(2026-08-12已给用户，结果未知；⚠️脚本同日稍晚已更新，
若之前已下载过旧版，需要重新curl覆盖)**：
```bash
cd /home/scratch/yinmt202607/gene/scripts
curl -O https://raw.githubusercontent.com/Inmpain/rice_adna_pipeline/codex/ecotype-pca-panel/scripts/ecotype_pca/map_besthit_to_irgsp.sh
chmod +x map_besthit_to_irgsp.sh

BESTHIT_DIR=/home/scratch/yinmt202607/gene/results/oryza_competitive_mapping/besthit
IRGSP_FA=/home/scratch/yinmt202607/db/asian_rice_panel_index/irgsp.fa
OUT_DIR=/home/scratch/yinmt202607/gene/results/ecotype_pca/bam_irgsp

nohup bash map_besthit_to_irgsp.sh "$BESTHIT_DIR" "$IRGSP_FA" "$OUT_DIR" \
  LV6000619499 LV6000619917 LV6000620016 LV6000620032 \
  LV6000620166 LV6000620172 LV6000654686 LV6000654698 \
  LV7008416272 LV7008416280 LV7008416294 LV7008416329 \
  LV7008416339 LV7008416349 LV7008416379 LV7008416407 \
  > map_besthit_to_irgsp.log 2>&1 &

tail -80 map_besthit_to_irgsp.log
```
16个样本名抄自besthit分支`ORYZA_BESTHIT_HANDOFF.md`第6.3节。`OUT_DIR`是
新建目录，之前不存在。**每个样本跑完脚本会打印mapped/duplicates_flagged/
mapq30/mapq20四个数字，这些数字要看一眼是否合理(不应该是0，也不应该
接近besthit输出的原始read数——besthit那16个样本的kept_reads在2445-12924
之间，见besthit分支6.3节，比对之后应该接近这个量级，不应该差太多)。**

**②做什么**：给定一个样本的BAM和某个panel的`.snp`文件，只在样本reads
真实覆盖到的panel位点上做pseudo-haplotype调用(随机抽一条覆盖read的
碱基当基因型，永远是纯合0或2，不会是1/杂合)，输出行数/顺序跟panel的
`.snp`文件完全一致(没覆盖到的位点填9)。**两个关键方法学决定**：
- 默认排除transition位点(A/G、C/T)的调用(`--transversions-only`默认
  开启)，规避古DNA末端损伤(C→T/G→A)误判，因为本项目自己的损伤窗口
  校准还没定论(besthit分支7.5节)。**操作上应该对每个样本每个panel跑
  两次(TV默认版+`--no-transversions-only`的ALL版)，双轨对照，见执行
  计划第7节。**
- **⚠️REF/ALT哪一列是哪个方向，权威判断依据是leave-one-out模拟
  (执行计划第6节)，不是check_ref.py的FASTA匹配率**——已知720这个
  panel的匹配率只有91.5%不是100%(3.1节)，但这不能直接等同于0/2编码
  方向错误，两者是不同问题，见执行计划2.1节的详细区分。`--swap-ref-alt`
  只应该在leave-one-out模拟证实需要翻转之后才加。

**③做什么**：把②对多个古代样本的调用结果，作为额外的列，一次性拼接进
panel自己的`.eigenstratgeno`文件(不是分开建文件——smartpca的
`-lsqproject`机制要求现代和古代样本在同一份基因型文件里，靠`.ind`文件
最后一列的群体标签区分谁参与建轴、谁被投影，见④)。**2026-08-12同日
稍晚已加固**：行数不匹配现在是硬退出+清理临时文件(不会再有IndexError
崩溃或静默写出不完整文件的情况)，新增字符校验和重复ID检查，输出改为
原子写入，见执行计划2.3节。

**④做什么，以及为什么现在填不了**：`lsqproject: YES` + `poplistname`
只列现代群体标签(不含古代样本的占位标签，默认叫`Ancient`)，这样
smartpca只用`poplistname`里列出的现代样本计算特征向量，`Ancient`标签
的古代样本会被投影而不参与建轴。**`poplistname`文件的具体内容依赖真实
的现代群体标签已经合并进`.ind`文件——这件事三个panel都还没做，是待办
1/2/2c，也是现在唯一真正卡住往下走的事**，所以④现在只是个带占位符的
模板，故意没有编一份假的`poplistname`去填(填假的会让smartpca"看起来
能跑"但结果毫无意义，比报错更危险)。模板还应补充
`numchrom:12`/`numthreads:8`/`numoutlieriter:0`，见执行计划第9节。

**三个panel的服务器路径速查(截至2026-08-11晚全部确认过)**：
```
29M_3k:   /home/scratch/yinmt202607/db/29M_3k/
          NB_final_snp.eigenstratgeno / .snp / .ind
6.7M_720: /home/scratch/yinmt202607/db/6.7M_720/
          asn720.6m.geno / .ind / .snp
Civáň:    /home/scratch/yinmt202607/db/paper1/
          civan_snp.eigenstratgeno / .snp / .ind
```

**群体标签三项待办，具体核对命令（当前唯一的真实阻塞项，建议顺序
2→1→2c）**：

1. **29M_3k(待办2，建议最先做)**——直接用main分支
   `docs/references/3k_rice_genomes_project/rice_line_metadata_20141029.xlsx`
   里的IND/AUS/ARO/TRJ/TEJ/ADM标签按样本ID匹配进`NB_final_snp.ind`。
   这份元数据表描述的对象就是3K RG这3024份材料本身，理论上应该接近
   100%覆盖，不需要先核实覆盖率这一步，可以直接写匹配脚本。
2. **720(待办1)**——先核实`asn720data/asn720.pop.fam`的`OrA-OrF`标签
   能覆盖`asn720.6m.ind`里多少样本：
   ```bash
   awk '{print $2}' /home/scratch/yinmt202607/asn720data/asn720.pop.fam | sort > /tmp/pop_ids.txt
   awk '{print $1}' /home/scratch/yinmt202607/db/6.7M_720/asn720.6m.ind | sed 's/_merged$//' | sort > /tmp/ind_ids.txt
   comm -12 /tmp/pop_ids.txt /tmp/ind_ids.txt | wc -l    # 重叠数
   comm -23 /tmp/ind_ids.txt /tmp/pop_ids.txt | head -20  # asn720.6m.ind里匹配不上的样本，抽样看
   ```
   同时还要读`db/wild_rice_pangenome_README.txt`(从没人读过，见1.2节)
   确认`OrA-OrF`的精确定义。
3. **Civáň(待办2c)**——从`db/paper1/Table_S1.csv`按样本ID匹配栽培/
   野生/亚群标签进`civan_snp.ind`，注意这张表表头有大量尾随空列(见
   1.3节)，解析时要跳过或忽略这些空列，不要被列数吓到。

**四步流水线设计要点回顾(精华摘要，不用重读5.0-5.7节全部历史)**：
- ①bwa映射到irgsp.fa(三个panel共用一份BAM，已对齐主流程去重链路)
- ②pseudo-haplotype调用：只在真实覆盖到的位点调用，默认排除transition
  位点规避末端损伤，输出跟panel `.snp`文件行数/顺序完全一致，
  TV/ALL双轨、report字段已拆分no_coverage/allele_mismatch
- ③把②的结果作为额外列拼进panel自己的`.eigenstratgeno`(不是分开建
  文件)，已加固为硬检查+原子写入
- ④smartpca：`lsqproject: YES` + `poplistname`只列现代群体标签(不含
  `Ancient`)

**下一步顺序**：见`docs/ECOTYPE_PCA_EXECUTION_PLAN.md`第1节的10步
清单（Phase 0覆盖普查→群体标签→样本专属子集+leave-one-out冒烟→
铺开到全部16样本×3个panel→bootstrap点云增量设计）。
