# 研究总纲：吴哥沉积物古DNA水稻判定——证据阶梯与实施路线图

> 本文档是横跨 `codex/oryza-competitive-mapping`(besthit) 和
> `codex/ecotype-pca-panel` 两条工作线之上的**项目级科研框架**，回答
> "我们到底要证明什么、现在能证明到哪一级、还缺什么"。两个分支各自的
> `ORYZA_BESTHIT_HANDOFF.md`/`ECOTYPE_PCA_PANEL.md`仍是各自技术细节的
> 权威文档，本文档不重复其中的参数/bug记录，只负责把两条线放进同一个
> 科学论证结构里。
>
> 内容整合自：用户与GPT讨论产出的科研框架草案 + 本仓库两条分支的实际
> 技术现状。2026-08-07首次写下。

---

## 0. 一句话定位

我们不是要回答"这是不是水稻"一个问题，而是要在**证据阶梯**上尽量往上爬：
从"是不是古DNA、是不是Oryza属"，到"栽培还是野生"，到"哪个栽培群/野生谱系"，
到"什么生态型(旱稻/水稻)"，每一级都需要独立的证据类型，**不能用低层级的
证据直接跳答高层级的问题**（比如不能因为PCA投影靠近某个野生簇，就直接说
"这是旱稻"——见第3节的关键澄清）。

---

## 1. 证据阶梯（回答问题的层级结构）

| 层级 | 问题 | 主要证据 | 本项目对应环节 |
|---|---|---|---|
| 0 | DNA是否真的古老且属于该层位？ | 片段长度、末端损伤、空白对照、层位重复 | 未系统做，见第4节工作线1 |
| 1 | 是否为Oryza属？ | 与其他物种竞争比对 | besthit分支：131数据库竞争mapping |
| 2 | 属内哪个主要谱系？ | 核基因组+叶绿体竞争比对 | besthit分支：目前只做二元Oryza/非Oryza，未分谱系 |
| 3 | 栽培(sativa)还是野生(ORSC)？ | 大量核SNP群体似然 | ecotype-pca-panel分支：PCA祖源投影 |
| 4 | 哪个栽培群或野生簇？ | 参考panel投影、f-statistics | ecotype-pca-panel分支：目前设计到这一级 |
| 5 | 什么生态型(旱稻/水稻/深水稻等)？ | 群体归属+功能位点+沉积环境证据 | **目前没有专门资源覆盖这一级**，见第3节 |

**当前项目实际进度**：besthit分支在做层级1（部分涉及2，因为竞争库本身覆盖
18个Oryza种，但KEEP判定目前只二元切成"3个目标种 vs 其他"，没有细分谱系）；
ecotype-pca-panel分支设计到层级3-4；**层级0(古DNA真实性QC)和层级5(生态型)
目前都还没有专门的资源/流程**，是两个真正的缺口。

---

## 2. 三类Panel资源——现状审计

### A. besthit竞争比对库（回答层级1-2）

- 路径：129个WGS真核生物shard(`/home/database/ref20250728/cph_euk/wgs_eukaryota.{1..129}.fas.gz`)
  + 亚洲水稻panel(`db/asian_rice_panel_index/asian_rice_panel.fa`)
  + IRGSP(`db/asian_rice_panel_index/irgsp_bt2idx`)，共131个数据库
- taxid映射：`db/asian_rice_panel_index/all_wgs_asian_irgsp.acc2taxid`——
  **这是三个来源合并后的统一映射表**，不是只统计WGS部分（见
  `docs/ORYZA_BESTHIT_HANDOFF.md`第3/6.1节）
- Oryza属18个已知种在此文件里都有覆盖，contig数从12到905不等（*O.
  longistaminata*905最多，*O. meridionalis*12最少），但**目前只有3个种
  (rufipogon/sativa/nivara)算作KEEP判定的"目标Oryza"**，其余15个种即使
  比对上也只进非Oryza审计表
- ⚠️GPT指出的缺口，我们认同：这批参考是NCBI WGS contig级别（零散片段），
  不是染色体级别；**但contig边界本身未必是古DNA短read(几十bp)误配的最大
  来源**——更大风险是物种间参考数量失衡、重复/近重复contig、错误物种标签。
  升级到染色体级参考时不能只是"换个文件"，还要做等量化、去冗余、诊断位点
  验证，否则参考数量失衡问题可能原样带到新库里

### B. PCA祖源panel（回答层级3-4，更准确的名字应是population-ancestry PCA panel）

| 资源 | 样本/SNP | 当前用途 | 状态 |
|---|---|---|---|
| `db/29M_3k/` | 3024份3K RG驯化稻，~2900万SNP，PLINK | 建立栽培稻坐标(IND/AUS/ARO/TRJ/TEJ/ADM) | convertf转EIGENSTRAT已在服务器验证跑通(29,635,224 SNP/3024样本无损) |
| `db/6.7M_720/` | 720份(以野生稻为主)，~670万SNP，EIGENSTRAT | 野生稻基因型 | 待与29M_3k求交集，mergeit未验证；样本ID两种风格混合(ERR开头+B0xx_merged) |
| `db/asn720data/` | 疑似同一批720份，仅9.5万SNP | **只用其标签列**：FID是`OrA-OrF`群体标签，键在ERR风格ID上 | 2026-08-07推翻"已弃用"的旧判断；标签覆盖率(尤其B0xx风格样本)未查 |

三者关系：29M_3k给栽培稻坐标，6.7M_720给野生稻基因型，asn720data(目前)
只贡献标签，三者必须先用**共享accession/原始run ID/样本元数据**建立可审计
的crosswalk，才能形成统一panel——不能假设"720"这个数字在三个文件里指的是
同一批个体。

### C. 尚未接入流程的高质量组装原料

- `db/3k/wild/`：140+个染色体级野生稻组装，**下载了但身份完全没确认过**
- *A pangenome reference of wild and cultivated rice*(Guo et al. 2025,
  Nature 642:662-671)：145个染色体级组装(129份*O. rufipogon*+16份*O.
  sativa*)，配套`wild_sample_info.csv`，谱系分组`Or-Ia/Or-Ib/Or-II/
  Or-IIIa/Or-IIIb/Or-unspecific`
- `db/wild_rice_pangenome_README.txt`：刚发现，**内容还没读**，很可能说明
  上面两项是不是同一个东西

这三项不是"第四个PCA panel"，而是**组装+元数据原料库**，要先做身份审计、
变异提取、坐标统一，才能变成正式可用的panel（对应besthit数据库或PCA panel
的扩充素材）。

---

## 3. 关键澄清：野生稻有三套互不相同的分类体系，不能互译

这是本次整理里最容易踩的坑，必须写清楚：

| 体系 | 标签 | 来源 | 关系 |
|---|---|---|---|
| `asn720data`自带 | `OrA/OrB/OrC/OrD/OrE/OrF`、`OrADM` | 本项目720号panel自己的聚类结果，出处/聚类脚本未知 | **来源未知**，不能假设和下面两套是一回事 |
| Huang 2012 + Guo 2025泛基因组 | `Or-Ia/Or-Ib/Or-II/Or-IIIa/Or-IIIb/Or-unspecific` | 两篇论文的谱系分组(2025年在2012年基础上细化+部分样本重新归类) | 体系内部可比较，Guo 2025提出Or-IIIa最接近japonica祖先、Or-Ia接近indica祖先、Or-Ib接近aus祖先 |
| Kim et al. 2016 | `W1-W6` + admixed | 独立的另一批聚类 | 不能仅凭编号映射到Or-Ia等；W1接近indica、W4接近aus、W6接近japonica(但可能是近期基因流而非直接祖先) |

**字母/数字编号只是各自研究坐标系里的簇名**，不同研究选择的样本、SNP、
过滤标准、聚类数都不同。在拿到`OrA-OrF`的accession元数据、或确认
`asn720data`和某篇论文同源之前，**不能把`OrB`翻译成`Or-Ia`或`W1`或某个
生态型**。这是`docs/ECOTYPE_PCA_PANEL.md`第3.2节待办第1条要先查清楚的事。

**另一个更根本的结构性问题**：旱稻(upland)和水稻(paddy)**不是独立的系统
发育分支**——indica和japonica内部都有旱稻品种。也就是说，就算完全查清楚
`OrA-OrF`/`Or-Ia-IIIb`/`W1-W6`这几套遗传谱系标签，**它们回答的是"这个样本
遗传上更接近哪个野生/栽培谱系"，不直接等于"这个样本是旱稻还是水稻"**。
真正的生态型判定需要一份独立的**耕作方式/生活史元数据**（栽培型：
upland/rainfed lowland/irrigated/deepwater；野生型：annual/perennial、
淹水时长等），这份元数据目前项目里完全没有，原计划靠SNP-Seek passport表
(已下线)，是一个**第四类资源缺口**，不是靠PCA panel能自动补上的。

---

## 4. 五条工作线

### 工作线1：古DNA真实性与样品单位（对应证据阶梯层级0，目前空白）

问题：16个古稻样本的reads是否真正古老、属于目标层位？一个沉积样品是否
混合了多个植株/时期？

需要补充：每个library的read数/长度分布/去重率/末端C→T/G→A损伤指标
（`results/02.irgsp/mapdamage_summary_bwa.tsv`理论上有这批数据，但main
分支文档记录这批BWA版mapDamage从未真正跑完过——见`docs/3krgp_integration_
and_simulation_prep.md`）、extraction/library blank对照、`angkor_robot_
library.txt`里的层位/年代模型细节。

### 工作线2：Oryza属与物种确认（对应层级1-2，besthit分支主体）

现有工具：131数据库+KEEP规则(见第2节A)。

优先补充：
1. 导出131个数据库的正式manifest(数据库名/accession/物种/taxid/contig数/
   总长度/N50/版本)——目前只有零散的contig计数，没有正式清单
2. `db/3k/wild/`与Guo 2025 pangenome做checksum/accession/sample ID对照，
   确认是否是同一批数据（见第2节C）
3. 评估建立平衡、去冗余的Pan-Oryza竞争参考集(用染色体级组装，控制每物种
   参考量，避免contig多的物种获得不公平优势)
4. 重新设计KEEP的分级标签(`Oryza-confirmed`/`target-AA-complex`/
   `O.sativa-like`/`ORSC-like`/`other-Oryza`)，而不是现在这种二元判定

### 工作线3：栽培/野生及遗传祖源（对应层级3-4，ecotype-pca-panel分支主体）

现有工具：见第2节B。

优先补充：
1. 建立`sample_crosswalk.tsv`：把29M_3k、6.7M_720、asn720data三边的样本
   ID(原始ID/ERR/SRR/B0xx/accession)、以及可能的Guo 2025谱系标签，统一
   对齐到一张表
2. 不只靠字符串匹配去重——`CX382`这类疑似重复(见`docs/ECOTYPE_PCA_PANEL.md`
   1.2节)最好用共享SNP做fingerprint/IBS确认是否为同一个体
3. 完成`29M_3k ∩ 6.7M_720`求交集(mergeit)，做downsampling测试古代样本
   低覆盖度投影偏差
4. 除PCA外，补充f3/f4、D-statistics等正式统计检验，给置信区间，不要只
   凭"投影位置靠近某簇"下结论

### 工作线4：生态型与耕作方式（对应层级5，目前完全空白，需要新建资源）

问题：是upland、rainfed lowland、irrigated、deepwater/floating，还是其他
适应型？

**当前结论：现有PCA panel回答不了这个问题**（见第3节）。需要另建一个独立
的`ecology_metadata_panel`：样本稳定ID(要能跟29M/6.7M的SNP数据对上)、
cultivation type、水深/淹水时长/播种方式/成熟期/光周期、annual/perennial、
数据来源与置信等级。可以用功能位点(SUB1A耐淹、SNORKEL1/2深水伸长、sh4/
PROG1/Rc/Bh4驯化位点)作第二证据层，但不能替代全基因组+环境证据。

### 工作线5：吴哥考古—水文综合解释（把遗传结果落到考古叙事）

需要补充：样本精确出土地点(寺庙居住层/粮食加工区/水库/壕沟/农田)、年代
及沉积速率、稻属花粉/植硅体/炭化籽粒等植物考古证据、吴哥及洞里萨湖周边
现代野生/栽培/杂草稻地理参考。

**新找到的直接相关文献**：Cobo Castillo et al. 2020, *The Khmer did not
live by rice alone: Archaeobotanical investigations at Angkor Wat and Ta
Prohm*, Archaeological Research in Asia 24:100213——**专门做吴哥窟和塔普伦
寺的植物考古学**，是目前唯一直接对应本项目遗址的文献，之前完全没查到过，
值得优先精读。

---

## 5. 分阶段实施计划

| Phase | 任务 | 交付物 | 通过标准 |
|---|---|---|---|
| **0：资源审计与冻结** | 生成全部数据库/panel/样本/标签的manifest，核对坐标版本、物种名、taxid | `reference_manifest.tsv`、`sample_crosswalk.tsv`、`panel_overlap_report.md`、`keep_taxonomy_rules.tsv` | 720是否独立个体、B0xx与ERR的关系、`db/3k/wild/`与Guo2025的关系都有可复查答案 |
| **1：物种竞争库升级** | 平衡去冗余的Pan-Oryza库、模拟短损伤reads测试、多级KEEP | 物种混淆矩阵、每物种precision/recall | 能量化O.sativa/ORSC/其他Oryza间的误分率 |
| **2：统一祖源Panel** | 合并29M_3k+6.7M_720、补完720标签、建东南亚子面板 | 固定SNP集、RAY投影、f-statistics | OrA-OrF每个标签都有明确来源和地理分布，对应不上就明写"不可互译" |
| **3：生态型元数据层** | 接入passport/论文补充表/本地采样信息 | `ecology_metadata_panel.tsv` | 元数据质量不足时，明确停在"遗传群"层，不制造生态型结论 |
| **4：吴哥综合解释** | 合并古DNA+年代+沉积学+宏遗存证据 | 每个古稻样本的证据卡(真实性/物种层级/祖源/生态型支持/置信等级) | 每条结论可追溯到数据层，"未检出"/"数据不足"/"支持排除"/"支持归属"严格分开 |

当前项目实际所处阶段：besthit分支在Phase 1（竞争库已建但未做manifest/
混淆矩阵），ecotype-pca-panel分支在Phase 0-2之间（坐标核对完成、求交集
未跑通、标签crosswalk未建），Phase 3/4完全未开始。

---

## 6. 当前立即优先级(P0)

1. **确认720号个体身份与标签覆盖**——解释`OrA-OrF`的第一阻塞点，见
   `docs/ECOTYPE_PCA_PANEL.md`3.2节第1条
2. **确认`db/3k/wild/`与Guo 2025 pangenome的关系**——避免重复建库，释放
   现成染色体级资源，见`docs/ORYZA_BESTHIT_HANDOFF.md`7.2节
3. **导出131数据库+KEEP的正式manifest**——先量化参考不均衡，再谈要不要
   替换/扩充参考
4. 【次优先】验证`29M_3k`∩`6.7M_720`的mergeit求交集
5. 【暂缓，等1-4有结论后再做】另建生态元数据层——不要继续把`OrA-OrF`或
   `IND/TRJ`直接当upland/deepwater标签用
6. 样品侧并行补齐古DNA QC(工作线1)——否则即使群体投影做得再漂亮，也证明
   不了它代表吴哥时期的真实信号

1-3都是**只需要在服务器上`cat`/`ls`/`grep`就能拿到答案**的信息核对工作，
不需要新写分析脚本，是当前性价比最高的下一步。

---

## 7. 与两个分支HANDOFF文档的关系

本文档不重复`ORYZA_BESTHIT_HANDOFF.md`(besthit技术细节、bug修复历史、
smoke test具体数据)和`ECOTYPE_PCA_PANEL.md`(convertf/mergeit踩坑记录、
CX382重复线索)里已有的内容，只负责把两条线的进度对应到统一的证据阶梯和
实施阶段上。**任何具体技术决策/参数/bug仍以两个分支各自的HANDOFF文档为
准**；本文档的P0清单里第1/2条如果查出结果，应该分别写回
`ECOTYPE_PCA_PANEL.md`和`ORYZA_BESTHIT_HANDOFF.md`，本文档只更新"现在在
哪个Phase"这个层面的总结。

---

## 8. 主要参考文献

- Huang et al. 2012, *A map of rice genome variation reveals the origin
  of cultivated rice*. Nature 490:497-501.
- Kim et al. 2016, *Population Dynamics Among six Major Groups of the
  Oryza rufipogon Species Complex*. Rice 9:56.
- Wang et al. 2018, *Genomic variation in 3,010 diverse accessions of
  Asian cultivated rice*. Nature 557:43-49.（即本项目`29M_3k`/`db/3k/`
  数据来源的3K RGP项目论文）
- Gutaker et al. 2020, *Genomic history and ecology of the geographic
  spread of rice*. Nature Plants 6:492-502.
- Cobo Castillo et al. 2020, *The Khmer did not live by rice alone:
  Archaeobotanical investigations at Angkor Wat and Ta Prohm*.
  Archaeological Research in Asia 24:100213.（**直接对应本项目遗址的
  考古植物学文献，工作线5优先精读**）
- Guo et al. 2025, *A pangenome reference of wild and cultivated rice*.
  Nature 642:662-671.（即本项目`wild_sample_info.csv`/145组装的来源）
