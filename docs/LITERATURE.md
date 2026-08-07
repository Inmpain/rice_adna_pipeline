# 文献追踪：我们关心的论文与用途

> 这份文档只做一件事：记录"我们为什么关心这篇文章、这篇文章能不能直接支撑
> 我们现在的分析"。原始数据文件（能下载到的）放在
> `docs/references/`，本文档只放引用信息+用途判断，不重复贴大段原文。
>
> 本文档只在 `main` 分支维护（跨分支共用的参考文献信息，不在各工作分支
> 重复），各工作分支的具体分析结论仍写在各自的 handoff 文档里，本文档
> 只负责"这篇文章说了什么、能不能用、边界在哪"。

最后更新：2026-08-08

---

## 0. ⏰ 待办提醒：明天优先读这两篇（ORSC相关）

1. **Genome evolution and diversity of wild and cultivated rice species**
   （见1.1节——⚠️目前只有标题，作者/期刊/年份/DOI 待从Zotero补充确认）
2. **Phenotypic Variation and the Impact of Admixture in the Oryza
   rufipogon Species Complex (ORSC)**（见1.2节，引用信息完整）

这两篇已经下载到Zotero。读完后请回来更新1.1节的完整引用信息，以及两篇
文章各自读后对panel/生态型判定的具体启发（如果和1.2节已经写的初步判断
不一样，直接改这份文档，不用另外开新文档）。

---

## 1. 待读 —— ORSC（Oryza rufipogon Species Complex）两篇

### 1.1 Genome evolution and diversity of wild and cultivated rice species

⚠️**引用信息不完整**：目前只有用户从Zotero给出的标题，没有作者/期刊/年份/
DOI——这里不猜测填充，避免记错。**读文档的人请在读这篇文章时把完整引用
信息（作者、期刊、卷期页、DOI/PMID）补进来**，同时补充这篇文章对我们panel
构建/物种鉴定/生态型判定分别有什么具体启发。

### 1.2 Phenotypic Variation and the Impact of Admixture in the Oryza rufipogon Species Complex (ORSC)

**引用**：Frontiers in Plant Science, 2022. DOI:
[10.3389/fpls.2022.787703](https://www.frontiersin.org/journals/plant-science/articles/10.3389/fpls.2022.787703/full)

**这篇文章实际提供什么**：
- 分析了240份ORSC材料的44个形态/生活史性状，222份在IRRI完成主要表型鉴定
- 使用已有的**GBS**（genotyping-by-sequencing，简化基因组测序，只测基因组
  中限制性酶切位点附近的一部分序列）数据，共113,739个SNP
- 根据表型区分**多年生型 / 一年生型 / 中间(混合)型**
- 结合群体结构、栽培稻基因渗入（admixture）和叶绿体单倍型
- 补充表给出：材料编号、来源地区、物种注释、表型分组、遗传亚群、SINE
  code、叶绿体单倍型

**关键限制——GBS ≠ 染色体级基因组组装**：GBS只测到限制性酶切位点附近的
稀疏序列用来拿SNP，**这篇文章不提供240套完整、染色体水平的FASTA组装**。

**能不能用这篇文章，逐项判断**：

| 用途 | 能不能用 |
|---|---|
| 判断一年生/多年生/中间型 | 可以 |
| 判断材料是否受到栽培稻渗入 | 可以 |
| 划分ORSC遗传亚群 | 可以 |
| 从中选择生态型代表 | 可以 |
| 下载完整染色体FASTA | 不可以 |
| 把GBS reads拼成可靠参考基因组 | 不建议 |
| 直接把GBS序列加入mapping panel | 不建议 |

**一个特别容易混淆的地方**：文章里的 `W1–W6` 是**六个ORSC遗传亚群的名称**，
不是样品编号。而我们`asian_rice_panel.fa`里的7个野生稻基因组编号
——`G25_ruf_W1214` `W0169` `W1750` `W3037` `W1536` `W1726` `W2064`——
是**具体材料编号**（来自Guo et al. 2025，见2.1节），不能因为都以`W`开头
就把`W1214`理解成"属于W1亚群"、`W2064`理解成"属于W2亚群"，需要查两篇
文章各自材料表里真正的对应关系（IRGC编号/采集坐标），才能确认某个具体
材料落在哪个ORSC遗传亚群、是哪种生态型。

**两篇文章（本节1.2 + 2.1的Guo 2025）应该怎么配合**：
1. 以Guo et al. 2025的129份染色体级野生稻组装作为候选池
2. 用Guo 2025论文的群体结构、地理坐标和遗传距离判断遗传覆盖
3. 用本文（1.2）定义生态维度：多年生型 / 一年生型 / 中间型 / 栽培稻高渗入型 /
   不同地理和生境
4. 通过共同的种质编号、IRGC编号或采集坐标，把两套数据连接起来
   ——**不能只凭`Wxxxx`名称做对应**
5. 检查我们当前7个W材料分别落在哪个生态/遗传类别；缺哪一类，再从
   Guo 2025的129套完整组装里针对性补充

**对我们`asian_rice_panel`野生稻部分的具体结论**：可以扩充，但候选基因组
应优先来自Guo et al. 2025的129套完整组装；本文（1.2）更适合作为"生态标签
体系"和选样标准，不是序列来源。建议覆盖：

| 生态/遗传类别 | 建议代表数 |
|---|---:|
| 多年生、低栽培稻渗入 | 2–3 |
| 一年生/*nivara*-like | 2–3 |
| 中间生态型 | 1–2 |
| 明显栽培稻渗入型 | 1–2 |
| 地理极端或独特遗传群 | 1–2 |

不一定需要下载129套全部加入线性panel。**下一步待办**：先判断当前7套
（`G25_ruf_W1214/W0169/W1750/W3037/W1536/W1726/W2064`）分别属于本文
定义的哪个遗传群/生态型/是否有栽培稻渗入；如果7套集中在同一群或同一
生态型，再从129套里有针对性补充（这一条也记在
`docs/RESEARCH_ROADMAP.md`和`codex/oryza-competitive-mapping`分支的
`docs/ORYZA_BESTHIT_HANDOFF.md`第7/8节里，三处保持一致）。

---

## 2. 已用于当前panel —— 野生-栽培稻pangenome

### 2.1 Guo et al. 2025, A pangenome reference of wild and cultivated rice

**引用**：Nature 642:662–671 (2025). PMID: 40240605

**这篇文章提供什么**：
- 129份 *O. rufipogon* 染色体水平组装 + 16份 *O. sativa* 染色体水平组装，
  合计**145套**野生—栽培稻pangenome参考
- 主要由PacBio HiFi构建，少数使用ONT
- 提供primary assembly，部分还有alternate haplotype assembly
- 提供FASTA、注释、VCF和图pangenome(graph pangenome)文件

**数据获取**：
- 论文：https://www.nature.com/articles/s41586-025-08883-6
- NGDC BioProject: [PRJCA024131](https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA024131)
- Figshare数据集：https://plus.figshare.com/articles/dataset/A_pangenome_reference_of_wild_and_cultivated_rice/25697817

**与我们panel的关系**：我们`asian_rice_panel.fa`里的7个`G25_ruf_W*`野生
稻基因组，就是从这篇文章的145套pangenome里选出来的（不是来自1.2节的
ORSC表型文章）。例如：
- `W1214`：GWH accession `GWHESIZ00000000`，染色体水平，HiFi约19.2×
- `W0169`：GWH accession `GWHESIE00000000`，染色体水平，HiFi约19.5×

也就是说，我们目前的7个`G25_ruf_W*`本身就来自更适合作为mapping
reference的完整基因组资源，而不是2022年ORSC表型文章里的GBS组装
——这也是1.2节强调"不能直接从ORSC表型文章下载完整基因组加入FASTA
panel"的原因：我们已经在用对的那篇文章的数据了，缺的是**生态型标签**，
1.2节的文章负责补这个标签，不是补序列。

**⚠️与`db/3k/wild/`的关系待确认**：`db/3k/wild/`下140+个
`{SampleID}.transfer.merge.chr.fasta`野生稻/近缘种染色体级组装，是否
与本文的145套pangenome同源/重叠，还没有核实——见
`docs/RESEARCH_ROADMAP.md`第2节C、`codex/oryza-competitive-mapping`
分支`docs/ORYZA_BESTHIT_HANDOFF.md`第7.2节，这是main分支和该分支共同
在追的同一个开放问题。

---

## 3. 3K Rice Genomes Project 系列（亚洲栽培稻遗传多样性背景）

三篇通常一起引用的奠基性文献：

1. (Ed.). (2014). *The 3,000 rice genomes project*. GigaScience, 3(1).
   https://doi.org/10.1186/2047-217x-3-7 (PMID: 24872877)
2. Li, J.-Y., Wang, J., & Zeigler, R. S. (2014). *The 3,000 rice genomes
   project: new opportunities and challenges for future rice research*.
   GigaScience, 3(1). https://doi.org/10.1186/2047-217x-3-8
   (PMID: 24872878)
3. Wang, W., Mauleon, R., Hu, Z., Chebotarov, D., Tai, S., Wu, Z., Li,
   M., Zheng, T., Fuentes, R. R., Zhang, F., Mansueto, L., Copetti, D.,
   Sanciangco, M., Palis, K. C., Xu, J., Sun, C., Fu, B., Zhang, H.,
   Gao, Y., et al. (2018). *Genomic variation in 3,010 diverse
   accessions of Asian cultivated rice*. Nature, 557(7703), 43–49.
   https://doi.org/10.1038/s41586-018-0063-9 (PMID: 29695866)

**与我们项目的关系**：3024份亚洲栽培稻的SNP矩阵是`db/3k/`（服务器）和
`db/29M_3k/`（`codex/ecotype-pca-panel`分支用，3024份材料29M biallelic
SNP，PLINK格式）里数据的直接来源，是品种命名规则（VARIETY_INDEX/NAME/
IRIS_ID体系）和群体结构背景的权威出处。

**原始项目文档已上传**（2026-08-08，从本机`/Volumes/SSD/Downloads/
3kgroup/`导入，见`docs/references/3k_rice_genomes_project/`，避免以后
要用的时候再到处找）：

| 文件 | 内容 |
|---|---|
| `readme.txt` | 官方文件清单说明 |
| `200001.txt` | GigaDB数据集本身的引用格式（区别于上面3篇论文引用） |
| `consortium_list.csv` | 3K项目联盟作者名单（CAAS/BGI/IRRI三方） |
| `rice_line_metadata_20141029.xlsx` | 3000份材料完整元数据，分IRRI/CAAS两个分表（431KB，二进制，git直接存的原始文件） |
| `seq_file_mapping_to_SRA.txt` | 每个原始测序文件到SRA run/experiment/sample/project accession的映射表（24MB，tab分隔） |

原始测序数据（3K项目本身，不是元数据）不在GigaDB托管，需要去SRA/ENA/DDBJ：
- Europe (ENA): http://www.ebi.ac.uk/ena/data/view/PRJEB6180
- USA (SRA): http://www.ncbi.nlm.nih.gov/sra/?term=PRJEB6180
- Asia (DDBJ): http://trace.ddbj.nig.ac.jp/DRASearch/study?acc=ERP005654

---

## 4. 文档索引

跨分支的完整文档目录见`file_path.md`第九节"全部文档索引"。
