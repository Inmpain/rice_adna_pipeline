# 亚洲栽培稻参考基因组讨论：完整对话（外部导入，供 besthit 分支参考）

> 本文档是用户与GPT关于 `asian_rice_panel.fa`(12个基因组: np7/mh63/X24_kas/
> azu/arc/liuxu + 7个G25_ruf_W*野生稻) 选取理由、Oryza属结构、参考基因组
> 质量评估标准、以及**taxid错标问题排查**的完整讨论记录。2026-08-08导入，
> 基本保留原文，未大改。
>
> **本次讨论直接产出的行动项**：确认`np7`(IRGSP/Nipponbare) 在
> `asian_rice_panel.acc2taxid`里被错标成4529(O. rufipogon)，应为4530
> (O. sativa)——修正脚本见`scripts/oryza_besthit/fix_asian_rice_panel_taxid.sh`，
> 详见`docs/ORYZA_BESTHIT_HANDOFF.md`相关章节。

> 导出说明：以下按当前会话顺序收录用户消息和助手正式回复，保留表格、链接、命令与终端输出；不包含系统指令、后台工具调用及检索日志。

---

## 用户

请育璐帮忙把这几个基因组打包，给薛者共享一下，名字也改为短名字.fasta。这些作为所有亚洲栽培稻可能reads来源的参考。我选取的标准是尽量cover亚洲栽培和野生稻的遗传多样性，存在一定冗余，考虑到你们的目标，冗余没关系，但是尽量要全：np7 mh63 X24_kas azu arc liuxu G25_ruf_W1214 G25_ruf_W0169 G25_ruf_W1750 G25_ruf_W3037 G25_ruf_W1536 G25_ruf_W1726 G25_ruf_W2064 。 为什么这个能覆盖遗传多样性？

## 助手

这套参考能覆盖较大遗传多样性，关键不在于“基因组数量多”，而在于它同时覆盖了栽培稻的主要亚群，以及普通野生稻中分化较深的多个谱系。

- 栽培稻部分覆盖主要遗传轴：若 `np7` 是 Nipponbare/MSU7、`X24_kas` 是 Kasalath，那么分别包括温带粳稻 Nipponbare、热带粳稻 Azucena、籼稻 MH63 和 Liu Xu、aus 型 Kasalath，以及 circum-Basmati 型 ARC10497。也就是说，indica–japonica 的主分化，以及 aus、Basmati 和籼/粳内部差异都有代表。[MH63](https://plants.ensembl.org/Oryza_sativa_mh63/Info/Index)、[Azucena](https://plants.ensembl.org/Oryza_sativa_azucena/Info/Annotation/)、[Liu Xu](https://plants.ensembl.org/Oryza_sativa_liuxu/Info/Index)、[ARC10497](https://plants.ensembl.org/Oryza_sativa_arc/Info/Index)、[Kasalath](https://academic.oup.com/dnaresearch/article/21/4/397/446522)。

- 野生稻部分不是随机挑选，而是跨系统发育支系和地理区域取样：  
  `W1726/W2064` 属 Or-Ib，`W0169/W1214` 属 Or-II，`W3037` 属 Or-IIIa，`W1536` 属 Or-IIIb，`W1750` 属 Or-unspecific；产地跨中国、泰国、菲律宾、孟加拉、印度和斯里兰卡。Or-IIIa 接近粳稻祖先，Or-I 与籼稻关系较近，而 Or-II、Or-IIIb 是与现生栽培稻距离更远的分支，因此能补充驯化过程中丢失或低频的序列、结构变异和抗性相关等位基因。[野生—栽培稻泛基因组研究及其 Supplementary Tables 1–2](https://pmc.ncbi.nlm.nih.gov/articles/PMC12176639/)。

- 多参考能降低单一 Nipponbare 参考造成的偏倚。该泛基因组研究发现，野生和栽培材料合计包含大量 Nipponbare 中不存在的序列及野生稻特异基因。因此，来自高度分化单倍型、插入序列或可变基因的 reads，在这套多参考中更容易被识别。

不过需要准确表述：这是“对主要遗传空间的代表性覆盖”，不是覆盖亚洲稻的全部单倍型。当前野生材料覆盖了论文划分的六个支系中的五个，缺少明确的 Or-Ia。Or-Ia 与 indica 接近且受籼稻基因流影响，MH63、Liu Xu 等能部分补偿，但如果目标真的是“尽量要全”，建议再加一个典型 Or-Ia 野生稻参考。

可以把理由概括成一句话：

> 该参考集在栽培稻层面覆盖了温带粳、热带粳、籼、aus 和 circum-Basmati 等主要亚群，在野生稻层面覆盖了多个深度分化的 O. rufipogon 支系及东亚、东南亚和南亚来源，因此能够同时代表驯化群体的主要遗传背景和野生祖先库中未被驯化保留的序列多样性。冗余会增加近缘参考之间的归属歧义，但有利于降低漏检和单一参考偏倚。

---

## 用户

他们的挑选的参考基因组从oryzagenus水平如何合理，比如oryza genus的组成状况，再说这个基因组的靠谱之处

## 助手

从 Oryza 属的尺度看，这套参考的合理性要加一个前提：它不是为了覆盖整个稻属，而是为了覆盖“亚洲栽培稻及其最近缘野生祖先”这一 AA 基因组谱系。在这个目标下选择是合理的；如果声称覆盖整个 Oryza genus，则明显不够。

## 1. Oryza 属的组成及这套材料的位置

不同分类体系对物种数略有差异。较新的基因组研究将 Oryza 描述为2个栽培种和约20–25个现存野生种，分属11种基因组类型：

| 层级 | 组成 | 与本项目关系 |
|---|---|---|
| 二倍体基因组 | AA、BB、CC、EE、FF、GG | 亚洲栽培稻属于 AA |
| 异源四倍体基因组 | BBCC、CCDD、HHJJ、HHKK、KKLL | 与亚洲栽培稻距离较远 |
| 栽培种 | 亚洲栽培稻 O. sativa；非洲栽培稻 O. glaberrima | 本项目针对前者 |
| O. sativa 的近缘祖先库 | 主要是 AA 型的 O. rufipogon species complex | 本套野生稻参考的来源 |
| 其他野生稻 | O. officinalis、O. australiensis、O. brachyantha、O. granulata 等 | 不是亚洲栽培稻的直接祖先，通常不作为阳性参考 |

近期稻属泛基因组研究列出了 AA、BB、CC、EE、FF、GG、BBCC、CCDD、HHJJ、HHKK 和 KKLL 这11种基因组类型，并显示非AA物种含有大量物种特异序列和较大的染色体结构差异。[Oryza 属泛基因组研究](https://www.nature.com/articles/s41467-024-54427-3)。

因此，这套参考可以定义为：

> Oryza 属内、以 AA 基因组支系为边界的亚洲栽培稻来源参考集。

这个边界是合理的，因为：

1. O. sativa 和 O. rufipogon 都是 AA 型二倍体，染色体组成相同、亲缘关系最近，后者是亚洲栽培稻的主要野生祖先库。
2. 亚洲栽培稻的 reads，绝大部分应由 O. sativa 的不同亚群解释；少数在现代栽培种中缺失或高度分化的序列，可由 O. rufipogon 的不同谱系补充。
3. BB、CC、EE、FF、GG及四倍体稻与 O. sativa 距离较远，不是一般亚洲栽培稻 reads 的直接来源。全部加入反而会增大保守基因、重复序列的多重比对问题。

换言之，选择 O. sativa + 多谱系 O. rufipogon，比机械地从每一种 Oryza genome type 各取一个，更符合“亚洲栽培稻 reads 来源”这个问题。

## 2. 栽培稻部分为什么合理

在假定 `np7` 是 Nipponbare/MSU7、`X24_kas` 是 Kasalath 的情况下，栽培稻材料覆盖了 O. sativa 的主要遗传分支：

| 参考 | 代表类型 |
|---|---|
| np7/Nipponbare | 温带粳稻，Geng-japonica temperate |
| Azucena | 热带粳稻，GJ-trop1 |
| MH63 | 混合型籼稻，Xian-indica-admixed |
| Liu Xu | XI-3B2 籼稻 |
| Kasalath | circum-Aus/aus |
| ARC10497 | circum-Basmati/aromatic |

因此它不仅覆盖最深的 indica–japonica 分化，还加入了 aus、circum-Basmati，以及粳稻、籼稻内部不同亚群。Azucena、ARC10497、Liu Xu 等本来就是按照3K水稻群体的 PCA/群体中心选择出来的代表材料，而不是凭品种知名度挑选。相关“platinum reference”项目用12套基因组代表了15个亚洲栽培稻亚群中的12个。[亚洲栽培稻群体参考基因组研究](https://www.nature.com/articles/s41597-020-0438-2)。

## 3. 野生稻部分为什么合理

所选7个 O. rufipogon 不是同一地区或同一支系的重复材料，而是覆盖多个系统发育分支：

| 材料 | 分支 | 地理来源 |
|---|---|---|
| W1726 | Or-Ib | 泰国 |
| W2064 | Or-Ib | 孟加拉 |
| W0169 | Or-II | 泰国 |
| W1214 | Or-II | 菲律宾 |
| W3037 | Or-IIIa | 中国 |
| W1536 | Or-IIIb | 斯里兰卡 |
| W1750 | Or-unspecific | 印度 |

其价值包括：

- Or-IIIa 接近粳稻祖先谱系；
- Or-I 类群与籼稻关系较近；
- Or-II、Or-IIIb 与现代栽培稻距离更远，可补充驯化瓶颈中丢失的序列；
- Or-unspecific 可覆盖不能被典型分支充分解释的混合或特殊遗传背景；
- 地理上同时覆盖东亚、东南亚和南亚。

这比从一个地区集中选多个 O. rufipogon 更能覆盖祖先单倍型、结构变异和存在—缺失变异。[野生—栽培稻泛基因组研究](https://www.nature.com/articles/s41586-025-08883-6)。

但需要指出一个边界：论文将 O. rufipogon 划分为 Or-Ia、Or-Ib、Or-II、Or-IIIa、Or-IIIb 和 Or-unspecific 六类，目前这7个材料覆盖了其中五类，缺少明确的 Or-Ia。Or-Ia 与 indica 很近，而且存在较多籼稻基因流，所以 MH63、Liu Xu 等可以部分补偿；若强调“野生稻谱系尽量完整”，最好再增加一个典型 Or-Ia。

## 4. 这些参考基因组技术上为什么靠谱

### 栽培稻参考

Azucena、ARC10497、Liu Xu 等“platinum reference”基因组具有以下特点：

- 使用超过100×覆盖的 PacBio 长读长数据；
- 使用 Illumina 数据进行碱基纠错；
- 使用 Bionano optical maps 独立验证染色体结构；
- 12套组装平均 contig N50 为23.10 Mb；
- 平均 BUSCO 完整度97.9%，重新校正后10套超过98%；
- 多数组装接近无缺口染色体级别。

Nipponbare 是长期人工校正的标准粳稻参考；MH63也已有高质量甚至无缺口版本。因此，前提是分享的确实是明确版本的 IRGSP/MSU、MH63RS 和 PSRefSeq 组装，而不是来源不明的旧文件。

### 所选野生稻参考

我核对了原论文 Supplementary Table 1。这7个材料全部采用 PacBio HiFi 组装，其指标范围为：

- 测序深度：18.7–24.4×；
- primary contig N50：约7.9–18.8 Mb；
- genome BUSCO：98.4%–98.7%；
- consensus QV：54.4–62.6；
- LAI：17.6–27.8；
- 均提供染色体级组装。

这些指标说明：

- QV超过50，碱基层面的随机错误率很低；
- BUSCO约98.5%，基因区完整度很高；
- Mb级 contig N50 足以可靠表示大多数结构变异；
- HiFi 对高度相似序列、杂合位点和复杂重复区的解析能力明显优于单纯短读长组装。

不过并非每个组装完全等质。W1750 的 LAI 约17.6，重复区连续性低于其余多数材料；W1214 和 W1750 的 alternate contigs 较多，反映较高杂合度。用于 reads 比对时建议以 primary/chromosome assembly 为主，不要把 primary 和全部 alternate haplotigs 无区别拼在一起，否则会增加重复比对。

另外，若 `X24_kas` 对应的是较早期 Kasalath 组装，其公开版本约覆盖91.1%的基因组、总长约330.6 Mb，质量可能低于其他长读长染色体级参考。因此需要先确认 `X24_kas` 的确切来源和版本；必要时可换成更新的 aus/circum-Aus 长读长参考。

## 5. 可以直接采用的表述

> 从 Oryza 属的系统发育结构看，本参考集并非试图覆盖稻属全部11种基因组类型，而是围绕亚洲栽培稻所在的 AA 基因组支系构建。栽培稻部分覆盖温带粳、热带粳、籼稻不同亚群、aus 和 circum-Basmati；野生稻部分覆盖 O. rufipogon 的多个深度分化谱系以及东亚、东南亚和南亚来源。该设计既代表了亚洲栽培稻现存的主要群体结构，又补充了其野生祖先库中因驯化瓶颈而丢失的序列和结构变异。所用参考主要为 PacBio/PacBio HiFi 长读长、染色体级组装，具有较高的 contig N50、BUSCO 完整度和共识序列质量，因此适合作为高灵敏度的 reads 来源参考集。

最后，改成短文件名时最好同时保留一个 manifest，记录原始品种名、群体/分支、组装 accession 和版本、primary/alternate 状态、是否含叶绿体和线粒体、原文件 MD5。否则 `.fasta` 本身可靠，但重命名后会失去最重要的来源可追溯性。

---

## 用户

在假定 `np7` 是 Nipponbare/MSU7、`X24_kas` 是 Kasalath 的情况下，栽培稻材料覆盖了 O. sativa 的主要遗传分支：

| 参考代表类型 | |
|---|---|
| np7/Nipponbare | 温带粳稻，Geng-japonica temperate |
| Azucena | 热带粳稻，GJ-trop1 |
| MH63 | 混合型籼稻，Xian-indica-admixed |
| Liu Xu | XI-3B2 籼稻 |
| Kasalath | circum-Aus/aus |
| ARC10497 | circum-Basmati/aromatic |

因此它不仅覆盖最深的 indica–japonica 分化，还加入了 aus、circum-Basmati，以及粳稻、籼稻内部不同亚群。Azucena、ARC10497、Liu Xu 等本来就是按照3K水稻群体的 PCA/群体中心选择出来的代表材料，而不是凭品种知名度挑选。相关“platinum reference”项目用12套基因组代表了15个亚洲栽培稻亚群中的12个。[亚洲栽培稻群体参考基因组研究](https://www.nature.com/articles/s41597-020-0438-2?utm_source=chatgpt.com)。 这个我没懂

## 助手

简单说，这段话想表达的是：亚洲栽培稻内部不是一个遗传上均一的群体，而是由若干差异明显的亚群组成。选这些品种，是为了给几个主要亚群各放一个“代表点”。

“按照3K水稻群体的 PCA/群体中心选择”可以这样理解：

- 研究人员比较3,000多个水稻品种的全基因组变异；
- 用 PCA 把每个品种画成遗传空间中的一个点；
- 遗传相似的品种会聚成一团，形成温带粳稻、热带粳稻、不同籼稻、aus、Basmati 等群体；
- 在每一团中，挑选最靠近该群体中心的品种，作为这个亚群的典型代表。

类似于要代表中国不同地区的人群，不是只选人数最多的地区，而是从不同地区分别选一个典型样本。

这6个材料分别承担不同作用：

| 参考 | 它主要代表什么 |
|---|---|
| Nipponbare | 温带粳稻，是经典粳稻参考 |
| Azucena | 热带粳稻，补充 Nipponbare 不能代表的热带粳稻变异 |
| MH63 | 籼稻，但具有一定混合遗传背景 |
| Liu Xu | 另一个明确的籼稻亚群 XI-3B2，补充籼稻内部差异 |
| Kasalath | aus/circum-Aus，与普通籼稻和粳稻都有明显差异 |
| ARC10497 | circum-Basmati，即香稻/Basmati 相关群体 |

所以，为什么不能只用 Nipponbare？

因为来自籼稻、aus 或 Basmati 的 reads，如果对应区域在 Nipponbare 中差异很大、发生结构变异，甚至完全缺失，就可能比对不上或错误比对。加入这些不同亚群的代表以后，相当于在亚洲栽培稻遗传空间的几个主要方向都放置了参考，能够减少只用一个粳稻参考造成的偏倚。

但我前面的表述有一点需要说得更准确：

> “12套基因组代表15个亚群中的12个”说的是原始 platinum reference 项目的完整参考集合，不是说目前列出的这6个材料覆盖了12个亚群。

目前这6个材料覆盖的是主要大类和几个重要内部亚群，属于“较粗粒度但范围较广”的代表性覆盖；它们并没有覆盖所有细分亚群。例如，亚热带粳稻、另一类热带粳稻和若干籼稻亚群没有单独的代表。

更准确的说明可以改成：

> 亚洲栽培稻由温带粳、热带粳、多个籼稻亚群、circum-Aus 和 circum-Basmati 等遗传群体组成。本参考集分别选择 Nipponbare、Azucena、MH63、Liu Xu、Kasalath 和 ARC10497，代表其中几个主要且分化明显的群体。它不能穷尽亚洲栽培稻的全部亚群和单倍型，但可以覆盖其主要遗传分化方向，从而显著降低只使用 Nipponbare 或单一籼稻参考造成的比对偏倚。

其中还要保留两个前提：先确认 `np7` 确实是 Nipponbare/MSU7，`X24_kas` 确实是 Kasalath；否则不能仅凭文件简称判断材料身份。

---

## 用户

那这和我们用重测序数据还是pangenome还是refseq上面下载的，还是文章测序自己组装的染色体水平的基因组，有很大的差别吗，每种的偏好是如何的呢？

## 助手

差别很大，但这四个概念并不是同一层面的四种东西：

- 重测序数据是原始 reads；
- 染色体水平基因组是 reads 经过组装得到的参考序列；
- pangenome 是把多个基因组组织起来的方法或数据结构；
- RefSeq 是保存、标准化和发布基因组的数据库体系。

一套文章组装的染色体水平基因组，后来可能同时被收入 GenBank/RefSeq，也可能成为某个 pangenome 的组成部分。

## 四类资源的区别

| 资源 | 本质 | 优点 | 主要问题 | 在你们项目中的用途 |
|---|---|---|---|---|
| 重测序数据 | 未组装的短 reads，通常10–30× | 样本多，适合描述群体频率和选择代表材料 | 不是完整参考；难以恢复插入、重复区和大型结构变异 | 用来验证参考集是否覆盖群体，不宜直接作为比对参考 |
| 短读长重测序组装 | 用 Illumina reads 拼出的基因组 | 成本低，样本数量可很多 | 碎片化、缺重复区、结构变异不可靠；参考引导组装还会继承参考偏倚 | 一般不推荐进入核心参考集 |
| 文章发表的长读长染色体组装 | PacBio/HiFi/ONT组装，通常加Hi-C或光学图谱 | 序列完整，可表示插入、PAV和SV | “染色体水平”不等于完全正确；质量取决于数据和验证 | 最适合作为参考集的基本单元 |
| RefSeq/GenBank/GWH中的组装 | 数据库发布的某一组装版本 | accession明确、可追溯、容易复现 | 数据库身份不自动代表最高质量；可能版本较旧 | 同一组装有官方 accession 时优先从数据库下载 |
| 多线性参考/pangenome集合 | 多个完整基因组并列 | 简单、兼容常规比对软件，便于追踪命中哪个基因组 | 重复序列多，multi-mapping严重，索引较大 | 适合来源筛查和按群体归属 |
| 图 pangenome | 把共有序列和不同单倍型整合成图 | 降低单一参考偏倚，mapping和变异检测更全面 | 软件和结果解释更复杂，不容易把read唯一归到某个品种 | 适合最大化比对率和变异检测 |

## 1. 重测序数据：适合“选参考”，不适合“充当参考”

3K水稻项目等重测序数据最大的价值，是告诉我们：

- 亚洲栽培稻有哪些亚群；
- 每个品种在遗传空间中的位置；
- 哪些品种接近群体中心；
- 目前选的参考是否遗漏某一大片遗传空间。

但短读长重测序数据本身不能直接等同于参考基因组。原因是短 reads 很难重新拼出：

- 转座子和重复序列；
- 着丝粒和近着丝粒区域；
- 大型插入、倒位和拷贝数变异；
- 该品种特异、但 Nipponbare 中不存在的序列。

如果将短 reads 比对到 Nipponbare 后生成“consensus genome”，通常只是在 Nipponbare 框架上替换 SNP 和小 indel，不能真正恢复这个品种自身的结构和新增序列。

因此，对你们来说，重测序数据更适合作为“测试集”：

> 将大量不同亚群的重测序样本映射到候选参考集，检查每个样本是否至少存在一个高匹配参考，并比较增加某个参考后能否显著提升 mapping rate 或 k-mer coverage。

## 2. Pangenome：质量取决于组成它的基因组

“pangenome”不自动意味着高质量。它可能是：

1. 从几千个短读长样本推断出来的基因集合；
2. 多个线性组装的简单集合；
3. 由多个高质量染色体组装构建的图基因组；
4. 只包含基因区域的 gene pangenome；
5. 覆盖全基因组序列和结构变异的 graph pangenome。

如果目标是 DNA reads 比对，应该优先选择：

> 由多个长读长、染色体水平、群体代表性基因组构建的全基因组 pangenome。

不应只使用“gene pangenome”，因为它可能只有基因/CDS，缺少基因间区、重复序列和结构变异序列。

图 pangenome 对提高 mapping recall 和变异检测最有优势，但若你们的核心问题是“这个 read 更可能来自哪个参考材料或哪个亚群”，多个线性参考进行 competitive mapping 往往更容易解释。图上一个 read 可能匹配多个共享节点，不一定能直接归属到具体品种。

## 3. RefSeq：优势是可追溯，不是天然比文章组装更好

需要区分：

- `GCF_...`：通常表示 NCBI RefSeq assembly；
- `GCA_...`：GenBank assembly；
- `GWH...`：Genome Warehouse assembly；
- 文章补充网站中的 FASTA：可能与数据库版本相同，也可能是作者更新版。

RefSeq的优势主要是：

- accession和版本明确；
- 文件格式较统一；
- 有基本质量和污染检查；
- 后续可以准确复现；
- 不容易出现“这个 fasta 到底是哪一版”的问题。

但RefSeq不保证一定比最新文章组装更完整。某个刚发表的 HiFi/T2T 组装可能比RefSeq中的旧版本好，只是还没有进入RefSeq。植物中也有不少高质量组装只有 `GCA_` 或 GWH accession。

因此正确偏好不是“RefSeq永远优先”，而是：

> 先选择生物学上合适、技术指标高的组装；如果同一个组装有正式数据库 accession，则优先从正式数据库下载并锁定版本。

## 4. “染色体水平”也分质量高低

仅仅标注 chromosome-level 还不够。它可能是：

- HiFi contigs直接接近整条染色体；
- PacBio contigs经过Hi-C或Bionano验证；
- 很碎的短读长 contigs，依靠Nipponbare共线性强行排列成12条伪染色体。

三者都可能写成“染色体水平”，但可靠性差别很大。

建议重点检查：

| 指标 | 理想情况 |
|---|---|
| 测序技术 | PacBio HiFi优先，其次为充分纠错的PacBio/ONT |
| 组装层级 | 12条主要染色体，另有少量 unplaced contigs |
| QV | ≥40较好，≥50很好 |
| Genome BUSCO | ≥95%可用，≥98%很好 |
| Contig N50 | 水稻中最好达到约10 Mb以上 |
| LAI | ≥20通常表示重复区达到较好的参考质量 |
| 结构验证 | 有Hi-C、Bionano、遗传图谱或独立共线性验证 |
| 污染检查 | 无明显微生物、接头或其他植物序列 |
| 单倍型处理 | primary assembly与alternate haplotigs清楚分开 |

你们挑选的7个 G25 O. rufipogon 材料属于 HiFi 染色体水平组装，BUSCO约98.4%–98.7%、QV约54–63，整体上很适合作为参考。[野生—栽培稻泛基因组研究](https://www.nature.com/articles/s41586-025-08883-6)。

Azucena、ARC10497、Liu Xu 等参考使用超过100× PacBio数据、Illumina纠错和Bionano验证，属于质量较高的栽培稻群体参考。[Platinum reference研究](https://www.nature.com/articles/s41597-020-0438-2)。

## 针对你们目标的推荐偏好

你们的目标是“尽量覆盖所有亚洲栽培稻可能的 reads 来源，允许一定冗余”。建议采用三层设计：

### 第一层：参考序列主体

选择每个主要栽培稻亚群和 O. rufipogon 主要谱系的高质量长读长、染色体水平 primary assembly。

优先级大致是：

1. HiFi/T2T + Hi-C或光学图谱验证；
2. 高覆盖PacBio/ONT + 短读长纠错 + 结构验证；
3. 普通长读长染色体组装；
4. 短读长 de novo assembly；
5. 基于Nipponbare的参考引导 consensus。

### 第二层：组织方式

如果更重视“read能否找到匹配来源”：

- 优先图 pangenome，reference bias最低。

如果更重视“read命中了哪个品种/亚群”：

- 使用多个线性参考做 competitive mapping；
- 对相同分数的多重比对，归到共同的亚群或 `O. sativa/O. rufipogon complex`，不要强行归到具体品种。

你们现在“多个短名字 `.fasta`”的做法更接近第二种，操作简单，也适合来源筛查。

### 第三层：用重测序数据验证覆盖度

用3K水稻或其他群体重测序数据测试：

- 每个亚群的mapping rate；
- unmapped reads比例；
- 每个样本到最近参考的遗传距离；
- 加入新参考后，新增的可解释 k-mers/read数量；
- 新增参考的边际收益是否趋于饱和。

这才是证明“参考集足够覆盖”的最好证据。

## 最实用的结论

对你们这个任务，最合适的不是在“重测序、pangenome、RefSeq、文章组装”中四选一，而是：

> 以论文发表并有正式 accession 的高质量长读长染色体组装作为基本参考；按主要栽培稻亚群和野生稻谱系构建多参考或图 pangenome；再用大规模重测序数据检验其群体覆盖度。

分享时建议主参考只放 primary chromosomes 和必要的 unplaced contigs。Alternate haplotigs、叶绿体、线粒体最好分别标记，不要全部无区别拼接。文件即使改成短名，也应附带 manifest 和 MD5，否则后续很难判断具体用了哪个组装版本。

---

## 用户

GWH 是什么？还有，你就比如The barley pan-genome reveals the hidden legacy of mutation breeding这个的pangenome我们下载下来了不论从捕捉沉积物内缘率非常低的reads还是后续的下游的分析，有什么用呢，还有短读长重测序组装 在这个方面的用

## 助手

## GWH 是什么？

GWH 是 Genome Warehouse（基因组仓库），隶属于中国国家基因组科学数据中心 NGDC/CNCB。它的功能类似 NCBI GenBank/RefSeq 和欧洲 ENA，用来保存：

- 基因组组装序列；
- BioProject、BioSample信息；
- 组装方法、测序深度和版本；
- GFF、CDS、蛋白等注释文件。

GWH accession 通常类似 `GWHESIZ00000000.1`。它说明该基因组有正式的数据来源和版本，但与 GenBank 一样，“存入 GWH”不等于一定是最高质量，仍需检查测序技术、组装水平、BUSCO、QV等指标。[GWH官方说明](https://ngdc.cncb.ac.cn/gwh/documents)。

## 这套大麦 pangenome 实际包含什么？

[2020年的这篇论文](https://www.nature.com/articles/s41586-020-2947-8)选择了20个大麦基因型，包括地方品种、现代栽培种和1个野生大麦，用来代表全球大麦的主要群体结构。

这里可能下载到两类完全不同的文件：

1. 20套完整染色体级基因组  
   每套约4–5 Gb，合起来接近100 Gb，包含大量重复序列。

2. single-copy pan-genome  
   作者从20套基因组中提取非重复、单拷贝区域，得到约638.6 Mb的非冗余序列，共约147万个序列簇。其中约402.5 Mb为20个材料共有，约235.9 Mb具有存在/缺失或拷贝数差异。

这两种文件在沉积物低内源DNA分析中的用途完全不同。

## 对沉积物中极低内源率 reads 有什么用？

### 1. 提高对遗传分化 reads 的召回率

只用 Morex 一个参考时，如果古代大麦或野生大麦的某段序列：

- 与 Morex 差异较大；
- 在 Morex 中缺失；
- 属于其他品种特有的 PAV；
- 携带不同单倍型；

这些短 reads 可能无法比对。20个参考增加了不同单倍型和存在—缺失序列，因此有机会找回一部分单参考漏掉的 reads。

但提升主要发生在可变、低拷贝和基因富集区域，不会把0.01%的内源率神奇地提高到很高。

### 2. 全基因组集合灵敏，但假阳性风险高

大麦基因组约5 Gb，且大部分是转座子和重复序列。如果把20套完整基因组简单拼接：

- 搜索空间接近100 Gb；
- 大量 read 会多重比对；
- 小麦、黑麦、其他 Hordeum 和环境植物的保守序列可能误配到大麦；
- 微生物或低复杂度 reads 也可能偶然命中；
- 加入参考越多，单纯依靠“比对 read 数”判断存在性的风险越高。

因此不能把“20基因组上多比对到几十条 reads”直接解释成沉积物中存在大麦。

### 3. single-copy pan-genome 更适合严格确认

论文的 single-copy pan-genome 去除了大量重复区域，因此特别适合：

- 筛选较可信的大麦核基因组 reads；
- 计算独立位点覆盖广度；
- 建立大麦特异 k-mer；
- 区分真实低覆盖信号与重复序列误配；
- 检测部分 PAV 和群体信息。

缺点是它只包含全基因组约10%的单拷贝序列，因此灵敏度会低于完整基因组。

比较合适的思路是：

> 完整 pangenome 做高灵敏度初筛，single-copy pan-genome 和近缘物种竞争比对做严格确认。

同时应把大麦叶绿体基因组单独作为一套参考。叶绿体拷贝数高，在低内源样品中通常更容易检出；但最终最好再用多个独立核基因组位点确认。

## 推荐的低内源分析框架

| 阶段 | 参考资源 | 目的 |
|---|---|---|
| 初筛 | 20套完整基因组或图 pangenome | 尽量找回可能的大麦 reads |
| 严格确认 | single-copy pan-genome | 去掉重复序列造成的假命中 |
| 物种鉴定 | 大麦、小麦、黑麦及其他 Hordeum 竞争参考 | 防止把近缘谷物误判成大麦 |
| 高灵敏检测 | 大麦叶绿体参考 | 提高存在性检测能力 |
| 古DNA认证 | 长度、末端C→T/G→A、重复率、阴性对照 | 判断是否为真实古代DNA |
| 结果统计 | 独立位点数和覆盖广度 | 不只看总比对reads数 |

在极低内源率下，覆盖到20个互不相邻的单拷贝位点，通常比200条全部落在同一个转座子上的 reads 更可信。

## 后续分析能做到什么？

取决于真实大麦 reads 的数量：

| 下游目标 | 极低内源数据的可行性 |
|---|---|
| 判断是否存在大麦/Hordeum DNA | 可行，需竞争比对和古DNA认证 |
| 区分大麦与小麦、黑麦 | 有条件可行，需要物种特异位点 |
| 区分野生与栽培大麦 | reads足够覆盖诊断位点时可尝试 |
| 推断大致群体或地理亲缘 | 低覆盖下使用 genotype likelihood、投影 PCA 或诊断SNP |
| 精确鉴定具体品种 | 通常不可靠 |
| 调用全基因组 SNP | 极低覆盖时不宜按常规二倍体方法调用 |
| 判断PAV/CNV | 需要多个位点和足够深度，极低覆盖时只能做候选证据 |
| 检测大型倒位 | 基本不可行，除非靶向富集或覆盖度较高 |
| 功能基因分析 | 可对命中基因注释，但不能轻易根据“未检出”判断基因缺失 |

尤其要注意：“没有 reads 命中某基因”在低内源沉积物中通常只是覆盖不足，不能解释为古代材料缺少该基因。

## 短读长组装在这里有没有用？

有用，而且对于几十 bp 的古DNA reads，短读长组装并非完全不合适。

2020年这套大麦 pangenome 的大多数材料主要通过 Illumina paired-end、mate-pair、10X和Hi-C构建染色体级组装；Morex另外有PacBio CLR长读长组装。它们的优势是：

- 基因区和单拷贝区局部序列通常较可靠；
- 可以增加不同材料的SNP、短单倍型和部分PAV；
- 对短古DNA read而言，局部碱基序列是否存在比整条染色体是否连续更重要；
- 论文中5,602个大于5 kb的Barke相对Morex缺失中，有5,446个得到遗传定位支持，说明相当一部分PAV是可信的。

但短读长组装的弱点在大麦中也很明显：

- 重复序列容易塌缩；
- 可能把组装缺口误判为生物学缺失；
- CNV拷贝数不准确；
- 新插入和复杂基因簇容易遗漏；
- 倒位和易位断点不够精确；
- 抗病基因等串联重复区域尤其容易出问题。

因此它适合：

- 低内源 reads 初筛；
- 单拷贝序列确认；
- SNP和部分PAV标记；
- 群体亲缘分析。

不适合单独承担：

- 精确CNV；
- 重复区分析；
- 复杂结构变异；
- 参考间完整PAV比较。

2024年已经发表了包含76个野生和栽培大麦长读长组装的扩展泛基因组，结构变异和野生多样性覆盖都优于2020年的20基因组版本。[76个长读长大麦泛基因组](https://www.nature.com/articles/s41586-024-08187-1)。比较实际的做法是保留2020年的 single-copy pan-genome 用于严格筛查，同时用2024年的长读长组装补充完整基因组和结构变异参考。

---

## 用户

那我还有问题，现在是捕捉尽量多的遗传多样性用的这些基因组，那不考虑地域吗，还是说等弄地域变异的的panel的时候再说？还是有其他更好的方法，比如admixture

## 助手

地域不能完全不考虑，但在“尽量捕捉 reads”这个阶段，遗传多样性应当优先于地域标签；到“推断来源地区”阶段，再建立专门的地域变异 panel。最合理的是两层设计，而不是用同一套参考同时解决两个问题。

## 第一层：高灵敏度捕捉参考

目标是回答：

> 这些沉积物 reads 中有没有稻/大麦，以及尽量不要因为参考偏倚而漏掉？

选参考时优先级建议是：

1. 系统发育分支或主要遗传亚群；
2. PCA/遗传距离中的边缘和中心材料；
3. 高质量染色体组装；
4. 地域覆盖；
5. 新增参考能够带来多少新的 k-mer 或可比对 reads。

地域在这一层是“约束条件”，不是主要分类标准。因为：

- 相距很远的两个地区，材料可能因传播和育种而遗传上很接近；
- 同一地区可能同时存在籼、粳、aus或野生稻等高度分化材料；
- 古代作物传播、贸易和迁徙会打破现代地域与遗传结构的对应关系。

例如，大麦2020年泛基因组的20个材料是根据约19,778份材料的 PCA 选择的；而该PCA的主要驱动因素本身就包括地理来源、二棱/六棱和春性/冬性。因此，遗传空间选择已经间接包含部分地域信息，但作者仍检查了不同地理和生态类型是否被覆盖。[大麦泛基因组研究](https://www.nature.com/articles/s41586-020-2947-8)。

所以不是“不考虑地域”，而是：

> 先按遗传结构选，之后检查每个主要地理区是否都有代表；发现遗传上独特但没有覆盖的地区，再补材料。

## 第二层：地域来源判定 panel

目标变为：

> 已经确认是目标作物的 reads，更像哪个地域群体？

这时不能只依赖那十几套参考基因组。每个地区一两个参考远远不足，因为单个品种不能代表整个地区。应该使用具有明确采集地信息的大规模重测序数据，例如：

- 野生种和地方品种；
- 传统农家种；
- 现代育成品种单独处理；
- 每个地区尽量平衡取样；
- 有条件时加入古代基因组。

从这些样本中建立地域 SNP/PAV panel，筛选：

- 地区间高 FST 的位点；
- allele-frequency差异大的位点；
- 地区富集或私有等位基因；
- 可唯一比对的单拷贝区域；
- 在古DNA短片段中仍容易覆盖的位点；
- 多个位于不同染色体、低连锁的独立位点。

对于内源率极低的沉积物样品，地域判断应基于多个独立位点的总体似然，而不能根据一两个“地区特异”read直接归属。

## ADMIXTURE 有什么作用？

ADMIXTURE有用，但它不是选参考或定义地域的唯一方法。它回答的是：

> 一个个体的基因组可以由几个统计学祖源成分以什么比例解释？

例如某个材料可能表现为：

- 80%籼稻相关成分；
- 15% aus相关成分；
- 5%粳稻相关成分。

它可以帮助：

- 确认主要群体结构；
- 找到祖源成分较纯、适合作为“端点参考”的材料；
- 找到明显混合的过渡材料；
- 防止把admixed材料错误地当作某个地区的纯代表；
- 评估现代地域群体之间的基因流。

但ADMIXTURE成分不等于地理来源。某个“成分1”只是数学上的遗传成分，必须结合样本产地、谱系、历史和其他分析解释。

此外，ADMIXTURE结果容易受到以下因素影响：

- K值选择；
- 各地区样本数不平衡；
- SNP连锁；
- 现代育种材料大量重复；
- 缺失率和低覆盖；
- 参考面板本身没有覆盖真正祖先群体。

## 更好的组合方法

建议不要在 PCA 和 ADMIXTURE 之间二选一，而是组合使用：

| 方法 | 主要作用 |
|---|---|
| PCA | 看遗传空间、发现边缘材料和离群点 |
| ADMIXTURE | 看祖源成分和混合程度 |
| 系统发育树/IBS距离 | 确认亲缘和深度分支 |
| FST | 筛选地区间分化位点 |
| k-mer增益 | 判断新基因组是否补充新序列 |
| mapping增益 | 判断新增参考是否找回更多真实reads |
| 地理分层 | 防止某些区域完全没有代表 |
| 交叉验证 | 检验地域panel能否预测未参与建模的样本 |

一个实用的选择规则是：

1. 每个主要遗传亚群或野生谱系至少选一个高质量参考；
2. 在每个亚群内，优先选遗传距离较远的材料；
3. 检查东亚、东南亚、南亚等主要地区是否均有代表；
4. 选择ADMIXTURE祖源成分较纯的材料作为亚群“端点”；
5. 冗余允许时，加入少量典型admixed材料；
6. 计算每增加一个参考带来的新增 k-mer 和 mapping reads；
7. 当边际收益趋近饱和时停止增加。

要注意，admixed材料通常是已有单倍型的重组，不一定贡献很多新序列。因此，为“捕捉多样性”选参考时，遗传端点往往比大量中间混合材料更有价值；但如果某个混合群体经过长期独立演化并含有自己的PAV或特殊单倍型，就值得保留。

## 对你们项目最适合的架构

建议最终分成三套而不是一套：

- `capture_reference`：高质量、多亚群、多谱系基因组，用于最大化 reads 召回；
- `geographic_marker_panel`：来自大量有地理信息重测序样本的高信息 SNP/PAV，用于地域归属；
- `competitive_decoy_reference`：近缘作物、野生种和其他禾本科，用于排除错误比对。

因此，现在挑参考时需要考虑地域，但只需确保不存在明显地域空白；不必为了每个省或国家都加入一个完整基因组。精细的地域推断应该留给第二阶段的大样本重测序 panel，并结合 PCA、ADMIXTURE、FST 和基因型似然共同完成。

---

## 用户

### 第二层：组织方式

如果更重视“read能否找到匹配来源”：

- 优先图 pangenome，reference bias最低 你说的reference bias是什么意思

## 助手

`reference bias`（参考基因组偏倚）是指：

> reads 是否能够成功比对，以及后续识别出的等位基因和变异，会偏向于“与所选参考基因组更相似”的序列。

它并不是样品本身的生物学差异，而是由选择某一个参考造成的技术偏差。

例如只使用 Nipponbare 作为水稻参考：

- 粳稻 reads 与 Nipponbare 相似，容易比对；
- 籼稻、aus、Basmati reads 差异较多，比对率可能较低；
- 如果某段序列只存在于籼稻、但 Nipponbare 中完全没有，这些 reads 根本没有位置可比对；
- 最后可能错误地认为样品更接近粳稻，或者认为该区域不存在。

## 三种典型的参考偏倚

### 1. SNP和小 indel 导致的偏倚

假设参考序列是：

```text
参考：ACGTACGTACGT
样品：ACGTTCGTACGT
```

样品有一个真实 SNP。一般仍可比对，但得分会稍低。

如果一个很短的古DNA read 同时具有：

- 几个真实群体变异；
- 一个小 indel；
- 末端C→T损伤；
- 测序错误；

它可能超过比对软件允许的错配数，最终无法比对。这样保留下来的 reads 会偏向参考等位基因。

### 2. 参考中不存在的序列完全无法捕捉

假设某个籼稻或野生稻具有一段20 kb插入，但 Nipponbare 没有：

```text
Nipponbare：A────────B
籼稻：      A──新增20 kb──B
```

来自这20 kb插入的 reads 在 Nipponbare 上没有对应位置，因此会：

- 完全不比对；
- 错误比对到相似重复序列；
- 被当作污染或非水稻 reads 丢弃。

这也是泛基因组最重要的价值之一：把单一参考缺失的序列也加入搜索空间。

### 3. 变异频率偏向参考等位基因

假设某个位点真实存在两种等位基因：

- A：与参考相同；
- G：与参考不同。

带A的 reads容易比对，带G的 reads因为附近还有其他变异而更容易丢失。最后观察到的数据可能是：

```text
真实：A 50%，G 50%
比对后：A 70%，G 30%
```

这样会影响：

- 杂合位点判断；
- allele frequency；
- PCA和ADMIXTURE；
- 群体亲缘关系；
- 古代与现代群体的遗传距离；
- 地域来源推断。

## 为什么图 pangenome 能降低这种偏倚？

单一线性参考只能提供一条路径：

```text
A → Nipponbare单倍型 → B
```

图 pangenome 可以同时提供多个路径：

```text
        → 粳稻单倍型 →
A  →    → 籼稻单倍型 →    → B
        → aus单倍型  →
```

read不必强行贴到 Nipponbare 路径上，而是可以选择与自己最相似的单倍型路径。因此：

- 含非参考 SNP/indel 的 reads 更容易比对；
- 品种特异插入序列有对应节点；
- alternative allele 不容易系统性丢失；
- 不同群体的 mapping rate 更均衡；
- 后续 SNP、PAV和群体分析受单一参考的影响更小。

## 在低内源沉积物数据中更明显

古DNA reads通常短、损伤多，而沉积物中目标DNA比例又低。一个50 bp read只要同时包含两三个群体变异和一个损伤位点，就可能无法比对到单一参考。

因此图 pangenome有机会找回单参考遗漏的真实 reads。但它也有代价：

- 多个单倍型高度相似，read可能无法唯一归属；
- 搜索空间增大，假阳性和多重比对增加；
- “成功比对”不等于能够判断具体品种或地域。

所以更准确的表述是：

> 图 pangenome更适合提高目标 reads 的召回率并降低单参考偏倚；但若要判断 reads 来自哪个品种、亚群或地区，还需要单拷贝序列、诊断SNP/PAV和近缘物种竞争比对。

例如，你们可以先用图 pangenome尽量捕捉水稻或大麦 reads，再将捕捉到的 reads 投射到统一坐标，并用地域 marker panel 做 PCA、ADMIXTURE或基因型似然分析。

---

## 用户

### 第三层：用重测序数据验证覆盖度

用3K水稻或其他群体重测序数据测试：

- 每个亚群的mapping rate；
- unmapped reads比例；
- 每个样本到最近参考的遗传距离；
- 加入新参考后，新增的可解释 k-mers/read数量；
- 新增参考的边际收益是否趋于饱和。

这才是证明“参考集足够覆盖”的最好证据。这是什么意思，还有np7的taxid应该是什么

## 助手

“用重测序数据验证覆盖度”就是：先不要凭参考基因组的名字判断它们是否足够，而是拿大量真实水稻样本做一次模拟考试，看这套参考能否公平地解释各亚群的 reads。

## 具体怎么验证

假设最开始只有 Nipponbare，一个籼稻样本可能只有较低比例的 reads 能高质量比对。依次加入 MH63、Kasalath、Azucena 等参考后，观察结果是否改善。

下面数字只是示意：

| 参考集 | 籼稻样本高质量 mapping rate |
|---|---:|
| Nipponbare | 72% |
| + MH63 | 89% |
| + Liu Xu | 92% |
| + Kasalath | 92.5% |
| + 第三个相近籼稻 | 92.6% |

前几个参考带来明显提升，后面每增加一个相近参考，只提高0.1%，说明参考集开始趋于饱和。

### 1. 每个亚群的 mapping rate

分别取温带粳、热带粳、籼稻不同亚群、aus和Basmati的重测序样本，映射到候选参考集。

重点不是总体平均值，而是检查：

- 是否所有亚群都能较好比对；
- 某个亚群是否系统性低于其他亚群；
- 加入该亚群参考后是否明显改善；
- 高MAPQ、唯一比对和覆盖广度是否同时提高。

如果粳稻95%，但aus只有75%，说明参考集对aus覆盖不足。

### 2. unmapped reads比例

检查没有比对上的 reads 是什么：

- 低质量序列；
- 接头或低复杂度序列；
- 叶绿体、线粒体；
- 微生物污染；
- 目标基因组中存在、但当前参考缺失的序列；
- 与当前参考差异过大的真实水稻序列。

不能把所有 unmapped reads 都解释为“遗漏的遗传多样性”。应该先去污染、去低质量，再判断加入新参考后是否有一部分原来 unmapped 的 reads 变成了高质量、合理分布的水稻比对。

### 3. 每个样本到最近参考的遗传距离

将3K水稻样本和候选参考一起做 PCA、IBS距离或 k-mer距离。

对每个重测序样本，寻找遗传上距离最近的参考：

```text
样本A → 最近参考 Nipponbare，距离较小
样本B → 最近参考 MH63，距离较小
样本C → 所有参考都很远
```

样本C所在群体就是潜在缺口。如果很多 circum-Basmati 样本都离所有参考很远，就需要加入 ARC10497之类的参考。

目标不是每个样本都有完全相同的参考，而是每个主要亚群都至少有一个距离足够近的参考。

### 4. 新增的可解释 k-mers/read

比如只用 Nipponbare时，某个籼稻样本中有100万个高质量 k-mers 无法在参考中找到。

加入MH63后：

- 其中30万个能在MH63中找到；
- 且位于正常的水稻染色体区域；
- 说明MH63补充了Nipponbare缺少的序列。

再加入一个非常接近MH63的籼稻参考，只新增5,000个 k-mers，说明这个参考的边际贡献较小。

这里“可解释”不是“随便匹配上”，而是要求：

- 匹配质量高；
- 尽量是唯一或低拷贝序列；
- 排除微生物和近缘植物污染；
- 分布在多个独立基因组区域；
- 不能全部来自同一个转座子。

### 5. 边际收益是否饱和

把参考一个一个加入，记录新增的：

- 高质量 mapped reads；
- unique k-mers；
- 被覆盖的PAV；
- 新增单拷贝序列；
- 原来距离很远、现在有近邻参考的样本数量。

如果曲线逐渐变平：

```text
参考数增加：  1   2   3   4   5   6   7
覆盖率：     72  88  93  95  96  96.2 96.3%
```

说明主要遗传空间已经覆盖。最后几个参考可能仍有价值，但主要作用是冗余保险，而不是补充新的大类群。

实际分析中最好改变加入顺序或采用“每次加入增益最大的参考”的方法，因为饱和曲线会受加入顺序影响。

对于你们的沉积物古DNA场景，还可以把现代重测序 reads 截短到30–80 bp，并模拟末端损伤后再测试。这样比直接使用正常150 bp现代 reads更接近真实应用。

## np7 的 TaxID

如果 `np7` 确实是 Nipponbare 的 MSU Release 7 / Os-Nipponbare-Reference-IRGSP-1.0，那么建议使用：

```text
TaxID: 39947
Taxon: Oryza sativa Japonica Group
Cultivar: Nipponbare
Assembly: IRGSP-1.0
Assembly accession: GCA_001433935.1 / GCF_001433935.1
Annotation: MSU Release 7
```

NCBI将 Nipponbare/IRGSP-1.0 放在 `Oryza sativa Japonica Group`，对应 [NCBI:txid39947](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=39947&mode=info)。IRGSP-1.0的组装 accession 为 [GCA_001433935.1](https://www.ebi.ac.uk/ena/browser/view/GCA_001433935.1)。

另外两个可能用到的层级是：

| 层级 | TaxID | 何时使用 |
|---|---:|---|
| Oryza sativa | 4530 | 只做物种水平分类时 |
| O. sativa Japonica Group | 39947 | Nipponbare参考序列的推荐TaxID |
| temperate japonica subgroup | 1736657 | 仅当数据库明确要求亚群级标签时 |

`np7` 本身不是正式分类名，只是内部文件简称。“7”一般指MSU annotation Release 7，并不是Nipponbare的TaxID或一个独立基因组物种。

如果用于 Kraken/Centrifuge 一类分类数据库，建议给Nipponbare序列标 `39947`，然后在最终汇总时可以向上合并到 `O. sativa` 的 `4530`。由于不同 O. sativa 亚群之间大量 reads 完全相同，不要把共享 reads 强行解释为 Nipponbare或温带粳稻特异；细分到亚群需要诊断SNP/PAV或LCA规则。

---

## 用户

```text
grep "np" asian_rice_panel.acc2taxid

np7.Chr1 np7.Chr1 4529
np7.Chr2 np7.Chr2 4529
np7.Chr3 np7.Chr3 4529
np7.Chr4 np7.Chr4 4529
np7.Chr5 np7.Chr5 4529
np7.Chr6 np7.Chr6 4529
np7.Chr7 np7.Chr7 4529
np7.Chr8 np7.Chr8 4529
np7.Chr9 np7.Chr9 4529
np7.Chr10 np7.Chr10 4529
np7.Chr11 np7.Chr11 4529
np7.Chr12 np7.Chr12 4529
np7.ChrUn np7.ChrUn 4529
np7.ChrSy np7.ChrSy 4529
np7.ChrM np7.ChrM 4529
np7.ChrC np7.ChrC 4529
```

所以我觉得他这个应该是错了，可以抽一下他的染色体看和irgsp是不是一样的。但是我也不知道他用的这个np7是不是我们经常用的irgsp

## 助手

你的判断基本正确：`4529` 是普通野生稻 Oryza rufipogon 的 TaxID。若 `np7` 真是 Nipponbare/IRGSP-1.0，这个 acc2taxid 标注就是错的。下一步应先确认序列身份，再改 TaxID，不能只凭 `np7` 这个简称。

`ChrSy` 和 `ChrUn` 强烈提示它来自 MSU/RGAP Release 7，因为这两个伪分子名称出现在该数据集中。但最可靠的方法仍是比较染色体长度和实际序列。

假设组合参考是 `asian_rice_panel.fa`：

```bash
samtools faidx asian_rice_panel.fa
awk '$1 ~ /^np7\./ {print $1, $2}' asian_rice_panel.fa.fai
```

然后查看官方IRGSP参考：

```bash
samtools faidx irgsp.fa
awk '{print $1, $2}' irgsp.fa.fai
```

可以提取单条染色体并比较除header外的序列MD5：

```bash
samtools faidx asian_rice_panel.fa np7.Chr1 |
    sed '1d' | tr -d '\n' | md5sum

samtools faidx irgsp.fa Chr1 |
    sed '1d' | tr -d '\n' | md5sum
```

若12条染色体全部一致，就能确认是同一套IRGSP/Nipponbare序列。也可以用全基因组比对：

```bash
dnadiff -p np7_vs_irgsp IRGSP_nuclear.fasta np7_nuclear.fasta
```

如果确认是Nipponbare，最具体的官方TaxID是39947；若整个panel统一只报告到物种级，则可以使用4530。4529无论如何都不适用于Nipponbare。

---

## 用户

```text
samtools faidx irgsp_bt2idx.2.bt2
[E::fai_build_core] Format error, unexpected character at line 1
[faidx] Could not build fai index irgsp_bt2idx.2.bt2.fai
(/home/usr/yinmt/.local/mamba/snakemake) [yinmt@login01 asian_rice_panel_index]$
(/home/usr/yinmt/.local/mamba/snakemake) [yinmt@login01 asian_rice_panel_index]$ samtools faidx asian_rice_panel.fa
(/home/usr/yinmt/.local/mamba/snakemake) [yinmt@login01 asian_rice_panel_index]$ awk '$1 ~ /^np7\./ {print $1, $2}'

asian_rice_panel.fa.fai

(/home/usr/yinmt/.local/mamba/snakemake) [yinmt@login01 asian_rice_panel_index]$ ls
all_wgs_asian_irgsp.acc2taxid  asian_rice_panel.fa.rev.1.bt2l              irgsp_bt2idx.rev.2.bt2
asian_rice_panel.acc2taxid     asian_rice_panel.fa.rev.2.bt2l              irgsp.fa
asian_rice_panel.fa            asian_rice_panel.fa.sa                      irgsp.fa.amb
asian_rice_panel.fa.1.bt2l     blastdb                                     irgsp.fa.ann
asian_rice_panel.fa.2.bt2l     IRGSP-1.0_representative                    irgsp.fa.bwt
asian_rice_panel.fa.3.bt2l     IRGSP-1.0_representative_2026-02-05.tar.gz  irgsp.fa.fai
asian_rice_panel.fa.4.bt2l     irgsp.acc2taxid                             irgsp.fa.pac
asian_rice_panel.fa.amb        irgsp.asm10.mmi                             irgsp_lowcomplexity.bed
asian_rice_panel.fa.ann        irgsp_bt2idx.1.bt2                          irgsp_lowcomplexity.sorted.bed
asian_rice_panel.fa.bwt        irgsp_bt2idx.2.bt2                          irgsp_lowcomplexity.txt
asian_rice_panel.fa.fai        irgsp_bt2idx.3.bt2                          lowcomplexity_qc
asian_rice_panel.fa.gz         irgsp_bt2idx.4.bt2                          wgs_bt2l_size_and_memory.tsv
asian_rice_panel.fa.pac        irgsp_bt2idx.rev.1.bt2                      
(/home/usr/yinmt/.local/mamba/snakemake) [yinmt@login01 asian_rice_panel_index]$ grep "np" asian_rice_panel.fa
>np7.Chr1
>np7.Chr2
>np7.Chr3
>np7.Chr4
>np7.Chr5
>np7.Chr6
>np7.Chr7
>np7.Chr8
>np7.Chr9
>np7.Chr10
>np7.Chr11
>np7.Chr12
>np7.ChrUn
>np7.ChrSy
>np7.ChrM
>np7.ChrC
```

## 助手

报错原因很明确：`irgsp_bt2idx.2.bt2` 是 Bowtie2 的二进制索引，不是 FASTA，不能用 `samtools faidx`。你应该比较现成的 `irgsp.fa` 和 `asian_rice_panel.fa`；两者的 `.fai` 都已经存在。

先运行：

```bash
awk '$1 ~ /^np7\./ {print $1, $2}' asian_rice_panel.fa.fai
```

你刚才的 `awk` 命令没有把文件名放在同一条命令里，所以没有输出。

再查看IRGSP的序列名和长度：

```bash
awk '{print $1, $2}' irgsp.fa.fai
```

或者只看header：

```bash
grep '^>' irgsp.fa | head -20
```

如果IRGSP的header也是 `Chr1` 到 `Chr12`，可以直接比较长度：

```bash
diff -u \
  <(awk '$1 ~ /^np7\.Chr([1-9]|1[0-2])$/ {
      name=$1
      sub(/^np7\./, "", name)
      print name, $2
    }' asian_rice_panel.fa.fai | sort -V) \
  <(awk '$1 ~ /^Chr([1-9]|1[0-2])$/ {
      print $1, $2
    }' irgsp.fa.fai | sort -V)
```

没有输出表示12条染色体的名称归一化后长度完全一致。

## 做最可靠的逐染色体MD5比较

如果 `irgsp.fa` 的header是 `Chr1`、`Chr2` 等：

```bash
for c in {1..12}; do
    panel_md5=$(
        samtools faidx asian_rice_panel.fa "np7.Chr${c}" |
        sed '1d' |
        tr -d '\n\r' |
        md5sum |
        awk '{print $1}'
    )

    irgsp_md5=$(
        samtools faidx irgsp.fa "Chr${c}" |
        sed '1d' |
        tr -d '\n\r' |
        md5sum |
        awk '{print $1}'
    )

    if [[ "$panel_md5" == "$irgsp_md5" ]]; then
        status="IDENTICAL"
    else
        status="DIFFERENT"
    fi

    printf "Chr%s\t%s\t%s\t%s\n" \
        "$c" "$status" "$panel_md5" "$irgsp_md5"
done
```

如果12条全部显示 `IDENTICAL`，就可以确认 `np7.Chr1–Chr12` 与这个 `irgsp.fa` 的序列完全一致，不只是“相似”。

## ChrUn、ChrSy、ChrM和ChrC

这些额外序列不一定出现在你下载的 `irgsp.fa` 中：

- `ChrUn`：未定位序列；
- `ChrSy`：MSU/RGAP中的特殊未定位伪分子；
- `ChrM`：线粒体；
- `ChrC`：叶绿体。

因此先比较1–12号核染色体即可。`ChrUn + ChrSy` 的存在本身已经很像MSU Release 7数据包。

你也可以先检查IRGSP是否包含这些序列：

```bash
grep -E '^(ChrUn|ChrSy|ChrM|ChrC)[[:space:]]' irgsp.fa.fai
```

另外，检查FASTA header应使用：

```bash
grep '^>np7\.' asian_rice_panel.fa
```

比 `grep "np"` 更严格，避免将来误匹配序列内容。

一旦12条染色体MD5全部相同，`np7` 就可以确认是该IRGSP/Nipponbare参考，原来的 `4529` 必须修改：具体标注用 `39947`；如果整个分类库只保留物种级TaxID，则用 `4530`。

---

## 用户

```text
np7.Chr1 43270923
np7.Chr2 35937250
np7.Chr3 36413819
np7.Chr4 35502694
np7.Chr5 29958434
np7.Chr6 31248787
np7.Chr7 29697621
np7.Chr8 28443022
np7.Chr9 23012720
np7.Chr10 23207287
np7.Chr11 29021106
np7.Chr12 27531856
np7.ChrUn 633585
np7.ChrSy 592136
np7.ChrM 490520
np7.ChrC 134525

chr01 43270923
chr02 35937250
chr03 36413819
chr04 35502694
chr05 29958434
chr06 31248787
chr07 29697621
chr08 28443022
chr09 23012720
chr10 23207287
chr11 29021106
chr12 27531856
```

`asian_rice_panel.acc2taxid` 中 `np7.Chr1` 至 `np7.ChrC` 均标记为4529。所以原本的这个4529是个啥玩意？

## 助手

`4529` 是 NCBI Taxonomy 中的：

```text
Oryza rufipogon
普通野生稻
NCBI:txid4529
```

所以它不是Nipponbare的编号，而是被错误赋给 `np7` 的野生稻 TaxID。

你现在的结果已经提供了很强的证据：12条染色体长度全部精确一致，加上 `ChrUn`、`ChrSy` 的MSU7特征，基本可以判断：

> `np7` 是 Nipponbare MSU Release 7/IRGSP体系的序列，而不是 O. rufipogon。

如果要最终确认序列完全相同，可以适配 `chr01` 命名做MD5：

```bash
for i in $(seq 1 12); do
    panel_chr="np7.Chr${i}"
    irgsp_chr=$(printf 'chr%02d' "$i")

    panel_md5=$(
        samtools faidx asian_rice_panel.fa "$panel_chr" |
        sed '1d' | tr -d '\n\r' | md5sum | awk '{print $1}'
    )

    irgsp_md5=$(
        samtools faidx irgsp.fa "$irgsp_chr" |
        sed '1d' | tr -d '\n\r' | md5sum | awk '{print $1}'
    )

    [[ "$panel_md5" == "$irgsp_md5" ]] &&
        status="IDENTICAL" ||
        status="DIFFERENT"

    printf "%s\t%s\t%s\t%s\n" \
        "$panel_chr" "$status" "$panel_md5" "$irgsp_md5"
done
```

如果12条均为 `IDENTICAL`，就能确认 `np7.Chr1–Chr12` 与 `irgsp.fa` 完全相同。此时 `np7` 的4529标注确定是错的，应改为39947；若你们统一使用物种级TaxID，则改为4530。

---

## 用户

运行逐染色体MD5比较后，IRGSP一侧提示：

```text
[W::fai_get_val] Reference Chr1 not found in FASTA file, returning empty sequence
[faidx] Failed to fetch sequence in Chr1
Chr1 DIFFERENT db9bc24dfc44fdfa8a1e603465a35ca2 d41d8cd98f00b204e9800998ecf8427e
```

Chr2至Chr12同样提示 `Reference ChrN not found`，IRGSP一侧MD5均为 `d41d8cd98f00b204e9800998ecf8427e`。

## 助手

这不是序列不同，而是IRGSP的染色体名称写错了。

你的 `irgsp.fa.fai` 中是：

```text
chr01
chr02
...
chr12
```

脚本却查询：

```text
Chr1
Chr2
...
Chr12
```

FASTA名称区分大小写，而且01不能写成1。因此 `samtools` 没有取到IRGSP序列。你看到的：

```text
d41d8cd98f00b204e9800998ecf8427e
```

是空字符串的MD5，不代表IRGSP染色体。

先测试一条：

```bash
samtools faidx irgsp.fa chr01 | head
```

然后使用下面修正后的循环：

```bash
for c in $(seq 1 12); do
    panel_ref="np7.Chr${c}"
    irgsp_ref=$(printf 'chr%02d' "$c")

    panel_md5=$(
        samtools faidx asian_rice_panel.fa "$panel_ref" |
        sed '1d' |
        tr -d '\n\r' |
        tr '[:lower:]' '[:upper:]' |
        md5sum |
        awk '{print $1}'
    )

    irgsp_md5=$(
        samtools faidx irgsp.fa "$irgsp_ref" |
        sed '1d' |
        tr -d '\n\r' |
        tr '[:lower:]' '[:upper:]' |
        md5sum |
        awk '{print $1}'
    )

    if [[ "$panel_md5" == "$irgsp_md5" ]]; then
        status="IDENTICAL"
    else
        status="DIFFERENT"
    fi

    printf "%s\t%s\t%s\t%s\t%s\n" \
        "$panel_ref" "$irgsp_ref" "$status" \
        "$panel_md5" "$irgsp_md5"
done
```

预期输出类似：

```text
np7.Chr1  chr01  IDENTICAL  db9bc...  db9bc...
np7.Chr2  chr02  IDENTICAL  311f5...  311f5...
```

如果12条均为 `IDENTICAL`，就能确认 `np7.Chr1–Chr12` 与 `irgsp.fa` 完全相同。此时 `np7` 的4529标注确定是错的，应改为39947；若你们统一使用物种级TaxID，则改为4530。

---

## 用户

我们就到4530就行，不用再细了，然后同np7 mh63 X24_kas azu arc liuxu G25_ruf_W1214 G25_ruf_W0169 G25_ruf_W1750 G25_ruf_W3037 G25_ruf_W1536 G25_ruf_W1726 G25_ruf_W2064 。 这些的taxid，要是irgsp/np7都是4530的话，其余的taxid应该是什么？

## 助手

如果只分类到物种层级，那么规则很简单：

- 所有亚洲栽培稻 Oryza sativa：`4530`
- 所有普通野生稻 Oryza rufipogon：`4529`

对应关系如下：

| 参考前缀 | 材料 | TaxID |
|---|---|---:|
| `np7` | Nipponbare/IRGSP，O. sativa | 4530 |
| `mh63` | Minghui 63，O. sativa | 4530 |
| `X24_kas` | 如果确实是Kasalath，O. sativa | 4530 |
| `azu` | Azucena，O. sativa | 4530 |
| `arc` | ARC10497，O. sativa | 4530 |
| `liuxu` | Liu Xu，O. sativa | 4530 |
| `G25_ruf_W1214` | O. rufipogon | 4529 |
| `G25_ruf_W0169` | O. rufipogon | 4529 |
| `G25_ruf_W1750` | O. rufipogon | 4529 |
| `G25_ruf_W3037` | O. rufipogon | 4529 |
| `G25_ruf_W1536` | O. rufipogon | 4529 |
| `G25_ruf_W1726` | O. rufipogon | 4529 |
| `G25_ruf_W2064` | O. rufipogon | 4529 |

其中 `ChrM`、`ChrC` 等线粒体和叶绿体序列也沿用来源生物的TaxID，不需要给细胞器另设物种TaxID。

可以用下面的命令批量修正：

```bash
cp asian_rice_panel.acc2taxid \
   asian_rice_panel.acc2taxid.before_taxid_fix

awk 'BEGIN{FS=OFS="\t"}
{
    if ($1 ~ /^(np7|mh63|X24_kas|azu|arc|liuxu)\./) {
        $3=4530
    }
    else if ($1 ~ /^G25_ruf_(W1214|W0169|W1750|W3037|W1536|W1726|W2064)\./) {
        $3=4529
    }
    print
}' asian_rice_panel.acc2taxid \
> asian_rice_panel.acc2taxid.fixed
```

唯一仍需确认的是 `X24_kas`：如果序列确实属于栽培稻品种 Kasalath，就应为4530；若名称表示杂交材料或组合参考，则需回看原始FASTA来源。不过从目前的命名和使用背景看，按Kasalath/O. sativa标为4530是合理的。
