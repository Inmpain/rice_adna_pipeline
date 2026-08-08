# 旱稻/水稻生态型 PCA 分析——数据库与设计草案

> 本文档记录一条与 `codex/oryza-competitive-mapping`(besthit)平行、下游依赖
> 其输出的新分析线：把 besthit 过滤后确认为 Oryza 的古代 reads，投影到现代
> 群体 PCA 空间，判断每个古代样本更接近旱稻还是水稻生态型。
> 首次写下: 2026-08-05

---

## 0. 现状一句话总结

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
图上。详见第2节新设计。

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

**⚠️当前整条线的硬性阻塞：besthit（`codex/oryza-competitive-mapping`分支）
还没有跑完**——这条PCA线的输入`<sample>.besthit_oryza.fastq.gz`要等besthit
产出，目前besthit那边还卡在acc2taxid taxid修正的`--apply`没跑（见该分支
`docs/ORYZA_BESTHIT_HANDOFF.md`第0.6节），16个古代样本里只有4个跑完全量
besthit。**这份文档目前只能做设计和脚本准备，不能真正跑通端到端**——
①A/①B/②A/②B可以先针对已经跑完besthit的4个样本(LV6000619499/619917/
620016/620032)写脚本、测试流程，但不要假设16个样本全部就绪。`③A`(3K RGP
官方标签，见3.2第2条)不依赖besthit，可以完全独立先做；`③B`(`OrA-OrF`
标签核实，见3.2第1条)同理。

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

### 1.3 两个panel的关系（历史设计，已废弃，存档）

~~分析设计里两者是互补而非替代关系：3K(29mio) 覆盖驯化稻的谱系多样性，
6.7M_720 补充野生稻/近缘种一侧——古稻样本理论上落在"驯化-野生"谱系的某个
位置，两个panel合起来才能给出有意义的PCA背景。~~ **2026-08-08已放弃**：
两个panel密度差23倍(chr1上305万 vs 78万)，真实交集只有一两千位点量级，
撑不起"合并成一个统一底盘"这个设计，改为两个独立PCA，见第0/2节。

## 2. 分析设计（2026-08-08修订版：两个独立PCA，不再求交集）

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
现在的方案牺牲了这一点，换来两条独立、密度更高(古代样本覆盖概率更高)的
路径。这个取舍是2026-08-08晚上根据mergeit实测交集量级（chr1上只有132个
真实重叠位点）做出的决定，不是从一开始就这样设计的——如果以后有更好的
方式统一坐标系(比如两边都投影到同一个更大的参考面板)，可以重新评估要不要
合并，但不是当前优先级。

**①A/①B"先subset到该样本实际覆盖的位点"这一步的意义，不只是逻辑上必须
(古代样本本来就不可能在panel全部位点上都有read覆盖)，还直接决定②/③两步
的计算量**——不管是29M_3k还是6.7M_720，②pseudo-haplotype调用和③smartpca
投影都只需要在"该样本reads实际覆盖到的那一小撮位点"上跑，不需要在panel
全部2900万/670万个位点上跑。16个古代样本、每个样本都远比整个panel稀疏，
**每个样本应该各自生成一份自己专属的、小得多的位点子集**，而不是对所有
样本都套用同一份大面板重复计算——这是2026-08-08用户明确要求写清楚的设计
原则，之前①A/①B的示意图里已经隐含了这一点，这里补一句显式说明，避免以后
实现的时候图省事直接对全量panel跑pseudo-haplotype/smartpca，浪费算力也
拖慢整条流程。具体实现上，这意味着②/③两步要做成"每个样本一份独立小面板"
的循环/并行结构，而不是"先固定一份面板，再套所有样本"的结构。

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
  smartpca只看自己panel内部的REF/ALT一致性，不涉及跨panel比较。

### 3.2 待办(按优先级，2026-08-08晚间因放弃合并而重排)

1. **`6.7M_720`独立PCA的群体标签来源**——确认`asn720data/asn720.pop.fam`
   的FID列(`OrA-OrF`)能覆盖多少`asn720.6m.ind`里的样本：
   a) 按IID(样本ID)把`asn720.pop.fam`的标签匹配到`asn720.6m.ind`上，
      算出能覆盖多少比例的720号样本(尤其`B0xx_merged`风格的样本有没有
      对应条目，目前只在`asn720.pop.fam`里见过`ERR`风格的ID)
   b) 读`db/wild_rice_pangenome_README.txt`(新发现，见1.2节)，确认
      `OrA-OrF`的精确定义，以及是否与`A pangenome reference of wild
      and cultivated rice`(Nature 2025)论文里的`Or-Ia/Or-Ib/Or-II/
      Or-IIIa/Or-IIIb`分组是**同一套体系还是两套不同的编号方案**——
      字母+罗马数字 vs 单字母，命名习惯不同，不能想当然认为是一回事
2. **`29M_3k`独立PCA的群体标签**——直接用3K RGP官方元数据
   (`docs/references/3k_rice_genomes_project/rice_line_metadata_20141029.xlsx`，
   main分支)里的IND/AUS/ARO/TRJ/TEJ/ADM标准亚群标签，**不需要再等
   `OrA-OrF`**，这条相对简单，可以先做。
3. **两条独立PCA各自的pseudo-haplotype调用脚本 + smartpca具体参数**（尤其
   `-lsqproject`相关配置）都还没写。PCA-A(29M_3k)现在就可以开始写，不用
   等标签问题解决——投影/建PC空间不需要标签，标签只在最后解读阶段才用得
   到；PCA-B(6.7M_720)同理，标签核实(待办1)可以跟脚本编写并行推进。
4. **核查`asn720.6m.ind`与`NB_final_snp.ind`之间的样本ID重叠**(即
   `CX382`疑似重复案例)——去掉`_merged`后缀后系统性比对两份样本列表，
   确认重叠规模。这原本是"合并两个panel前必须先去重"的问题，现在两个
   panel各自独立跑PCA，紧迫性降低（同一个体如果分别出现在两个panel里，
   不会互相污染对方的PCA结果），但如果最终要综合解读PCA-A和PCA-B的结果，
   知道两边有没有重复个体仍然有参考价值，优先级降到最后。
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
