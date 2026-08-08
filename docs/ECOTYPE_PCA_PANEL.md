# 旱稻/水稻生态型 PCA 分析——数据库与设计草案

> 本文档记录一条与 `codex/oryza-competitive-mapping`(besthit)平行、下游依赖
> 其输出的新分析线：把 besthit 过滤后确认为 Oryza 的古代 reads，投影到现代
> 群体 PCA 空间，判断每个古代样本更接近旱稻还是水稻生态型。
> 首次写下: 2026-08-05

---

## 0. 现状一句话总结

两个参考 panel 已经到位(见第1节实际路径)。坐标系核对已完成(见3.1)：染色体
命名一致、REF/ALT方向虽相反但mergeit可自动处理。convertf 这步(29M_3k
PLINK→EIGENSTRAT)已在服务器上验证跑通(2026-08-07，见3.2节step 2)：
29,635,224个SNP、3024个样本转换前后完全一致，无丢失。**⚠️2026-08-07重大
更正**：`asn720data`此前被判定"密度太低、直接弃用、不用查跟6.7M_720的关系"
——这个判断错了。`asn720data/asn720.pop.fam`的**FID列(第一列)是真实的群体
标签(`OrA`/`OrB`/`OrC`/`OrD`/`OrE`/`OrF`...)**，键在跟`asn720.6m.ind`同一套
`ERR0685xx`风格样本ID上——而`asn720.6m.ind`目前把这批样本的群体标签记成
笼统的`control`。这很可能就是第3.2节第3条一直找不到的旱稻/水稻标签来源，
比指望SNP-Seek或翻论文补充材料直接得多。详见1.2/3.2。

**⚠️2026-08-08晚间更新（mergeit已完成，但交集为0，找到了很可能的原因）**：
`par.MERGE`已在服务器上跑完（`##end of mergeit: 7329.262 seconds cpu
46421.983 Mbytes in use`），但输出的`Histogram of checkmatch return codes`
显示`total: 0`——**两个panel之间一个SNP都没匹配上**，`db/merged_29M3k_6M7_720/`
下产出的`merged.{geno,snp,ind}`预期是空的或接近空的（还没让用户贴
`wc -l`确认，见待办）。

**很可能的原因（未100%确认，见下面的验证步骤）**：`mergeit`按**SNP ID
字符串**匹配两个数据集里的"同一个SNP"，不是按染色体+物理位置。而这两个
panel的SNP ID命名方式完全不同：
- `29M_3k`(convertf转换后的`NB_final_snp.snp`)：ID是**纯数字、就是物理
  位置本身**，例如`1026`、`1033`、`1047`（跟第4列物理位置完全相同的数字，
  没有染色体前缀）——这也解释了convertf那次"`first snp 1026 is number`"
  提示的来源
- `6.7M_720`(`asn720.6m.snp`)：ID是`{染色体}np{物理位置}`格式，例如
  `1np1409`、`1np1422`

`"1026"`和`"1np1026"`作为字符串永远不可能相等，所以哪怕两个panel在同一条
染色体同一个物理位置真的都有SNP（大概率是有的——两边分别抽查了chr1的位点
数：29M有3,057,565个、720有775,775个，量级上完全应该有交集），mergeit按
ID字符串比对时也会判定"没有共同SNP"，跟convertf那次"文件名后缀识别"是
同一类"工具的默认假设跟我们数据的实际格式不匹配"的坑。

**这个诊断还没有100%坐实**，下一个session接手时应该先做两件事，而不是
直接冲去改ID重跑（可能要跑很久，蒙对了省时间，蒙错了浪费几个小时）：
1. **要用户贴mergeit的完整log（不只是尾巴），而不是只看`##end of mergeit`
   这几行**——看有没有更具体的"未匹配原因"提示，或者mergeit是不是其实
   有一个按位置匹配的开关/参数被漏掉了
2. 确认`merged.{geno,snp,ind}`到底是空文件还是有内容但行数为0/接近0

**如果诊断成立，标准修复思路**（不需要动`.geno`/`.ind`和原始文件，只是
把两份`.snp`文件的ID列临时统一命名后再合并一次）：
```bash
# 生成ID统一为"{染色体}_{物理位置}"格式的.snp副本，不改.geno/.ind、不改行序
awk 'BEGIN{OFS="\t"} {$1=$2"_"$4; print}' NB_final_snp.snp > NB_final_snp.idfix.snp
awk 'BEGIN{OFS="\t"} {$1=$2"_"$4; print}' \
  /home/scratch/yinmt202607/db/6.7M_720/asn720.6m.snp > asn720.6m.idfix.snp
# 把par.MERGE里的snp1/snp2改指向这两个idfix文件，其余(geno1/ind1/geno2/ind2)不变，重跑mergeit
```

当前最优先的事三条并行：①**（今晚新发现，提到最前面）确认mergeit
0匹配的真实原因、决定要不要用上面的ID重命名方案重跑**；②核实`asn720data`
标签能覆盖`asn720.6m.ind`里多少样本、读`db/wild_rice_pangenome_README.txt`
（新发现，内容还没看）搞清楚`OrA-OrF`具体定义；③注意：③号任务在另一个
会话窗口(疑似Codex CLI)同步推进
中，该文档`db/wild_rice_pangenome_README.txt`内容那部分是那边在等用户贴，
不要重复索取。

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

### 1.2 6.7M_720 panel(野生稻为主，已上传)

- 路径：`/home/scratch/yinmt202607/db/6.7M_720/`
  ```
  asn720.6m.geno
  asn720.6m.geno.gz    (压缩版，跟未压缩版并存)
  asn720.6m.ind
  asn720.6m.snp
  ```
- **格式是EIGENSTRAT**(`.geno/.ind/.snp`三件套)，不是PLINK——这正好是
  smartpca的原生输入格式，后续两个panel求交集后，最终喂给smartpca的应该
  统一转成这个格式，而不是PLINK。
- 内容：720份样本，约670万个SNP位点，以野生稻为主
- 来源/参考坐标系：**尚未确认是否也是 called vs Nipponbare MSU7/IRGSP1.0**——
  这是第3节里最优先要查的问题，不确认这一点，后面所有"求交集"的结果都
  可能是静默错误(代码能跑，位点对不上导致交集集合偏小或为空，但不会报错)。
- **⚠️`asn720data`不能弃用——它的`.fam`文件是群体标签的关键来源
  (2026-08-07更正，推翻此前"已决定弃用"的判断)**：
  `asn720data/asn720.pop.{bed,bim,fam}`(720份、94,974个SNP位点)本身的
  **基因型数据**密度确实太低，不适合拿来做PCA底层数据，这一点原判断没错；
  但`asn720.pop.fam`的**FID列(第一列)记录着`OrA`/`OrB`/`OrC`/`OrD`/`OrE`/
  `OrF`这样的群体标签**，例如：
  ```
  OrD    ERR068594    0    0    0    1
  OrD    ERR068597    0    0    0    1
  OrA    ERR068598    0    0    0    1
  OrF    ERR068600    0    0    0    1
  OrC    ERR068604    0    0    0    1
  ```
  这批`ERR0685xx`风格的样本ID，跟`asn720.6m.ind`里那批被记成笼统`control`
  标签的样本**是同一套ID体系**(见下面"样本ID观察")。也就是说：即使基因型
  数据只用`asn720.6m.*`(6.7M位点，密度更高)，**群体标签应该从
  `asn720data/asn720.pop.fam`按样本ID(IID)匹配过去**，而不是用
  `asn720.6m.ind`自己那个占位符式的`control`。**这批`OrA-OrF`标签，很可能
  就是第3.2节第3条一直在找的旱稻/水稻(或至少野生稻谱系)标签来源**——
  且与本工作流之前看到的一张smartpca投影图(图例里正好有`OrA`至`OrF`分类，
  以及`ADM/ARO/AUS/IND/RAY/TEJ/TRJ`等3K RGP标准亚群编号)用的是同一套
  `OrA-OrF`编号，强烈暗示两者同源，值得优先核实覆盖率而不是继续等
  SNP-Seek或啃论文补充材料。
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

### 1.3 两个panel的关系

分析设计里两者是**互补而非替代**关系：3K(29mio) 覆盖驯化稻的谱系多样性，
6.7M_720 补充野生稻/近缘种一侧——古稻样本理论上落在"驯化-野生"谱系的某个
位置，两个panel合起来才能给出有意义的PCA背景。**注意**：这个"互补"的前提
是两边样本互不重叠；1.2节记录的`CX382`疑似重复如果查实存在系统性重叠，
这个假设需要重新评估。

## 2. 分析设计(草案)

```
besthit 过滤后的 Oryza reads (来自 codex/oryza-competitive-mapping 分支产出:
  <sample>.besthit_oryza.fastq.gz / 或对应BAM)
        │
        ▼
① 位点合集 = 29M_3k(PLINK) ∩ 6.7M_720(EIGENSTRAT)
   两者格式不同，需要先统一格式再求交集(染色体命名已确认一致、都是裸数字
   1-12，见3.1；REF/ALT方向相反但mergeit能自动处理，见3.1)：
   a) 用 EIGENSOFT `convertf` 把 `29M_3k`(PLINK bed/bim/fam，用原始裸数字
      染色体版`NB_final_snp.bim.orig`)转成EIGENSTRAT(跟 6.7M_720 已经是的
      格式统一，也是smartpca的原生输入格式)
   b) 用 EIGENSOFT `mergeit` 在两个EIGENSTRAT数据集之间按染色体+物理位置
      求交集，自动纠正两边A1/A2顺序颠倒的问题，并剔除strand-ambiguous
      (A/T、C/G)位点，产出合并后的panel
        │
        ▼
② 对每个古代样本单独求交集：
   该样本最终使用的panel子集 = ① ∩ (该样本 besthit 后 reads 实际覆盖到的位点)
   —— 每个样本的panel子集不同，因为每个样本的read覆盖位置是随机、稀疏的
        │
        ▼
③ 在②确定的位点上，对每个古代样本做 pseudo-haplotype 调用
   (每个位点随机抽一条覆盖该位置的read，取其碱基作为该样本在该位点的等位型
   —— 标准古DNA pseudo-haploid做法，不做常规diploid genotype calling)
        │
        ▼
④ smartpca：用现代样本(3K+720子集，完整基因型)建PCA参考空间，
   古代pseudo-haplotype样本用 -lsqproject 投影模式投上去
        │
        ▼
⑤ 看古代样本在PC空间里落在哪个/哪些现代亚群附近，结合亚群的旱稻/水稻标签
   (标签来源见第3节——新找到候选来源`asn720data/asn720.pop.fam`)判断
   生态型归属
```

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
  暂不处理，等mergeit真正跑起来后对照它自己的匹配/翻转/丢弃日志再回头看
  这17个位点会不会落在"丢弃"名单里。
- **两个panel的REF/ALT方向确实相反，但这不阻塞mergeit**：EIGENSOFT的
  `mergeit`按坐标匹配后会自行比较两边的等位基因集合(不依赖谁被记成A1/A2)，
  能自动识别并纠正这种顺序颠倒，不需要手动交换720的两列。真正无法靠工具
  自动解决的只有A/T、C/G这类strand-ambiguous位点，标准做法是合并前直接
  剔除，见3.2。

### 3.2 待办(按优先级)

1. **重新打开：确认 `6.7M_720` 和 `asn720data` 的关系(2026-08-07重新打开，
   之前误判为"已决定弃用、不用查")**——现在明确知道`asn720data/
   asn720.pop.fam`的FID列是群体标签(`OrA-OrF`等)，需要：
   a) 按IID(样本ID)把`asn720.pop.fam`的标签匹配到`asn720.6m.ind`上，
      算出能覆盖多少比例的720号样本(尤其`B0xx_merged`风格的样本有没有
      对应条目，目前只在`asn720.pop.fam`里见过`ERR`风格的ID)
   b) 读`db/wild_rice_pangenome_README.txt`(新发现，见1.2节)，确认
      `OrA-OrF`的精确定义，以及是否与`A pangenome reference of wild
      and cultivated rice`(Nature 2025)论文里的`Or-Ia/Or-Ib/Or-II/
      Or-IIIa/Or-IIIb`分组是**同一套体系还是两套不同的编号方案**——
      字母+罗马数字 vs 单字母，命名习惯不同，不能想当然认为是一回事
2. **格式统一 + 求交集**——par文件和执行脚本已写好并提交到
   `scripts/ecotype_pca/`，服务器上直接跑：
   ```bash
   cd /path/to/rice_adna_pipeline/scripts/ecotype_pca
   bash run_convert_merge.sh
   ```
   脚本依次做三件事(对应下面a/b/c)：
   - `par.PLINK.EIGENSTRAT`：29M_3k 用原始(裸数字染色体)bim 转 EIGENSTRAT。
     **✅ 2026-08-07已在服务器验证跑通**，29,635,224个SNP、3024个样本转换
     前后完全一致。真实坑点记录(避免以后重踩)：convertf **没有**
     `inputformat:`这个参数(不存在，写了会被静默忽略)，它是**靠文件名
     后缀自动识别格式**——PACKEDPED格式要求snp文件以`.pedsnp`/`.map`/
     `.bim`结尾、indiv文件以`.pedind`/`.ped`结尾(见convertf自带README)。
     我们的原始文件是`NB_final_snp.bim.orig`/`NB_final_snp.fam`，后缀不
     符合，convertf识别不出PLINK格式，退回成EIGENSTRAT原生`.snp`格式的
     解析器去读(列顺序是ID在前、chrom在后，跟PLINK bim的chrom在前顺序
     相反)，导致把每个SNP的数字ID误当成染色体号，最终段错误。
     修复：脚本会自动建软链接`NB_final_snp.rawchrom.bim`→`.bim.orig`、
     `NB_final_snp.rawchrom.pedind`→`.fam`(不复制大文件)，par文件指向
     这两个软链接名。`.ind`输出里群体标签是`???`是预期行为(`.fam`第6列
     本来就没有群体信息)，不是这步转丢的。转换耗时约4.4小时CPU时间、峰值
     内存约40.8GB(29,635,224 SNP × 3024样本的EIGENSTRAT文本输出体量巨大，
     ~90GB量级，这不是卡住，是数据量真实需要这么久，接手人排队列/申请
     资源时可以参考这个数字)。
   - 转换后自动跑一次 `check_ref.py`，核对convertf本身没有把A1/A2顺序转错
   - `par.MERGE`：mergeit 合并 29M_3k(转换后) 与 6.7M_720，剔除
     strand-ambiguous(A/T、C/G)位点，取共享位点，产出写到
     `db/merged_29M3k_6M7_720/merged.{geno,snp,ind}`。
     **⚠️2026-08-08已在服务器执行完成(`mergeit -p par.MERGE`)，`##end of
     mergeit`正常输出，但`Histogram of checkmatch return codes`显示
     `total: 0`——两个panel之间0个SNP匹配上，产出的`merged.{geno,snp,ind}`
     预期是空的**。参数文件本身没问题(没有"unrecognized parameter"报错)，
     真正的问题是两边`.snp`文件SNP ID命名方式不同(29M_3k是纯数字物理位置
     如`1026`，720是`{chrom}np{pos}`如`1np1409`)，字符串永远对不上，很可能
     是mergeit按ID匹配而非按坐标匹配导致的——详见0节的完整诊断和建议的ID
     统一重命名修复方案。**这个诊断还没有100%坐实，下一个session先要用户
     贴mergeit完整log（不只是`##end`那几行）确认，再决定要不要按0节的方案
     重新生成ID统一的`.snp`副本重跑**——不要凭这段猜测直接改代码重复跑
     7000+秒的mergeit，先确认诊断对不对。
   - 需要`convertf`/`mergeit`在PATH里、`python3`能`import pysam`
     (`check_ref.py`用到)
3. **旱稻/水稻(或至少 indica/aus/japonica/aromatic 亚群)标签来源——
   找到强力候选(2026-08-07)，待核实覆盖率**：不再单纯指望SNP-Seek(官网
   前端下线)或翻论文补充材料——`asn720data/asn720.pop.fam`的FID列
   (`OrA-OrF`)很可能就是答案，见1.2节和3.2第1条。备用/交叉验证来源仍然
   保留：*A pangenome reference of wild and cultivated rice*(Nature 2025，
   PMID 40240605)里129份*O. rufipogon*的`Or-Ia/Or-Ib/Or-II/Or-IIIa/
   Or-IIIb`谱系分组(该文摘要提到"所有驯化位点都来自japonica祖先谱系
   Or-IIIa")，以及`db/wild_rice_pangenome_README.txt`(新发现，未读，
   **另一个会话窗口正在处理，等对方回复cat的内容，不要重复索取**)。
   三个来源不确定是否互相独立还是同源，第1条查清楚之后再确定最终用哪个/
   要不要交叉验证。这一步不阻塞2，可以并行。
4. pseudo-haplotype 调用脚本、smartpca 具体参数(尤其是 `-lsqproject` 相关
   配置)还没写，等1-3项确认完再动手，避免在错误坐标系/面板/标签上重复
   返工。
5. **核查`asn720.6m.ind`与`NB_final_snp.ind`之间的样本ID重叠**(见1.2节
   `CX382`案例)——去掉`_merged`后缀后系统性比对两份样本列表，确认重叠
   规模；如果存在系统性重复，第2节①步求交集前需要先定去重策略。当前
   优先级低于1/2/3，不阻塞mergeit先跑通，但正式产出结果前必须处理。

## 4. 与 besthit 分支的关系

本分析线的输入直接依赖 `codex/oryza-competitive-mapping` 分支产出的
`<sample>.besthit_oryza.fastq.gz`(见该分支 `docs/ORYZA_BESTHIT_HANDOFF.md`
第5.3节)。besthit 那边的分析(第7/8节:acc2taxid 覆盖度诊断、野生稻组装
物种身份诊断)结果，间接也会影响本分析线——如果 besthit 阶段发现需要扩充
competitive mapping 的参考库(把 `db/3k/wild/` 那140+野生稻组装加进去)，
"确认是Oryza"的read集合本身会变化，本分析线第2节①②步要用的输入也要
跟着重跑。两条线目前互相独立推进，但下游会汇合，需要留意上游变动。
