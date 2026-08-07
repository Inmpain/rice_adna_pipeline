# 旱稻/水稻生态型 PCA 分析——数据库与设计草案

> 本文档记录一条与 `codex/oryza-competitive-mapping`(besthit)平行、下游依赖
> 其输出的新分析线：把 besthit 过滤后确认为 Oryza 的古代 reads，投影到现代
> 群体 PCA 空间，判断每个古代样本更接近旱稻还是水稻生态型。
> 首次写下: 2026-08-05

---

## 0. 现状一句话总结

两个参考 panel 已经到位(见第1节实际路径)。坐标系核对已完成(见3.1)：染色体
命名一致、REF/ALT方向虽相反但mergeit可自动处理；`asn720data`已决定弃用
(密度太低)。convertf 这步(29M_3k PLINK→EIGENSTRAT)已在服务器上验证跑通
(2026-08-07，见3.2节step 2)：29,635,224个SNP、3024个样本转换前后完全一致，
无丢失。当前最优先的事是接着跑 `run_convert_merge.sh` 剩下的 mergeit 那步
(29M_3k转换结果 与 6.7M_720 求交集)，再动手写pseudo-haplotype调用脚本。

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
- **`asn720data`(94,974个SNP)已决定弃用，不再需要确认两者关系**：
  `file_path.md`记录的 `asn720data/asn720.pop.{bed,bim,fam}`(720份现代/
  近现代品种PLINK面板，94,974个SNP位点)密度太低，撑不起本分析"panel密度
  换古代reads覆盖概率"的诉求，直接不用，只用`asn720.6m.*`这份6.7M位点的
  数据。两者是否为同一批720份材料这个问题不再重要，不用再花时间核实。
- **样本ID观察**：`asn720.6m.ind`里样本ID是两种风格混合——`ERR068594`一类
  (ENA测序run编号，群体标签`control`)和`B011_merged`一类(跟3K RG自己的
  品种编号风格接近，见besthit分支`ORYZA_BESTHIT_HANDOFF.md`第6节SV数据里
  出现过的`B071`类样本名)。说明这720份很可能是从不同批次/不同命名体系
  拼合而成，不是单一来源；暂不影响使用，仅作记录。

### 1.3 两个panel的关系

分析设计里两者是**互补而非替代**关系：3K(29mio) 覆盖驯化稻的谱系多样性，
6.7M_720 补充野生稻/近缘种一侧——古稻样本理论上落在"驯化-野生"谱系的某个
位置，两个panel合起来才能给出有意义的PCA背景。

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
   (标签来源见第3节，目前还没有)判断生态型归属
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

1. ~~确认 `6.7M_720` 和 `asn720data` 的关系~~——**已决定弃用`asn720data`
   (94,974个SNP密度太低)，只用`asn720.6m.*`，这条不再需要，见1.2节。
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
     本来就没有群体信息)，不是这步转丢的，对应第3条待办。
   - 转换后自动跑一次 `check_ref.py`，核对convertf本身没有把A1/A2顺序转错
   - `par.MERGE`：mergeit 合并 29M_3k(转换后) 与 6.7M_720，剔除
     strand-ambiguous(A/T、C/G)位点，取共享位点，产出写到
     `db/merged_29M3k_6M7_720/merged.{geno,snp,ind}`。
     **⏳尚未在服务器验证**——`par.MERGE`里的参数名(`geno1/snp1/ind1`、
     `genooutfilename`等)是按EIGENSOFT mergeit的标准约定写的，鉴于
     convertf那步"标准约定"就踩过一次坑(文件后缀识别)，这步大概率也会
     报错，报错的话把完整输出贴回来，跟convertf那次一样查server上
     mergeit的README/源码定位，不要凭README字面直接猜下一次修复
   - 需要`convertf`/`mergeit`在PATH里、`python3`能`import pysam`
     (`check_ref.py`用到)
3. **旱稻/水稻(或至少 indica/aus/japonica/aromatic 亚群)标签来源仍未解决**——
   这是 `docs/3krgp_integration_and_simulation_prep.md` 待办第4条就已经标注
   的老问题("passport分类表下载(SNP-Seek)—— PCA路径最大卡点")，SNP-Seek
   官网前端目前下线，需要找替代来源(比如 GigaDB 上的 3K RGP 论文补充材料)。
   没有这份标签，PCA做出来也没法解释哪个方向是旱稻哪个是水稻。这一步不
   阻塞1-2，可以并行找。
4. pseudo-haplotype 调用脚本、smartpca 具体参数(尤其是 `-lsqproject` 相关
   配置)还没写，等1-2项确认完再动手，避免在错误坐标系/面板上重复返工。

## 4. 与 besthit 分支的关系

本分析线的输入直接依赖 `codex/oryza-competitive-mapping` 分支产出的
`<sample>.besthit_oryza.fastq.gz`(见该分支 `docs/ORYZA_BESTHIT_HANDOFF.md`
第5.3节)。besthit 那边的分析(第7/8节:acc2taxid 覆盖度诊断、野生稻组装
物种身份诊断)结果，间接也会影响本分析线——如果 besthit 阶段发现需要扩充
competitive mapping 的参考库(把 `db/3k/wild/` 那140+野生稻组装加进去)，
"确认是Oryza"的read集合本身会变化，本分析线第2节①②步要用的输入也要
跟着重跑。两条线目前互相独立推进，但下游会汇合，需要留意上游变动。
