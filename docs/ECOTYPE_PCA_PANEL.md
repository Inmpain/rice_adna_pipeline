# 旱稻/水稻生态型 PCA 分析——数据库与设计草案

> 本文档记录一条与 `codex/oryza-competitive-mapping`(besthit)平行、下游依赖
> 其输出的新分析线：把 besthit 过滤后确认为 Oryza 的古代 reads，投影到现代
> 群体 PCA 空间，判断每个古代样本更接近旱稻还是水稻生态型。
> 首次写下: 2026-08-05

---

## 0. 现状一句话总结

两个参考 panel 已经到位(见第1节实际路径)，但**格式不统一**(3K是PLINK压缩
格式，720是EIGENSTRAT格式)，且坐标系一致性、与`asn720data`的关系都还没
核实。当前最优先的事是把第3节列的开放问题确认掉、把两个panel转成同一格式
求交集，再动手写pseudo-haplotype调用脚本。

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
- **与 `asn720data` 的关系——文件名前缀`asn720`基本实锤了**：`file_path.md`
  第二节记录的 `asn720data/asn720.pop.{bed,bim,fam}`(720份现代/近现代品种
  PLINK面板，94,974个SNP位点)，跟这里的`asn720.6m.*`共享同一个`asn720`
  命名前缀，几乎可以确定是**同一批720份材料，这次用更高密度(6.7M vs
  94,974，约70倍)重新call了SNP**。仍需核实一遍(比如对比`.ind`文件里的
  样本ID列表跟`asn720.pop.fam`是否完全一致)，但不用再当成两个独立数据源
  重新调研关系。

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
   两者格式不同，需要先统一格式再求交集：
   a) 用 EIGENSOFT `convertf` 把 `29M_3k`(PLINK bed/bim/fam)转成EIGENSTRAT
      (跟 6.7M_720 已经是的格式统一，也是smartpca的原生输入格式)
   b) 用 EIGENSOFT `mergeit` 在两个EIGENSTRAT数据集之间按染色体+物理位置
      (以及等位基因方向)求交集，产出合并后的panel
   (前提：两边坐标系必须一致，见1.2的待确认项——mergeit本身不会告诉你
   "坐标系不一致"，只会给出一个异常小的交集，容易被误当成"正常结果")
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

1. **【最优先】确认两个panel坐标系一致**：
   ```bash
   zcat /home/scratch/yinmt202607/db/29M_3k/NB_final_snp.bim.gz | head
   head /home/scratch/yinmt202607/db/6.7M_720/asn720.6m.snp
   ```
   看染色体命名(是否都是`chr01`风格、还是一个`chr01`一个`1`)、位置数值量级
   是否都落在IRGSP1.0/MSU7对应染色体长度范围内。不确认这点，`mergeit`求
   交集不会报错，只会默默给出一个异常小的交集，容易被当成"正常结果"接受。
2. **确认 `6.7M_720` 和 `asn720data` 的关系**：
   ```bash
   head /home/scratch/yinmt202607/db/6.7M_720/asn720.6m.ind
   head /home/scratch/yinmt202607/asn720data/asn720.pop.fam
   ```
   对比样本ID列表是否一致，确认是否真的是同一批720份材料的更高密度版本
   (见1.2节，文件名前缀已经很提示，但没有逐条核实过)。如果确认一致，
   `file_path.md`里"与16个angkor样本的关系尚未最终确认"这条注记直接继承
   过来即可，不用重新走一遍确认流程。
3. **格式统一 + 求交集**(1-2确认无误后再做，顺序不能反)：
   ```bash
   # a) 29M_3k 解压 + convertf 转 EIGENSTRAT
   cd /home/scratch/yinmt202607/db/29M_3k
   gunzip -k NB_final_snp.bed.gz NB_final_snp.bim.gz NB_final_snp.fam.gz
   # 写 convertf 的 par 文件(genotypename/snpname/indivname 指向PLINK三件套,
   # outputformat EIGENSTRAT), 然后:
   convertf -p par.PLINK.EIGENSTRAT

   # b) mergeit 合并 29M_3k(转换后) 与 6.7M_720，取共享位点
   mergeit -p par.MERGE
   ```
4. **旱稻/水稻(或至少 indica/aus/japonica/aromatic 亚群)标签来源仍未解决**——
   这是 `docs/3krgp_integration_and_simulation_prep.md` 待办第4条就已经标注
   的老问题("passport分类表下载(SNP-Seek)—— PCA路径最大卡点")，SNP-Seek
   官网前端目前下线，需要找替代来源(比如 GigaDB 上的 3K RGP 论文补充材料)。
   没有这份标签，PCA做出来也没法解释哪个方向是旱稻哪个是水稻。这一步不
   阻塞1-3，可以并行找。
5. pseudo-haplotype 调用脚本、smartpca 具体参数(尤其是 `-lsqproject` 相关
   配置)还没写，等1-3项确认完再动手，避免在错误坐标系/面板上重复返工。

## 4. 与 besthit 分支的关系

本分析线的输入直接依赖 `codex/oryza-competitive-mapping` 分支产出的
`<sample>.besthit_oryza.fastq.gz`(见该分支 `docs/ORYZA_BESTHIT_HANDOFF.md`
第5.3节)。besthit 那边的分析(第7/8节:acc2taxid 覆盖度诊断、野生稻组装
物种身份诊断)结果，间接也会影响本分析线——如果 besthit 阶段发现需要扩充
competitive mapping 的参考库(把 `db/3k/wild/` 那140+野生稻组装加进去)，
"确认是Oryza"的read集合本身会变化，本分析线第2节①②步要用的输入也要
跟着重跑。两条线目前互相独立推进，但下游会汇合，需要留意上游变动。
