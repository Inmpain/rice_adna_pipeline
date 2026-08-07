# rice_adna_pipeline — 跨分支现状汇总
更新时间: 2026-08-07（本文档由 Claude Code 整理，来源为 GitHub 仓库全部4个分支的现有文档，未涉及服务器直接操作）

> 本仓库有4个分支，3个未合并进main。这份文档是四条线的横向汇总入口；
> 每条线的完整细节仍以各自分支的权威文档为准（下面逐条标注）。

---

## 全局依赖关系

```
codex/oryza-screen-merge (上游工具)
        │ 合并 shotgun+panel1+panel2 候选Oryza FASTQ
        ▼
codex/oryza-competitive-mapping (核心: besthit)
        │ 131数据库竞争比对 + 古DNA损伤校正过滤，产出干净Oryza reads
        ▼
codex/ecotype-pca-panel (下游: PCA生态型判定)
        │ 投影到现代群体PCA，判断旱稻/水稻

main (独立主线，不依赖上面三条)
   57开花基因命中分析、BWA vs Bowtie2决策、3K RGP交叉验证
```

---

## 分支①: `main` — 57基因命中主线

**权威文档**: `file_path.md`(完整路径地图) → `PROJECT_STATUS.md` → `docs/decisions_log.md`

**当前状态**: 已定论，暂无进行中的开放任务；`PROJECT_STATUS.md`本身已知曾被截断，只剩最后一次追加内容。

**关键决策**（已用9格矩阵测试验证）:
- 提取阶段(阶段①)必须用 **BWA**（不是Bowtie2）——决定性变量
- 定量比对阶段(阶段②) BWA 或 Bowtie2新参数(`-N1`) 均可接受
- 证据: `docs/09_extraction_mapping_matrix_final.md`

**57基因命中数据质量问题**（`docs/flank1kb_msa_exploration.md`）:
- 42条命中reads里(Bowtie2版数据)，35.7%落在低复杂度重复区域，可信度存疑
- 已建立可复用QC流程: dustmasker + BLAST交叉验证 + bedtools批量筛查
- 三个原定优先基因重新排序: **OsGI证据最干净(3/3)** > DTH8/Ghd8部分可信(3/5) > DTH7/OsPRR37已被证伪(0/1)
- ±1kb侧翼扩展**没有**捕捉到调控区数据，最高覆盖率仅2.985%，深度瓶颈是根本限制

**3K RGP交叉验证**（`docs/3krgp_integration_and_simulation_prep.md`）:
- 42条命中reads里6条精确匹配3K RG已知SNP（不含三个最初优先基因）
- DTH8/Ghd8发现~1116-1359bp大片段缺失簇，163/3024份材料(5.4%)携带，与文献吻合
- DTH8断点检验: 16个古稻样本在该坐标全部零覆盖（零信息，非阴性证据）
- 研究重心已从"逐基因深挖"转向"群体尺度"分析（PCA/f-statistics），对应 `ecotype-pca-panel` 分支

**待办**（NGSNGS模拟环境搭建相关）:
1. `-mf`损伤文件格式转换（mapDamage2输出 → NGSNGS格式）未写
2. `readlen_dist_bwa`/`mapdamage_out_bwa`需要用BWA数据源重新跑（sbatch并行中）
3. 16基因组染色体命名确认后正式跑NGSNGS
4. **passport分类表下载(SNP-Seek)** —— 这是PCA路径最大卡点，与`ecotype-pca-panel`分支的卡点是同一个问题
5. 163份DTH8缺失型材料需要亚群标签才能验证aus/tropical japonica富集

---

## 分支②: `codex/oryza-screen-merge` — 候选FASTQ合并工具

**权威文档**: `oryza_screen_merge/README.md`

**功能**: 极简Snakemake工作流，仅做一件事——把每个样本的 shotgun / capture panel1 / capture panel2 三份候选Oryza FASTQ合并成一份压缩文件。**不**跑BWA、不建/读BAM、不重比对IRGSP、不去重、不做MAPQ过滤、不call变异、不算基因命中。

**输入**（服务器路径，已在config.yaml写死）:
```
/home/scratch/yinmt202607/tests/param_matrix_bt2_vs_bwa/01.reads_combined/bwa/
├── <sample>.prefiltered.IRGSP1.mapped.fq
├── <sample>_RicePanel1.bwa.primary_mapped.fastq.gz
└── <sample>_RicePanel2.bwa.primary_mapped.fastq.gz
```

**输出**: `/home/scratch/yinmt202607/gene/results/oryza_candidates_combined/<sample>.oryza_candidates.combined.fastq.gz`

**部署方式**（注意：README明确写了不要碰主repo工作副本，单独下载到gene目录）:
```bash
mkdir -p /home/scratch/yinmt202607/gene/scripts/oryza_screen_merge
# curl 三个文件: Snakefile / config.yaml / README.md
# from codex/oryza-screen-merge 分支的 raw.githubusercontent.com URL
snakemake -n -p   # 先dry-run
snakemake --cores 1 --printshellcmds
```

**状态**: 这是下游 `oryza-competitive-mapping` 分支输入的来源，功能单一，无已知问题记录。

---

## 分支③: `codex/oryza-competitive-mapping` — Best-hit 过滤（核心工作）

**权威文档**: `docs/ORYZA_BESTHIT_HANDOFF.md`（36KB，2026-08-05更新，最详细）

**脚本**: `scripts/oryza_besthit/`
- `oryza_besthit_damage_filter.py`（~550行）— 核心算法
- `submit_oryza_besthit.sh`（~430行）— SLURM提交脚本，5种模式(check/smoke/local/submit/merge)
- `submit_oryza_competitive_mapping.sh` — 上游mapping提交脚本（⚠️服务器上有v2/v3等多个版本号文件，未确认哪个和git同步，接手时需要先diff）

**核心算法**: 每条alignment算 `adjusted_NM = NM - terminal_damage_count`（只在读长两端`--damage-window`默认5bp内、且方向正确的C→T/G→A替换才被扣除）。按species分组取每个species最优alignment，Oryza最优hit的adjusted_NM ≤ 非Oryza最优hit → KEEP，否则REJECT。

**参数含义**:
| 参数 | 默认值 | 含义 |
|---|---|---|
| `--damage-window` | 5 (bp) | 读长两端多少bp内的替换有资格被当作末端损伤扣除 |
| `--top-n` | 10 | 非Oryza物种审计表里每条read记录的候选数（Oryza命中不占名额，必录） |
| `--oryza-taxids` | "4529 4530 4536" | rufipogon/sativa/nivara三个taxid |
| `JOB_MEM_MB` | 16000 | SLURM内存申请，⚠️未用`seff`实测校准过 |

**已修复的3个bug**（服务器实跑中发现）:
1. `module load python/ 2>/dev/null` 只重定向了stderr，集群把错误写到stdout → 改成 `>/dev/null 2>&1`
2. **最关键**: sbatch把脚本拷到spool目录执行，`BASH_SOURCE`推导路径会指向错误目录 → 改用`--export="ALL,PY_SCRIPT=..."`显式传绝对路径
3. 缺少批量模式 → 新增 `submit all` / `local all`

**当前进度**（截至2026-08-05）:
- mapping阶段: 3个样本完成(LV6000619499/619917/620016)，其余13个仍在SLURM队列
- besthit阶段: 仅1个样本做过1000-read smoke test

**Smoke test结果**（LV6000619499, 前1000条reads）:
```
input=1000, with_alignment=1000, with_oryza_hit=956, kept=103,
rejected_nonoryza_better=853, rejected_no_oryza=44, unclassified=0
```

**当前最优先的开放问题（第7.5节）**: `--damage-window=5bp`可能偏窄
- 最初假设"数据库野生稻缺失"已被证伪：Oryza属18个种普遍有几十到几百条contig，*sativa*(111条)反而覆盖度偏低，*rufipogon*(717条)更好
- REJECT名单里赢家物种(桉树/蜗牛/深海海绵/骡鹿等)生物学上无一致性，更像随机撞库
- **新证据**: 差1点惜败REJECT的541条read平均长度(~68.5bp)明显比守住KEEP的103条(~62.3bp)更长 → 支持"固定5bp窗口对长读长不利"的假说
- **下一步待验证**: 这批angkor样本有没有做过标准mapDamage风格损伤衰减曲线分析？衰减到背景水平大概多少bp？没有的话需要从BAM直接估算

**待办优先级顺序**（HANDOFF文档第9节）:
1. 【最优先】损伤衰减曲线诊断（第7.5节）
2. 如证实窗口偏窄，调大`--damage-window`重跑smoke test
3. 【降级，非阻塞】`db/3k/wild/`140+野生稻组装物种身份、taxid rank校验、内存实测
4. 只有排除损伤窗口问题后才评估加margin参数
5. 剩余13个样本随mapping进度陆续跑besthit（不依赖1-4，随时可做）
6. 参考基因组体系整理（优先级已降低）
7. 8.2选择信号扫描、8.3的57基因SV判生态型——都还没开始

---

## 分支④: `codex/ecotype-pca-panel` — PCA生态型判定（下游）

**权威文档**: `docs/ECOTYPE_PCA_PANEL.md`（2026-08-05首次写下）

**依赖**: 直接依赖`oryza-competitive-mapping`分支产出的`<sample>.besthit_oryza.fastq.gz`

**脚本**: `scripts/ecotype_pca/check_ref.py`（坐标系核对用）

**两个参考panel**:
| Panel | 路径 | 格式 | 内容 |
|---|---|---|---|
| 29M_3k | `db/29M_3k/` | PLINK bed/bim/fam | 3024份3K RG驯化稻材料，~2900万biallelic SNP，vs Nipponbare MSU7坐标系 |
| 6.7M_720 | `db/6.7M_720/` | EIGENSTRAT(.geno/.ind/.snp) | 720份样本(以野生稻为主)，~670万SNP位点 |

**已确认（坐标系核对，用check_ref.py完成）**:
- 染色体命名一致（都是裸数字1-12，不是chr01风格）——之前手动改的chr01版本(`.bim.chrfix`)是多余的，不要用
- REF/ALT方向相反：29M_3k是A2=REF/A1=ALT(干净)；6.7M_720是A1=REF/A2=ALT占91.5%(183/200)，但方向相反不阻塞mergeit(工具自动纠正)，只有A/T、C/G这类strand-ambiguous位点需要手动剔除
- `asn720data`(94,974个SNP，main分支`file_path.md`提到的旧数据)**已决定弃用**，密度太低，不再需要确认它和`asn720.6m.*`的关系

**分析设计（草案，5步）**:
1. `29M_3k`(convertf转EIGENSTRAT) ∩ `6.7M_720`(mergeit求交集，剔除strand-ambiguous位点)
2. 每个古代样本单独: 上一步结果 ∩ 该样本besthit后reads实际覆盖到的位点
3. pseudo-haplotype调用（每位点随机抽一条覆盖read取碱基，标准古DNA做法，不做常规diploid calling）
4. smartpca：现代样本建PCA空间，古代样本用`-lsqproject`投影模式
5. 看古代样本落在哪个现代亚群附近判断生态型

**当前状态**: 坐标系核对完成，convertf+mergeit方案已确定但**还没跑**，是当前最优先要执行的步骤。pseudo-haplotype脚本和smartpca具体参数都还没写。

**最大卡点**: 旱稻/水稻(或至少indica/aus/japonica/aromatic亚群)标签来源仍未解决——SNP-Seek官网前端下线，需要找替代来源(如GigaDB上3K RGP论文补充材料)。这不阻塞前两步，可以并行找。

**与besthit分支的联动**: 如果besthit那边（第7/8节）发现需要扩充competitive mapping参考库（比如把`db/3k/wild/`140+野生稻组装加进去），"确认是Oryza"的read集合会变化，这条线的第1/2步要跟着重跑。目前两条线独立推进，下游会汇合。

---

## 已知的软件/参数速查表（跨分支汇总）

| 软件 | 用于 | 关键参数 |
|---|---|---|
| BWA aln | 提取阶段(阶段①)、定量比对(阶段②) | 主线选择，决定性因素在提取阶段 |
| Bowtie2 | 历史/对比用 | `-N1`新参数可接受，旧参数最差 |
| bowtie2 (competitive mapping) | besthit上游 | `-k 100 -L 22 -i S,1,1.15 --mp 1,1 --rdg 0,1 --rfg 0,1 --score-min L,0,-0.1 --no-unal` |
| dustmasker | 低复杂度区域标记 | 全基因组扫描 |
| bedtools intersect | 低复杂度过滤/SNP交叉 | 用整条read跨度而非单点坐标(v2版，35.7%，比v1版14.3%更准) |
| mapDamage2 | 古DNA损伤统计 | 需sbatch并行(单线程，MCMC贝叶斯估计慢) |
| NGSNGS | 模拟古DNA reads | v0.9.2.2；`-lf`长度CDF、`-mf`损伤频率(格式转换未完成)、`-seq SE` |
| samtools consensus | flank1kb序列提取 | 零覆盖不补N，只输出有覆盖片段(踩坑点) |
| pysam | besthit脚本核心依赖 | `get_aligned_pairs(with_seq=True)`拿逐位点比对 |
| convertf/mergeit (EIGENSOFT) | PCA panel格式转换/求交集 | 用于ecotype-pca-panel分支 |
| smartpca | PCA投影 | `-lsqproject`模式，还未跑 |

---

## 本次整理说明

本文档由 Claude Code 在会话中读取 GitHub 仓库全部4个分支后整理，未接触服务器，
所有服务器路径/参数信息均转述自仓库内已有文档，未做实测验证。
按 `github-repo-protocol` skill 的Rule 2要求生成，作为跨分支的横向索引，
不替代各分支内的权威文档（已在各节标注）。

后续维护：这份文档只放在 main 分支（不在其余3个分支各放一份），
因为它本身就是跨分支的横向索引，内容天然应该只有一个权威副本，
避免像 `PROJECT_STATUS.md` 那样出现多处不同步的历史。
