# Oryza competitive mapping / best-hit 接手说明

更新时间：2026-08-05（本次更新：把 best-hit 阶段从设计到服务器实跑、debug、
结果分析、以及下游三条工作线的讨论全部并入，供交接给另一个 AI / 协作者使用）

> 读这份文档前提：你没有服务器直接执行权限，所有服务器端命令的输出都需要
> 人类协作者手动跑了贴回来。不要假设你能直接看到服务器文件系统。

## 0. 现状一句话总结

Best-hit 脚本（`oryza_besthit_damage_filter.py` + `submit_oryza_besthit.sh`）
已经在服务器上实际跑通（3个样本里的1个做过1000-read smoke test，其余2个mapping
已完成但besthit还没跑全量）。脚本本身没有已知 bug 了（过程中发现并修复了3个）。
**当前卡点不是代码，是方法论**：best-hit 的判定完全依赖 acc2taxid 数据库里
Oryza 属参考基因组的覆盖面，而现在这个覆盖面大概率只有 *O. sativa*（野生稻
*rufipogon*/*nivara* 有没有、有多少条 contig，这个诊断还没做完，见第7节）。
smoke test 的实际数字（第6节）显示"输赢边界"高度集中在差1个编辑距离，这个
模式本身就是"参考单一化"的间接证据。**在把这个诊断做完、并且视情况把
`db/3k/wild/` 下140+个野生稻/近缘种组装并入数据库之前，不建议去调
KEEP/REJECT 的判定阈值——那是在一个已知有偏的比较基础上调参，治标不治本。**

## 1. 项目背景

古 DNA / 环境 DNA 水稻项目（angkor，16个考古样本，年代跨度约公元1160-1981年，
元数据见服务器 `/home/scratch/yinmt202607/angkor_robot_library.txt`）。候选水稻
reads 已经由 shotgun、capture panel1、panel2 三份 FASTQ 合并完成。当前正在把
每个样本分别比对到 WGS 真核数据库、亚洲水稻 panel 和 IRGSP（competitive
mapping），然后执行 best-hit，排除更像非水稻物种的 reads，产出"确认是 Oryza"
的干净 FASTQ 供下游使用。

## 2. 关键路径

### GitHub

仓库：https://github.com/Inmpain/rice_adna_pipeline

当前分支：`codex/oryza-competitive-mapping`

本机（Mac）GitHub 工作副本：`/Users/inmpain/Documents/angkor/rice_adna_pipeline_publish`

本机参考代码（**不在这个 git 仓库里，是另外两个独立项目，只用来抄逻辑，不是
本项目依赖**）：
- 原 best-hit 项目（ngsLCA 风格通用分类器，最终**没有**直接复用，见第5节）：
  `/Users/inmpain/github/aeDNA_popgen`
- WGS mapping/动态内存参考：`/Users/inmpain/github/Pipeline_snakemake/new_single_multi`

### 服务器（`/home/scratch/yinmt202607/gene/scripts/`）

这个目录下当前实际存在的文件（用户在服务器上 `ls` 看到的，注意有历史版本
文件，不要跑错）：

```
merge_oryza_fastq.sh
oryza_besthit_damage_filter.py          # <- besthit 主脚本，本次工作的产出
oryza_screen_merge/
submit_oryza_besthit.sh                 # <- besthit 提交脚本，本次工作的产出
submit_oryza_competitive_mapping_old.sh
submit_oryza_competitive_mapping.sh     # <- 当前在用的 mapping 提交脚本
submit_oryza_competitive_mapping_v2.sh
submit_oryza_competitive_mapping_v3.sh
test_submit_oryza_competitive_mapping.sh
test_submit_oryza_competitive_mapping_v2.sh
test_submit_oryza_competitive_mapping_v3.sh
```

**未确认**：`submit_oryza_competitive_mapping.sh`（无版本号后缀）是否就是
git 仓库 `scripts/oryza_besthit/submit_oryza_competitive_mapping.sh` 的当前
版本，还是本地又手改过、和 git 不同步——git仓库里只有一份（无 v2/v3），服务器
上有多个版本号文件，说明 mapping 阶段的脚本在服务器上有过 git 之外的迭代。
接手时先 `diff` 一下确认，不要默认服务器和 git 一致。

服务器上下载 git 脚本的方式（这台机器 login node 能连外网 github.com，已验证
可用）：
```bash
cd /home/scratch/yinmt202607/gene/scripts
curl -fsSL -o <文件名> \
  https://raw.githubusercontent.com/Inmpain/rice_adna_pipeline/codex/oryza-competitive-mapping/scripts/oryza_besthit/<文件名>
```

服务器 python 环境：conda/mamba 环境 `/home/usr/yinmt/.local/mamba/snakemake`
（提示符显示 `(base)` 或 `(.../snakemake)`），`python3` 直接在 PATH 上，已确认
装了 `pysam`。**没有 `python/` modulefile**（`module load python/` 会报错，
脚本里已经用 `|| true` 吞掉了，见第4.1节）。

候选 FASTQ：`/home/scratch/yinmt202607/gene/results/oryza_candidates_combined/<sample>.oryza_candidates.combined.fastq.gz`

mapping 输出根目录：`/home/scratch/yinmt202607/gene/results/oryza_competitive_mapping/`
- 各数据库 BAM：`bam_by_database/<sample>/<sample>.<database>.bam`
- 每样本合并 BAM（besthit 的输入）：`by_sample/<sample>.competitive.name_sorted.bam`（+ `.finished` 标记）
- mapping 日志/提交记录/串行状态：`logs/` `submissions/` `series/`

besthit 输出根目录：`/home/scratch/yinmt202607/gene/results/oryza_competitive_mapping/besthit/`
（脚本会自动建，见第5节的输出说明）；smoke test 单独输出到 `besthit/smoke_test/`。

## 3. 数据库

WGS 真核数据库：`/home/database/ref20250728/cph_euk/wgs_eukaryota.{1..129}.fas.gz`（129个shard）

亚洲水稻 panel：`/home/scratch/yinmt202607/db/asian_rice_panel_index/asian_rice_panel.fa`

IRGSP：`/home/scratch/yinmt202607/db/asian_rice_panel_index/irgsp_bt2idx`

合并 accession→taxid：`/home/scratch/yinmt202607/db/asian_rice_panel_index/all_wgs_asian_irgsp.acc2taxid`
（标准 NCBI 三列格式：accession / accession.version / taxid）

NCBI taxonomy dump：`/home/database/ref20250728/taxonomy_CPH/ncbi/20250530/{nodes.dmp,names.dmp}`

目标 Oryza taxid（species rank，**已用 --oryza-taxids 默认值写死在脚本里，
但 rank 是不是真的 species 没有校验过，见第7节待办1**）：
- 4529：*Oryza rufipogon*（野生稻/普通野生稻）
- 4530：*Oryza sativa*（亚洲栽培稻）
- 4536：*Oryza nivara*（尼瓦拉野生稻）

### ⚠️ 潜在的额外野生稻参考资源（还没用上，见第7/8节）

`file_path.md`（`main` 分支）提到服务器上 `/home/scratch/yinmt202607/db/3k/wild/`
下有 **140+ 个 `{SampleID}.transfer.merge.chr.fasta`（野生稻/近缘种染色体级组装
基因组）**，来自 3K Rice Genome Project 相关资源，**目前完全不在
`all_wgs_asian_irgsp.acc2taxid` / competitive mapping 的131个数据库体系里**。
这些样本具体是什么物种、有没有已知的 taxid，仓库文档里没写，需要上服务器查
（命令见第7节）。如果这批组装真的能用，是目前看来解决"Oryza 参考单一化"
问题最直接的路径。

## 4. Mapping 阶段（competitive mapping，besthit 的上游）

每样本比对131个数据库（129 WGS shard + asian_rice_panel + IRGSP），Bowtie2：
```bash
-k 100 -L 22 -i S,1,1.15 --mp 1,1 --rdg 0,1 --rfg 0,1 --score-min L,0,-0.1 --no-unal
```

截至 2026-08-05，`by_sample/` 下已完成（`.finished` 标记齐全）的样本：
- LV6000619499（唯一做过 besthit smoke test 的样本，见第6节）
- LV6000619917
- LV6000620016

其余样本仍在串行队列（`.../series/`）里陆续跑。**besthit 阶段不需要等全部
mapping 完成**——`submit_oryza_besthit.sh submit all` / `local all` 会自动只
处理当前已有 `.finished` mapping BAM 的样本，可以随 mapping 进度反复重跑
（已完成 besthit 的样本会被跳过，见第5.4节）。

## 5. Best-hit / aDNA 损伤校正过滤（本次工作的主体）

### 5.1 为什么没有直接复用 aeDNA_popgen 的通用分类器

`aeDNA_popgen/01_TaxAssign/workflow/scripts/besthit_competitive.py` 是一个
通用的 ngsLCA 风格层级分类器：沿谱系从叶子往根走，每层比较"子树内最优 NM"
vs"子树外最优 NM"，赢的最深节点即为归属。这套逻辑**不区分古 DNA 末端损伤和
真实突变**，也不是"Oryza vs 非 Oryza 二元竞争"这个具体需求。用户明确要求
（原话）：

> 参考代码不能直接照搬，因为我们需要增加古DNA末端损伤校正，并且只判断
> Oryza 与非 Oryza 的竞争关系。

所以这里是专门重写的一套逻辑，核心是：**扣除末端损伤后的 `adjusted_NM`，
Oryza 最优 hit 必须与非 Oryza 最优 hit 打平或更好，才 KEEP**。

脚本落地：
- `scripts/oryza_besthit/oryza_besthit_damage_filter.py` —— 主统计脚本（约550行）
- `scripts/oryza_besthit/submit_oryza_besthit.sh` —— SLURM 提交脚本（约430行）

### 5.2 核心算法

每条 alignment 计算：

| 字段 | 含义 |
|---|---|
| `NM` | bowtie2 的 NM tag（原始编辑距离，含替换+indel） |
| `substitution_count` | 从 MD/CIGAR 精确统计的真实替换数（不含 indel），靠 `pysam.get_aligned_pairs(with_seq=True)` 拿到逐位点 (ref_base, read_base)，lowercase ref_base = 该位点是 mismatch |
| `terminal_damage_count` | 只在读长两端 `--damage-window`（默认5）范围内、且**换算回 read 自身 5'→3' 方向后**符合古 DNA 损伤特征（5'端 C→T，3'端 G→A）的替换才计入 |
| `adjusted_NM` | `NM - terminal_damage_count`（indel 永远不扣） |

**reverse-strand 的方向换算是全脚本最容易写错的地方**：BAM 里 SEQ 字段已经是
原始 read 的反向互补（SAM 规范），要拿到"read 自身 5'→3' 方向"，需要把
reference-forward 顺序的 (ref_base, read_base) 序列**倒序 + 互补**两步一起做——
只倒序不互补，会把参考链上的 G→A 信号错当成"3'端 G→A"留在错误的位置，
从而漏判或错判。已经用手工构造的 reverse-strand 测试 BAM 验证过这一步
（见5.3）。

按 species taxid 分组（BAM reference name → acc2taxid → 沿 nodes.dmp 上溯到
最近的 `rank=="species"` 祖先），每个 species 只保留最优 alignment，排序
优先级：`adjusted_NM` 最小 → `NM` 最小 → `substitution_count` 最小 → `AS`
最大 → reference 名称（稳定/可重复）。非 Oryza species 取前 `--top-n`
（默认10）写入审计表；Oryza（`--oryza-taxids`，默认4529/4530/4536）**不占
非 Oryza 的10个名额**，命中就必录。

判定（`best_nonoryza` / `best_oryza` 各自是排名第一的非Oryza/Oryza species）：

1. 没有 Oryza hit → **REJECT** (`no_oryza_hit`)
2. 有 Oryza hit、没有非 Oryza hit → **KEEP** (`oryza_only_no_competitor`)
3. `best_oryza.adjusted_NM <= best_nonoryza.adjusted_NM` → **KEEP** (`oryza_at_least_as_good`)
4. 否则 → **REJECT** (`nonoryza_better`)

不用 MAPQ（bowtie2 `-k 100` competitive mapping 下不可靠，`asian_rice_panel/
IRGSP` 那种多参考竞争场景下 MAPQ 基本是固定值，没有信息量）。BAM 没有去重、
也**不需要**在这一步去重——同一 read 的多行来自131个数据库 × `-k 100`，不是
PCR duplicate，去重逻辑放在更下游。

### 5.3 输出

每样本，全部 `.tmp` 写完原子 `os.rename`；**只有一致性检查通过才写
`<sample>.finished`**：

| 文件 | 内容 |
|---|---|
| `<sample>.besthit.top10_species.tsv.gz` | 审计表，一条 read 对应多行（非Oryza top10 + 必录的Oryza行），列：`read_name read_length species_hit_count top10_rank is_always_included species_taxid species_name reference_name NM substitution_count terminal_damage_count adjusted_NM AS` |
| `<sample>.oryza_filter.decisions.tsv.gz` | 每 read 一行的最终判定，列：`read_name best_nonoryza_taxid best_nonoryza_name best_nonoryza_NM best_nonoryza_damage best_nonoryza_adjusted_NM best_oryza_taxid best_oryza_name best_oryza_NM best_oryza_damage best_oryza_adjusted_NM decision reason` |
| `<sample>.besthit_oryza.fastq.gz` | 从候选 FASTQ 按 KEEP read name 抽出，每条只输出一次，序列/质量值原样保留 |
| `<sample>.summary.tsv` | 单样本一行：`sample input_reads reads_with_alignment reads_with_oryza_hit kept_reads rejected_nonoryza_better rejected_no_oryza unclassified_reads`，脚本内部自查 `input_reads == kept + rejected_nonoryza_better + rejected_no_oryza + unclassified_reads`，不一致会非零退出（但这个文件本身仍会被写出，方便debug；`.finished` 不会被写） |

**`unclassified_reads` 的定义要注意**：包含两类，(a) 候选 FASTQ 里的 read 但
BAM 里完全没有任何 alignment（bowtie2 `--no-unal` 直接没输出，这类 read
**不会**出现在 `decisions.tsv` 里，只能从 `summary.tsv` 的
`input_reads - reads_with_alignment` 间接算出数量，查不到具体是哪些
read）；(b) BAM 里有 alignment，但每一条都解析不到 taxid/species（会在
`decisions.tsv` 里留一行 `decision=UNCLASSIFIED, reason=no_resolvable_species`，
可以按名字查到）。

多样本汇总 `besthit_summary.tsv` **不**由 Python 脚本直接写（避免并行 SLURM
作业竞争同一个文件），而是 `submit_oryza_besthit.sh merge` 事后拼接每个样本
自己的 `<sample>.summary.tsv`。

### 5.4 用法（`submit_oryza_besthit.sh` 的5种模式）

```bash
cd /home/scratch/yinmt202607/gene/scripts

# 1. 校验路径 / python3+pysam / SLURM partition，不提交任何作业
bash submit_oryza_besthit.sh check

# 2a. SLURM 队列：单个样本前1000条 reads 做 smoke test
#     （独立输出目录 besthit/smoke_test/，不写 .finished，不影响正式产出）
bash submit_oryza_besthit.sh smoke LV6000619499

# 2b. 前台本地跑（不进 SLURM 队列，适合调试/小样本快速出结果；量产还是用3a/3b）
bash submit_oryza_besthit.sh local LV6000619499
bash submit_oryza_besthit.sh local all      # 全部已完成mapping的样本，依次前台跑

# 3a. SLURM 队列：全量跑（自动跳过 besthit/<sample>.finished 已存在的样本）
bash submit_oryza_besthit.sh submit LV6000619499 LV6000619917 LV6000620016

# 3b. 等价写法：自动发现所有已完成 mapping 的样本（可随 mapping 进度反复重跑）
bash submit_oryza_besthit.sh submit all

# 4. 所有作业跑完后，拼接每个样本的 summary.tsv
bash submit_oryza_besthit.sh merge
# -> .../oryza_competitive_mapping/besthit/besthit_summary.tsv
```

可覆盖的环境变量（默认值见脚本头部）：`BAM_DIR` `FASTQ_DIR` `ACC2TAXID`
`NODES` `NAMES` `ORYZA_TAXIDS`（默认 `"4529 4530 4536"`）`DAMAGE_WINDOW`
（默认5）`TOP_N`（默认10）`OUT_DIR` `SLURM_PARTITION`（默认`comp`）
`JOB_CPUS`（默认4）`JOB_MEM_MB`（默认16000，**未实测校准，见第7节待办4**）
`JOB_TIME`（默认`04:00:00`）。

若 `check` 报 `import pysam` 失败，需要先在这个 `python3` 所在环境装 pysam。

### 5.5 已完成的本地验证（本机 Mac，无服务器数据/无sbatch时做的）

用 `samtools view -bS` + `sort -n` 手工构造了一个最小 BAM（5条参考、5条候选
read，覆盖 KEEP-打平/KEEP-oryza独占/REJECT-惜败/UNCLASSIFIED/BAM里完全没有
alignment 这5种情况），在隔离 venv 里装 pysam 实测跑通全流程，重点验证：

- **reverse-strand 末端损伤校正确实会反转最终判定**：构造了一条 read，原始
  NM 下"非 Oryza 更优"应判 REJECT，扣除末端损伤后 `adjusted_NM` 打平，正确
  判成 KEEP（`oryza_at_least_as_good`）——这是本脚本存在的核心原因，是最关键
  的一条验证。
- `input_reads == kept + rejected_nonoryza_better + rejected_no_oryza +
  unclassified_reads` 自检通过；FASTQ 抽取数等于 `kept_reads`。
- `--limit-reads`（smoke test）路径确认不写 `.finished`。
- 重复跑同一个样本（模拟 SLURM 任务失败重提交）不会因为 `.tmp` 残留或
  rename 冲突出错（幂等）。
- `submit_oryza_besthit.sh run`（sbatch job body）、`local`、`merge` 三种
  调用路径都跑通。

### 5.6 Bug 修复历史（服务器实跑中发现，均已修复并推送）

| # | 问题 | 现象 | 修复 | commit |
|---|---|---|---|---|
| 1 | `module load python/ 2>/dev/null` 只重定向了 stderr，但这台集群的 `module` 命令把错误信息写到了 stdout | `check`/`smoke` 输出里混入 `ERROR: Unable to locate a modulefile for 'python/'`（无害噪音，不影响功能） | 改成 `>/dev/null 2>&1` | `a1738fc` |
| 2 | **最关键的一个**：`sbatch` 会把提交的脚本复制到每个作业专属的 spool 目录（`/var/spool/slurm/d/job<id>/`）并执行那份拷贝，不是原始文件；脚本内部用 `readlink -f "${BASH_SOURCE[0]}"` 推导自己所在目录来找同目录下的 `oryza_besthit_damage_filter.py`，这个推导在 SLURM job 里跑的时候会指向错误的 spool 目录 | smoke test 实际报错：`python3: can't open file '/var/spool/slurm/d/job1807625/oryza_besthit_damage_filter.py': No such file or directory` | 提交时把已经在提交端正确解析好的 `PY_SCRIPT` 绝对路径通过 `--export="ALL,PY_SCRIPT=${PY_SCRIPT}"` 显式传进 job 环境，job 内不再重新推导 | `a1738fc` |
| 3 | 无 `submit all` / `local` 模式，用户批量跑/本地调试不方便 | — | 新增 `submit all`（自动发现已完成mapping的样本）和 `local`/`local all`（前台不走 sbatch，复用同一套 worker + `.finished` 跳过逻辑） | `83a0f6f` |

**重要经验**：bug #2 这一类"脚本自我定位"的坑，本地测试完全测不出来（本机
没有 sbatch，只能测 `run`/`local` 这两个不依赖 `sbatch` 拷贝脚本的路径），
必须在真实 SLURM 环境里才会暴露。以后改这个脚本、涉及"脚本找自己同目录下的
文件"这类逻辑时要格外小心，优先考虑用 `--export` 显式传值，而不是依赖
`$0`/`BASH_SOURCE` 在 job 里重新推导。

## 6. 已经拿到的实际结果：smoke test 分析（LV6000619499，前1000条 reads）

```
sample        input_reads  reads_with_alignment  reads_with_oryza_hit  kept_reads  rejected_nonoryza_better  rejected_no_oryza  unclassified_reads
LV6000619499  1000         1000                  956                   103         853                       44                 0
```

decision/reason 分布：
```
853 REJECT  nonoryza_better
100 KEEP    oryza_at_least_as_good
 44 REJECT  no_oryza_hit
  3 KEEP    oryza_only_no_competitor
```

**REJECT(nonoryza_better) 的输赢差距分布**（`best_nonoryza_adjusted_NM -
best_oryza_adjusted_NM`，负值=Oryza更差）：
```
差1: 541 (63%)   差2: 288 (34%)   差3: 23   差4: 1
```

**KEEP 的差距分布**：
```
打平(差0): 98个   Oryza真正赢(差1): 2个
```

**REJECT(nonoryza_better) 里排名第一的非Oryza物种 top20**（taxid，还没转学名，
需要用 names.dmp 查）：
```
1711249(37) 145626(31) 519541(29) 13415(26) 223100(23) 2790670(19) 126911(19)
103762(18) 4682(13) 322858(13) 3077(12) 3759(11) 192012(9) 9872(8) 56866(8)
1709936(8) 988163(7) 641091(7) 980011(6) 39329(6)
```
（这份是 taxid，之前有一份对应学名的旧版本 top1 结果，含 *Eucalyptus
dawsonii* / *Zizania palustris*（近缘野生稻属）/ *Bradybaena similaris*（陆生
蜗牛）/ *Geodia barretti*（深海海绵）/ *Chamaecyparis obtusa* / *Mikania
micrantha* / *Pogostemon cablin* / *Orobanche coerulescens* / *Bidens
hawaiensis* / *Closterium sp.*（藻类）/ *Brachypodium stacei*（禾本科模式种）
/ *Ammopiptanthus mongolicus* / *Diplosoma virens*（海鞘）/ *Prunus
yedoensis* / 两种蛾类/石蛾 / *Gossypium australe* / *Allium sativum* / *Chlorella
vulgaris* 等——**taxid列表和学名列表是否完全对应没有逐条核实过**，只是同一批
分析里前后两次查询，接手时建议重新跑一次、taxid和学名一起输出，别假设两份
列表严格对应）

### 6.1 这些数字说明什么

1. **输赢边界高度集中在"差1"**：63%的REJECT只差1个编辑距离，98%的KEEP是
   刚好打平——说明多数候选read的"是不是Oryza"这个信号本身就很弱，卡在最容易
   受参考基因组完整度影响的区间。
2. **这个margin的杠杆效应极大**：如果把判定从"打平"放宽到"差1也算"
   （margin=1），KEEP数会从103跳到644（1000条里的64%）——这不是微调，是量级
   变化，说明选哪个边界现在还不是"调参数"层面的事，而是要先解决参考基因组
   覆盖面的问题（见下）。
3. **Top1非Oryza物种名单里混着两类**：一类生物学上说得通（Zizania=野生稻近
   缘属、Brachypodium=禾本科模式种、其他植物），另一类很奇怪（深海海绵、
   陆生蜗牛、海鞘、蛾类）——这类不像会是真污染，更可能是数据库/acc2taxid层面
   的系统性问题，或者是低复杂度/保守序列在多个物种间普遍匹配导致的噪声。
4. **核心假说（待验证）**：现在"跟Oryza比对得好不好"基本等于"跟*O.
   sativa*这一个参考比对得好不好"（因为数据库里大概率没有/很少野生稻参考）。
   末端损伤校正只扣"古DNA特有的末端脱氨基信号"，不会、也不能区分"这条
   read有真实的生物学分化（野生型血统/古代原始驯化稻/别的水稻谱系）"和
   "这条read根本不是水稻"——这两种情况现在被同等地算作普通错配、拉高NM。
   如果古代样本的水稻谱系本身跟现代sativa参考有真实遗传分化，这些read会被
   系统性地判成"输给非水稻竞争者"，即便它们其实是水稻。

## 7. 待确认的开放问题（还没拿到答案，是当前最优先的工作）

### 7.1 acc2taxid 里 Oryza 属各物种的实际覆盖度

```bash
ACC2TAXID=/home/scratch/yinmt202607/db/asian_rice_panel_index/all_wgs_asian_irgsp.acc2taxid
NODES=/home/database/ref20250728/taxonomy_CPH/ncbi/20250530/nodes.dmp
NAMES=/home/database/ref20250728/taxonomy_CPH/ncbi/20250530/names.dmp

awk -F'|' '{gsub(/[\t ]/,"",$1); gsub(/[\t ]/,"",$2); if($2=="4527") print $1}' "$NODES" > /tmp/oryza_species_taxids.txt
awk -F'|' '$0 ~ /scientific name/{gsub(/[\t ]/,"",$1); gsub(/^[\t ]+|[\t ]+$/,"",$2); print $1"\t"$2}' "$NAMES" | grep -Ff /tmp/oryza_species_taxids.txt
while read -r tid; do
  n=$(awk -F'\t' -v t="$tid" '$3==t' "$ACC2TAXID" | wc -l)
  echo -e "${tid}\t${n}"
done < /tmp/oryza_species_taxids.txt
```
**这条命令给了用户，用户还没跑/没贴回结果。这是当前第一优先级的诊断。**

### 7.2 `db/3k/wild/` 140+ 野生稻组装的物种身份

```bash
ls /home/scratch/yinmt202607/db/3k/wild/ | head -20
ls /home/scratch/yinmt202607/db/3k/wild/ | wc -l
ls /home/scratch/yinmt202607/db/3k/       # 找有没有说明文档/manifest
head -3 /home/scratch/yinmt202607/db/3k/tmp/3kall_variety_map.tsv
for f in $(ls /home/scratch/yinmt202607/db/3k/wild/ | head -5); do
  sid="${f%%.*}"
  echo "=== $sid ==="
  grep -w "$sid" /home/scratch/yinmt202607/db/3k/tmp/3kall_variety_map.tsv
done
```
**同样给了用户、还没有结果**。这份 `file_path.md` 里只有一行目录清单提及
（`db/3k/wild/{SampleID}.transfer.merge.chr.fasta  # 140+野生稻/近缘种染色体级
组装基因组`），仓库里没有任何地方写明这140+个样本具体是哪个物种/亚群——
需要用 `3kall_variety_map.tsv`（3K RGP官方样本元数据表，含
VARIETY_INDEX/NAME/IRIS_ID）交叉查询才能确认。**注意**：3K Rice Genome
Project官方主体是3024份**栽培稻**（indica/aus/aromatic/japonica等亚群），
不是野生稻——这批 `wild/` 目录既然单独命名和别的3K数据分开放，大概率是
额外补充的野生近缘种outgroup，但没有manifest佐证之前不要假设。

### 7.3 `--oryza-taxids` 是否真的是 species rank

```bash
grep -P '^4529\t|^4530\t|^4536\t' /home/database/ref20250728/taxonomy_CPH/ncbi/20250530/nodes.dmp
```
脚本启动时只校验这三个 taxid 存在于 nodes.dmp，**不校验 rank 是不是
species**——如果 asian_rice_panel/IRGSP 的 acc2taxid 把 contig 指到了某个非
species 的 Oryza 内部节点，`species_of()` 会往上找到别的 species 节点甚至
找不到。还没确认过。

### 7.4 MD tag 覆盖率 / 实际内存占用

Smoke test 日志里如果有 `[warn] N alignments had no usable SEQ/MD` 这一行，
需要看 N 的数值——用户贴过的 smoke test 输出片段里没有看到这行（只贴到
`[load] acc2taxid...` 就中断了，后续 `[warn]`/`[summary]` 完整输出没有再贴
过一次完整版），**没有正式确认过 MD tag 覆盖率是否有问题**，只是根据最终
1000条里956条有 alignment 走完全流程、没有报错退出，间接推测大概率没问题。
接手时建议重新要一次完整 log。

内存：`JOB_MEM_MB` 默认16000（16GB）是本次工作开始时凭经验猜的，从未用
`seff <job_id>` 或 `wc -l` nodes.dmp/names.dmp/acc2taxid 校准过实际占用——
理论上 taxonomy dump（几百万taxid，3个dict）+ acc2taxid dict 全量常驻内存
是这个脚本的内存大头，量级可能到GB级别，具体数字取决于这两个文件的真实大小，
需要 `seff` 或 `wc -l` 实测。

## 8. 后续三条工作线（讨论过，还没定最终方案）

这三条是用户提出的下一步方向，讨论了可行性和优先级，但都**还没有开始实际
写代码/跑分析**，需要先决定做哪个/什么顺序。

### 8.1（最优先，和第7节强相关）Oryza 参考基因组数据库整理/扩充

现状：`file_path.md`（main分支）里明确标注多处 ⚠️，说明参考基因组的整体情况
"还是很乱"（用户原话）：
- `db/16/`（资源组A，NCBI datasets 16基因组+Liftoff注释）和 `db/3k/`（资源组B，
  3K Rice Genome Project数据）是**两个不同来源**，关系没有理清。
- `asn720data/`（720份现代/近现代品种PLINK面板）跟16个angkor古代样本的关系
  "尚未最终确认"（是否包含、是否capture panel来源）。
- `db/3k/wild/` 这140+个野生稻/近缘种组装身份不明（见7.2）。
- 已经有一个独立的 `results/07.wild_rice_alignment/`（minimap2 asm10 预设，
  野生稻组装比对到IRGSP，用于paftools.js call变异，标注"进行中"）——这是
  另一条已经在跑的、跟competitive mapping平行的分析线，不要跟best-hit阶段
  混淆，但如果这些wild genome的minimap2比对已经产出了变异，可能可以复用来
  确认这些野生稻组装的物种/谱系身份。

**建议的下一步（不是最终决定，需要用户拍板）**：
1. 先做完第7.1/7.2节的诊断，拿到"acc2taxid里野生稻覆盖度"和"db/3k/wild/
   140+样本的物种身份"两组硬数据。
2. 如果确认这批野生稻组装可用，把它们**作为第132个（或更多个）competitive
   mapping数据库**加进去——这比重跑全部131个库便宜得多，只需要给已经跑完
   mapping的样本再补一次新数据库的比对，合并进现有的 `by_sample/*.name_sorted.bam`，
   然后重跑 besthit（不需要重新mapping其他131个库）。
3. 扩库之后再回头看第6节那个"打平/差1"的margin问题是否自然缓解——如果多数
   "惜败"read其实是能在扩库后找到更匹配的野生稻/近缘种参考，margin问题可能
   不需要再单独调参就基本解决。

### 8.2 时间序列选择信号扫描（selection scan）

背景（用户和同门的聊天记录）：想找"某个位点等位基因频率随时间**定向、快速**
变化"的信号，同门评价"很容易看""碰运气"。

**可行性评估**：
- 全项目16个古代样本，年代跨度约公元1160-1981年（近800年）——分到各个
  时间段可用个体数很少，是这类分析最大的硬约束。
- "碰运气/抽的多了总能碰上"这个说法点出了真实的方法论风险：不做多重检验
  校正/零分布模拟，genome-wide扫描出来的"候选位点"里假阳性会很多，纯粹的
  遗传漂变在小样本下也能产生看起来"剧烈"的频率波动。
- **正面信息**：项目之前已经搭建过 NGSNGS 模拟环境（git log:
  "NGSNGS simulation environment setup"，commit `bf107e84`），这个可以用来
  生成"纯漂变、无选择"下的零分布，给观测到的频率变化做显著性检验，而不是
  只看"变化快不快"这种没有统计标定的判断。
- **建议**：与其一上来做无先验的genome-wide扫描，优先聚焦在已经有生物学
  先验的候选基因（比如8.3提到的57个开花基因，本身就是驯化/生态型相关性状
  候选），信噪比会更好，也更容易向外解释结果。
- **这条线依赖第8.1**：genotype矩阵/等位基因频率的计算，前提是有干净、
  一致的比对/genotype pipeline，参考基因组整理不清楚会直接影响这里的
  可信度。

**还没做的事**：具体用哪些位点/genotype来源（是besthit过滤后的read直接
call变异，还是走已有的 `results/02.irgsp/` 主流程产出）、零分布模拟的具体
参数、时间分段方式，都没有讨论到细节，需要下一步单独设计。

### 8.3 57个开花基因：探针内部/两侧变异判定生态型（旱稻 vs 水稻）

背景：57个开花基因清单在 `db/gene/flower_gene.txt`（含 Level/Group/MSU_id/
RAPDB_id/Name），探针设计覆盖基因两侧（`db/gene/flower_gene.flank1kb.bed`，
±1kb）和基因内部。想法是找基因内部/附近的"典型插入"（结构变异）来判断
古代样本是旱稻还是水稻生态型。

**方法论上是站得住的**：这不是拍脑袋——项目自己已经有先例，git log里查过
DTH8/Ghd8经典大片段缺失（3024份现代品种里5.4%携带这个缺失，`commit
75f7e81`），这正是"用已知SV有无判定表型/生态型"这套方法在本项目里的先例。
现代水稻遗传学里也确实有一批性状对应已知的causal SV（比如落粒基因sh4、
半矮秆sd1、芒退化An-1等），不是凭空设计。

**最大的坑，项目自己已经踩过一次**：git log里同一批commit还有一条
"DTH8 breakpoint zero-coverage result"（commit `bf107e84`）——**哪怕是这么一个
研究得很透的明星基因，古DNA在这个位点的覆盖度都不够看**，这个已经提前预警了
"低覆盖度能不能可靠call出SV有无"这个问题在本项目里是真实存在的，不是理论
担忧。盲目对57个基因逐一去找"典型插入"，大概率在多数基因上撞到同样的墙。

**建议的下一步**：
1. **先做覆盖度QC**，不要直接冲进去分析。用现有 `results/02.irgsp/`（或
   `results/igv_package/bam_q30/`）BAM，对57个基因（含±1kb flank）逐样本
   算深度，看到底有几个基因、几个样本的覆盖度真的够格尝试call SV有无。这个
   半天工作量就能出结果，能立刻告诉你这条路现实不现实，比直接逐基因分析
   划算得多。
2. 确认这57个基因里**哪些在文献里真的有已知的旱稻/水稻区分性SV**（不是
   每个开花基因都有已知causal variant）——`flower_gene.txt` 的 Level/Group
   字段可能已经有部分线索，需要人工核实。
3. 覆盖度QC + 已知SV筛选两步都做完之后，能分析的基因数量可能远小于57，
   到时候再具体设计怎么从探针捕获的reads里判定SV有无（比如split-read/
   discordant-pair证据，还是简单的深度骤降，取决于具体是哪种SV）。

## 9. 给下一个接手者的具体待办（按优先级）

1. 【诊断，最优先】跑第7.1节的acc2taxid Oryza覆盖度诊断命令，拿到结果。
2. 【诊断】跑第7.2节的 `db/3k/wild/` 物种身份诊断命令，拿到结果。
3. 【诊断，可选但建议做】第7.3（taxid rank校验）、7.4（MD tag覆盖率完整
   log、`seff`内存实测）。
4. 根据1-2的结果决定：野生稻参考要不要扩库（见8.1），要扩的话具体怎么建
   索引、怎么接入现有的132个数据库体系。
5. 扩库（如果做）之后重新评估best-hit的KEEP/REJECT margin（见第6节的
   "打平vs差1"讨论），再决定要不要在 `oryza_besthit_damage_filter.py` 里
   加一个margin参数（当前是硬编码 `<=`，没有margin概念）。
6. 对剩余13个还在mapping队列里的样本，随mapping进度用 `submit_oryza_besthit.sh
   submit all` 陆续跑besthit。
7. 8.2（selection scan）和8.3（57基因SV判生态型）都还没真正开始，8.3的
   覆盖度QC是个独立、低成本、能快速出结论的子任务，可以和上面几条并行推进；
   8.2依赖更完整的genotype pipeline，优先级排后面。
