# Oryza competitive mapping / best-hit 接手说明

更新时间：2026-08-05

## 1. 当前目标

候选水稻 reads 已经由 shotgun、panel1、panel2 三份 FASTQ 合并完成。当前正在把每个样本分别比对到 WGS 真核数据库、亚洲水稻 panel 和 IRGSP。下一阶段需要根据所有比对结果执行 best-hit，排除更像非水稻物种的 reads。

## 2. 关键路径

### GitHub

仓库：

https://github.com/Inmpain/rice_adna_pipeline

当前分支：

`codex/oryza-competitive-mapping`

mapping 脚本：

`scripts/oryza_besthit/submit_oryza_competitive_mapping.sh`

测试脚本：

`scripts/oryza_besthit/test_submit_oryza_competitive_mapping.sh`

### 本机参考代码

原 best-hit 项目：

`/Users/inmpain/github/aeDNA_popgen`

WGS mapping/动态内存参考：

`/Users/inmpain/github/Pipeline_snakemake/new_single_multi`

本机 GitHub 工作副本：

`/Users/inmpain/Documents/angkor/rice_adna_pipeline_publish`

### 服务器

执行脚本：

`/home/scratch/yinmt202607/gene/scripts/submit_oryza_competitive_mapping.sh`

候选 FASTQ：

`/home/scratch/yinmt202607/gene/results/oryza_candidates_combined`

格式：

`<sample>.oryza_candidates.combined.fastq.gz`

mapping 输出根目录：

`/home/scratch/yinmt202607/gene/results/oryza_competitive_mapping`

各数据库 BAM：

`.../bam_by_database/<sample>/<sample>.<database>.bam`

每样本最终合并 BAM：

`.../by_sample/<sample>.competitive.name_sorted.bam`

日志：

`.../logs`

提交记录：

`.../submissions`

串行运行状态：

`.../series`

## 3. 数据库

WGS 真核数据库：

`/home/database/ref20250728/cph_euk/wgs_eukaryota.1.fas.gz`

至：

`/home/database/ref20250728/cph_euk/wgs_eukaryota.129.fas.gz`

亚洲水稻 panel：

`/home/scratch/yinmt202607/db/asian_rice_panel_index/asian_rice_panel.fa`

IRGSP：

`/home/scratch/yinmt202607/db/asian_rice_panel_index/irgsp_bt2idx`

合并 accession-taxid：

`/home/scratch/yinmt202607/db/asian_rice_panel_index/all_wgs_asian_irgsp.acc2taxid`

NCBI taxonomy：

`/home/database/ref20250728/taxonomy_CPH/ncbi/20250530`

目标 Oryza taxid：

- 4529：Oryza rufipogon
- 4530：Oryza sativa
- 4536：Oryza nivara

## 4. 当前 mapping 流程

每样本比对131个数据库：

- 129个 WGS shard
- asian_rice_panel
- IRGSP

Bowtie2 参数：

```bash
-k 100
-L 22
-i S,1,1.15
--mp 1,1
--rdg 0,1
--rfg 0,1
--score-min L,0,-0.1
--no-unal
```

截至 2026-08-05，`by_sample/` 下已经跑完并 name-sort 完成（`.finished` 标记齐全）的样本：

- LV6000619499
- LV6000619917
- LV6000620016

其余样本仍在串行队列里（见 `.../series`）。

## 5. best-hit / aDNA 损伤校正过滤（当前阶段）

**没有**复用 `aeDNA_popgen` 的通用 ngsLCA 风格层级分类器
（`besthit_competitive.py`）——那套做的是"沿谱系逐级上溯找最深胜出节点"，
不区分古 DNA 末端损伤和真实突变，也不是本项目需要的"Oryza vs 非 Oryza"
二元竞争。这里是专门为本项目重新写的一套逻辑，核心是：**扣除末端损伤后的
`adjusted_NM`，Oryza 最优 hit 必须与非 Oryza 最优 hit 打平或更好，才 KEEP**。

脚本落地在本仓库：

- `scripts/oryza_besthit/oryza_besthit_damage_filter.py` —— 主统计脚本
- `scripts/oryza_besthit/submit_oryza_besthit.sh` —— SLURM 提交脚本（check /
  smoke / submit / merge 四种模式，风格与 `submit_oryza_competitive_mapping.sh`
  一致）

### 5.1 核心算法

每条 alignment 计算：

- `NM`：bowtie2 的 NM tag（原始编辑距离）
- `substitution_count`：从 MD/CIGAR 精确统计的真实替换数（不含 indel），靠
  `pysam.get_aligned_pairs(with_seq=True)` 拿到逐位点的 (ref_base, read_base)
- `terminal_damage_count`：只在读长两端 `--damage-window`（默认5）范围内、且
  **换算回 read 自身 5'→3' 方向后**符合古 DNA 损伤特征（5'端 C→T，3'端 G→A）
  的替换才计入。reverse-strand alignment 的换算不是简单倒序——BAM 里 SEQ 已经
  是原始 read 的反向互补（SAM 规范），换算回 read 自身方向需要**倒序 + 互补**
  两步一起做，否则会把参考链上的 G→A 错当成非损伤位点漏掉。已经用手工构造的
  reverse-strand 测试 BAM 验证过这一步（见下方"已完成的本地验证"）。
- `adjusted_NM = NM - terminal_damage_count`（indel 永远不扣）

按 species taxid（BAM reference → acc2taxid → 沿 nodes.dmp 上溯到最近的
rank=="species" 祖先）分组，每个 species 只保留最优 alignment，排序优先级：
`adjusted_NM` 最小 → `NM` 最小 → `substitution_count` 最小 → `AS` 最大 →
reference 名称（保证稳定/可重复）。非 Oryza species 取前 `--top-n`（默认10）
写入审计表；Oryza（`--oryza-taxids`，默认 4529/4530/4536）不占非 Oryza 名额，
命中就必录。

判定：

- 没有 Oryza hit → REJECT (`no_oryza_hit`)
- 有 Oryza hit、没有非 Oryza hit → KEEP (`oryza_only_no_competitor`)
- `best_oryza.adjusted_NM <= best_nonoryza.adjusted_NM` → KEEP
  (`oryza_at_least_as_good`)
- 否则 → REJECT (`nonoryza_better`)

不用 MAPQ（bowtie2 `-k 100` competitive mapping 下不可靠）。BAM 没有去重、也
不需要在这一步去重（同一 read 的多行来自131个数据库 × `-k 100`，不是 PCR
duplicate）。

### 5.2 输出（每样本，`.tmp` 写完原子 rename；只有一致性检查通过才写 `.finished`）

- `<sample>.besthit.top10_species.tsv.gz`：审计表，一条 read 对应多行（非
  Oryza top10 + 必录的 Oryza 行），列含 `top10_rank`/`is_always_included`
- `<sample>.oryza_filter.decisions.tsv.gz`：每 read 一行的最终判定
  （`best_nonoryza_*` / `best_oryza_*` / `decision` / `reason`）
- `<sample>.besthit_oryza.fastq.gz`：从候选 FASTQ 里按 KEEP read name 抽出，
  每条只输出一次，序列/质量值原样保留
- `<sample>.summary.tsv`：单样本一行（`input_reads` 等于
  `kept + rejected_nonoryza_better + rejected_no_oryza + unclassified_reads`，
  脚本内部会自查，不一致直接非零退出，但不影响这个文件本身被写出，方便debug）

多样本一致性检查见脚本 docstring；`besthit_summary.tsv`（多样本汇总）由
`submit_oryza_besthit.sh merge` 事后拼接每个样本自己的 `.summary.tsv`
生成，刻意不在 Python 脚本里直接写共享文件，避免并行 SLURM 作业互相竞争。

### 5.3 已完成的本地验证

在没有服务器数据的情况下，用手工构造的最小 BAM（`samtools view -bS` +
`sort -n`，5条参考、5条候选 read）跑通了整条流程，重点验证：

- reverse-strand 末端损伤校正确实会**反转最终判定**：构造了一条 read，原始
  NM 下"非 Oryza 更优"应判 REJECT，扣除末端损伤后 `adjusted_NM` 打平，正确
  判成 KEEP（`oryza_at_least_as_good`）——这是本脚本存在的核心原因，必须验证。
- 一条读在候选 FASTQ 里、但 BAM 里完全没有 alignment 的 read，正确计入
  `unclassified_reads`，不出现在 `decisions.tsv`。
- `input_reads == kept + rejected_nonoryza_better + rejected_no_oryza +
  unclassified_reads` 自检通过；FASTQ 抽取数等于 `kept_reads`。
- `--limit-reads`（smoke test）路径**不**写 `.finished`（早期版本有个 bug 会
  在 smoke 模式下也写 `.finished`，已修复并重新验证——这个坑很关键，写错了会
  让批量提交脚本误判某样本"已完成"而跳过正式全量跑）。
- 重复跑同一个样本（模拟 SLURM 任务失败后重提交）不会因为 `.tmp` 残留或
  rename 冲突而出错。
- `submit_oryza_besthit.sh run`（内部 worker 模式，绕开 sbatch 直接调用）和
  `merge` 模式也跑通了完整路径。

没有验证的：真实服务器 BAM（131 库竞争比对后的实际 MD/AS 内容、acc2taxid 对
asian_rice_panel/IRGSP contig 的覆盖率）、以及 SLURM 提交本身（本机没有
sbatch）。下面几点在服务器上第一次跑时必须确认：

1. **taxid 是否都是 species rank**：`--oryza-taxids` 默认 4529/4530/4536，
   脚本会在启动时校验这三个 taxid 都存在于 `nodes.dmp`，但不校验 rank 是不是
   `species`——如果 asian_rice_panel/IRGSP 的 acc2taxid 把 contig 指到了某个
   非 species 的 Oryza 内部节点，会在 `species_of()` 里向上找到别的 species
   节点甚至找不到，建议先 `grep -P '^4529\t|^4530\t|^4536\t' nodes.dmp` 确认。
2. **MD tag 是否存在**：脚本假设 bowtie2 输出的 BAM 带 MD tag；如果没有，
   `substitution_count`/`terminal_damage_count` 会退化为 `NA`/`0`
   （`adjusted_NM=NM`，不报错，但也拿不到损伤校正），脚本会在结尾打印
   `[warn] N alignments had no usable SEQ/MD`——先用 smoke test 检查这个数字
   是否接近 0。
3. **acc2taxid 是否覆盖 asian_rice_panel / IRGSP 的 contig 命名**：`[load]`
   那行日志会报告"多少 reference contigs 没有 acc2taxid 命中"，数值应该很小
   （只有 129 个 WGS shard 的边角 contig 缺失是正常的，如果 asian_rice_panel/
   IRGSP 本身缺失就是 acc2taxid 合并有问题）。

### 5.4 用法

```bash
cd /home/scratch/yinmt202607/gene/scripts   # 或脚本所在目录

# 1. 校验路径/python3+pysam/SLURM partition，不提交任何作业
./submit_oryza_besthit.sh check

# 2. 先用一个样本的前1000条 reads 做 smoke test（输出到独立的 smoke_test 目录，
#    不写 .finished，不影响正式产出）
./submit_oryza_besthit.sh smoke LV6000619499

# 3. smoke test 结果没问题后，对已完成 mapping 的样本跑全量（自动跳过已有
#    .finished 的样本，可以放心重复执行）
./submit_oryza_besthit.sh submit LV6000619499 LV6000619917 LV6000620016

# 4. 所有作业跑完后，拼接每个样本的 summary.tsv
./submit_oryza_besthit.sh merge
# -> /home/scratch/yinmt202607/gene/results/oryza_competitive_mapping/besthit/besthit_summary.tsv
```

若 `check` 报 `import pysam` 失败，需要先在这个 `python3` 所在环境装 pysam
（`pip install --user pysam` 或专门建一个 conda env）。

### 5.5 下一步

1. 跑 `check` → `smoke`，人工核对 smoke test 的 `[warn]`/`[summary]` 输出（MD
   缺失比例、Oryza/非 Oryza 占比是否合理）。
2. 对 3 个已完成 mapping 的样本跑全量 `submit`，`merge` 后看
   `besthit_summary.tsv`：`unclassified_reads` 占比、`rejected_no_oryza` 占比
   是否符合预期。
3. 单样本结果符合预期后，再对剩余样本批量跑（mapping 串行队列跑完一个样本
   就可以对它单独 `submit`，不需要等全部样本 mapping 完成）。
4. 视情况检查 `top10_species.tsv.gz` 里非 Oryza 排名靠前的物种，确认是否有
   数据库 shard 或 acc2taxid 层面的系统性问题。
