# Oryza competitive mapping / best-hit 接手说明

更新时间：2026-08-08（第三次更新：发现并修正了`asian_rice_panel.acc2taxid`
里`np7`(IRGSP/Nipponbare)的taxid错标问题，详见新增的"⚠️0.5 taxid错标"节。
第二次更新记录：acc2taxid Oryza覆盖度诊断结果出来后，证伪了最早"数据库
野生稻缺失"的假设，改用"损伤窗口固定长度 vs 读长"这个新证据更充分的假说，
相应重排了第7/8/9节的优先级）

> 读这份文档前提：你没有服务器直接执行权限，所有服务器端命令的输出都需要
> 人类协作者手动跑了贴回来。不要假设你能直接看到服务器文件系统。

## 0. 现状一句话总结

Best-hit 脚本（`oryza_besthit_damage_filter.py` + `submit_oryza_besthit.sh`）
已经在服务器上实际跑通（3个样本里的1个做过1000-read smoke test，其余2个mapping
已完成但besthit还没跑全量）。脚本本身没有已知 bug 了（过程中发现并修复了3个，
见5.6）。**当前卡点不是代码**。

**方法论上走过一轮弯路，已经纠正，写在这里避免下一个人重复踩**：最早怀疑是
"acc2taxid 数据库里 Oryza 参考基因组覆盖不够、大概率只有 sativa"——**这个
假设已经用实测数据证伪**（见7.1）：Oryza 属18个已知种在数据库里普遍有几十到
几百条 contig，*rufipogon*(717条) 反而比 *sativa*(96+15=111条) 覆盖更好。
真正吻合观测数据的假说是**末端损伤校正窗口（`--damage-window`默认5bp）相对
读长是固定的，读长越长、两端5bp之外"没资格被校正"的中段就越大，累积一个
额外错配（无论是真实SNP、非末端损伤、还是测序错误）的概率也越高**——
smoke test里"差1个编辑距离惜败"的541条read，均值/中位数长度（~68bp）明显
比"打平守住KEEP"的103条read（~62bp）更长，方向和这个假说完全吻合（见6.2）。
**这个假说还差最后一步验证**：这批样本有没有做过标准的损伤衰减曲线分析
（末端N个位置的C→T/G→A频率随距离衰减的图，mapDamage风格），衰减到背景水平
大概在多少bp——如果实际损伤延伸远超5bp，说明`--damage-window`设得太窄，
这是下一个人接手后第一件要确认的事（见7.5）。

## 0.5 ⚠️新发现（2026-08-08）：`np7`(IRGSP/Nipponbare) 在 acc2taxid 里 taxid 标反了

用户和GPT的独立讨论（完整过程见`docs/asian_rice_panel_reference_design_
conversation.md`）里，用染色体长度+逐条染色体MD5比对，确认了
`asian_rice_panel.fa`里的`np7.Chr1-12`与`irgsp.fa`的`chr01-12`**完全一致
(IDENTICAL)**——即`np7` = Nipponbare/IRGSP-1.0 = *O. sativa*。但
`asian_rice_panel.acc2taxid`里`np7.*`(含Chr1-12、ChrUn、ChrSy、ChrM、ChrC
共16行)**全部被标成了4529(*O. rufipogon*，普通野生稻)**，应为4530
(*O. sativa*)。

**对已有分析的影响**：
- **不影响besthit的KEEP/REJECT二元判定**——`--oryza-taxids`默认
  `"4529 4530 4536"`，4529和4530都在白名单里，np7的reads无论被算成哪个
  taxid，仍然会被判定为"目标Oryza"，不会被误判成非Oryza。
- **但影响第6.1节的"各物种contig数"统计表**——那张表用来证伪"数据库野生稻
  缺失"假说的关键论据是"*sativa*只有111条、*rufipogon*有717条，sativa反而
  偏低"。np7的16行如果原本被错误计入了717这个rufipogon桶，说明**修正后
  *sativa*的真实contig数会比111更低（约95条）、*rufipogon*会比717更低
  （约701条）**——需要在服务器跑完修正脚本后重新统计这张表，确认"sativa
  偏低"这个结论的方向和幅度有没有变化（方向大概率不变，因为16条的调整量
  相对717/111的量级不大，但具体数字要更新）。

**修正方案**：`scripts/oryza_besthit/fix_asian_rice_panel_taxid.sh`——
本地已用合成数据测试过dry-run和--apply两种模式，逻辑是：
1. 重新给`asian_rice_panel.acc2taxid`里12个基因组(np7/mh63/X24_kas/azu/arc/
   liuxu → 4530；7个G25_ruf_W* → 4529)赋正确taxid
2. 在`all_wgs_asian_irgsp.acc2taxid`里**补丁式替换**这12个基因组对应的行
   （不是重新做三源合并——WGS+asian_rice_panel+IRGSP的原始合并脚本没有
   留存，不知道怎么拼的，所以只精确替换我们确认有问题的这一块，不动其余
   119个WGS shard和IRGSP各自的部分）

⚠️**只有np7做过这种逐染色体MD5的独立验证**。其余11个基因组(mh63/X24_kas/
azu/arc/liuxu/7个G25_ruf_W*)的物种归属是根据文件命名和来源判断的，沿用了
讨论里的结论，**没有再逐一验证**——如果后续发现某个命名对不上，这个脚本
和这份记录都需要跟着改。`X24_kas`尤其存疑（讨论里明确提到"如果确实是
Kasalath"这个前提未经验证）。

**待办**：
```bash
# 1. 先 dry-run，确认改动符合预期(不写文件)
curl -fsSL -o fix_asian_rice_panel_taxid.sh \
  https://raw.githubusercontent.com/Inmpain/rice_adna_pipeline/codex/oryza-competitive-mapping/scripts/oryza_besthit/fix_asian_rice_panel_taxid.sh
bash fix_asian_rice_panel_taxid.sh /home/scratch/yinmt202607/db/asian_rice_panel_index

# 2. 确认无误后加 --apply 真正落盘(会自动备份原文件)
bash fix_asian_rice_panel_taxid.sh /home/scratch/yinmt202607/db/asian_rice_panel_index --apply

# 3. 重新统计第6.1节的物种contig数表，更新这份文档
```

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
（12个基因组：np7/mh63/X24_kas/azu/arc/liuxu 6个O. sativa代表品种 + 7个
G25_ruf_W*的O. rufipogon野生稻，选取理由/覆盖的遗传谱系/质量评估见
`docs/asian_rice_panel_reference_design_conversation.md`）

IRGSP：`/home/scratch/yinmt202607/db/asian_rice_panel_index/irgsp_bt2idx`

合并 accession→taxid：`/home/scratch/yinmt202607/db/asian_rice_panel_index/all_wgs_asian_irgsp.acc2taxid`
（标准 NCBI 三列格式：accession / accession.version / taxid。**⚠️2026-08-08
发现np7标签错误，见0.5节，还有独立的`asian_rice_panel.acc2taxid`和
`irgsp.acc2taxid`两个源文件，只检查过前者**）

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
- `scripts/oryza_besthit/fix_asian_rice_panel_taxid.sh` —— taxid错标修正脚本（见0.5节，2026-08-08新增）

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
| 4 | `asian_rice_panel.acc2taxid` 里 `np7`(IRGSP/Nipponbare) 被错标成 4529(rufipogon)，应为4530(sativa) | 用户用染色体MD5比对发现 | 见0.5节，`fix_asian_rice_panel_taxid.sh`，**待服务器实跑验证** | 待补 |

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

**REJECT(nonoryza_better) 里排名第一的非Oryza物种 top20**（taxid → 学名，
已用 names.dmp 逐条核实）：
```
1711249  37  Eucalyptus dawsonii        桉树
145626   31  Bradybaena similaris       陆生蜗牛
519541   29  Geodia barretti            深海海绵
13415    26  Chamaecyparis obtusa       扁柏(针叶树)
223100   23  Orobanche coerulescens     列当(寄生植物)
2790670  19  Closterium sp. NIES-67     新月藻(藻类)
126911   19  Ammopiptanthus mongolicus  沙冬青(灌木)
103762   18  Zizania palustris          北美菰(野生稻近缘属，唯一说得通的)
4682     13  Allium sativum             大蒜
322858   13  Diplosoma virens           海鞘(海洋无脊椎动物)
3077     12  Chlorella vulgaris         小球藻
3759     11  Prunus yedoensis           樱花
192012    9  Mikania micrantha          薇甘菊(藤本)
9872      8  Odocoileus hemionus        骡鹿(哺乳动物)
56866     8  Luffa acutangula           丝瓜
1709936   8  Eurotiomycetes sp.         真菌
988163    7  Tiliacea citrago           蛾类
641091    7  Trididemnum clinides       海鞘
980011    6  Bidens hawaiensis          鬼针草
39329     6  Lavandula angustifolia     薰衣草
```
这份名单在生物学上没有一致性——深海海绵、两种海鞘、陆生蜗牛、北美骡鹿，跟
东南亚考古水稻沉积样本对不上号，除了 *Zizania palustris*（野生稻近缘属）之外
基本看不出真实污染的规律。**这是判断"惜败"read性质的关键证据**：如果真是
参考基因组缺失导致的误判，输给的应该是生物学上说得通的近缘种；现在这个
名单更像是"短片段靠概率撞上129个shard里某个不相关物种的某段序列"，见6.2。

### 6.1 acc2taxid 覆盖度诊断结果（第7.1节问题，已有答案）

⚠️**2026-08-08更新：下表统计于`np7`taxid错标(见0.5节)被修正之前，需要在
服务器跑完`fix_asian_rice_panel_taxid.sh --apply`后重新统计**。np7有16行
（Chr1-12+ChrUn+ChrSy+ChrM+ChrC）原本被错误计入了下面的4529(rufipogon)桶，
修正后*O. sativa*(4530)的真实数字会比111更低（约95）、*O. rufipogon*
(4529)会比717更低（约701）。**"sativa覆盖度偏低"这个结论的方向大概率不变
（16条的调整量相对717/111的量级不大），但下面的具体数字已知有误，重新统计
之前不要直接引用**。

Oryza 属18个已知种在 `all_wgs_asian_irgsp.acc2taxid` 里的 contig 数（用
`nodes.dmp` 里 parent==4527 找到全部种级taxid，逐个数 acc2taxid 里精确匹配
的行数；*sativa* 另外加了亚种级taxid 39946/39947/1050722/1080340/2998809
的计数，只有39947=japonica有15条，其余0条）：

```
4528  Oryza longistaminata    905
4529  Oryza rufipogon         717   <- 野生稻，覆盖很好 (⚠️待修正后重新统计)
4530  Oryza sativa             96 (+15 japonica = 111)  (⚠️待修正后重新统计)
4532  Oryza australiensis      60
4533  Oryza brachyantha        60
4534  Oryza latifolia          78
4535  Oryza officinalis        91
4536  Oryza nivara             26   <- 野生稻，有但不算多
4537  Oryza punctata           61
4538  Oryza glaberrima         60  (非洲栽培稻)
40148 Oryza glumipatula        47
40149 Oryza meridionalis       12
52545 Oryza alta               24
63629 Oryza minuta            189
65489 Oryza barthii            51
77588 Oryza coarctata         450
110451 Oryza schlechteri      119
127571 Oryza malampuzhaensis  198
```

**结论：这个假说被证伪**（⚠️此结论基于修正前的数字，方向大概率仍成立，
待重新统计后确认）。*sativa* (111条) 不但不缺，反而是覆盖度偏低的
那一档——*rufipogon*(717)、*longistaminata*(905)、*coarctata*(450) 等好几个
种覆盖度都远高于sativa。"数据库里野生稻/近缘种缺失"不是"差1惜败"这个模式
的解释，第8.1节原本"优先扩库"的建议已经不成立，改成"优先查损伤窗口"
（见6.2、7.5、8.1）。

### 6.2 差1惜败的read为什么会输：不是短，是长（新发现，推翻了最早的猜测）

最早怀疑"惜败"read是不是偏短（短片段信息量少、更容易被噪声压过），实测
结果**方向相反**：

```
                均值(bp)  中位数(bp)  n
差1惜败(REJECT)   ~68.5      ~68      541
KEEP(打平守住)    ~62.3      ~61      103
```

（原始 `top10_species.tsv.gz` 是一条read多行的审计表，第一次按行统计长度会
被"这条read命中了多少个物种"放大，必须先按read_name去重再统计——这是个
真实踩过的坑，两组的原始"行数"分别是541条read共约6778行、103条read共约
1225行，去重后总数才分别对上541/103，接手时如果要重跑这类分析记得先去重）

**长而不是短的read更容易惜败，机制上的解释**：`--damage-window` 固定只看
读长两端各5bp，超出这个范围的任何替换（不管是没修完的非末端损伤、真实的
生物学SNP、还是测序错误）都没资格被 `adjusted_NM` 扣除，会原样拉高分数。
读长越长，两端5bp加起来占全长的比例越小，"没资格被校正"的中段就越大，
累积一个额外错配、从而在adjusted_NM上被非Oryza竞争者反超的概率也就越高。
这个方向和实测数据完全吻合。

**这个假说还差一步验证**：损伤窗口5bp这个默认值有没有校准过？如果这批样本
的真实损伤信号（末端C→T/G→A频率）衰减到背景水平的距离明显超过5bp（比如
mapDamage风格分析常见的10-15bp甚至更远，取决于样本降解程度），那5bp的窗口
本身就设窄了，会系统性地对长读长不利——这不是"要不要加margin"的问题，是
损伤校正窗口参数本身需要重新校准。有没有现成的损伤曲线分析（不管是这批
angkor样本还是这个项目别的样本）？如果没有，可以从BAM直接算一个粗略版本
（末端N个位置的C→T/G→A频率 vs 距末端距离），不需要跑完整mapDamage2。
见7.5的具体诊断命令。

## 7. 待确认的开放问题

### 7.1 ~~acc2taxid 里 Oryza 属各物种的实际覆盖度~~ 已回答，见6.1

结论：Oryza属18个种普遍有几十到几百条contig，*sativa*(111条)不缺，反而是
偏低的一档。这个问题已解决，不再是待办。⚠️2026-08-08：数字待修正np7标签后
重新核实，见0.5/6.1节，结论方向预期不变。

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
**给了用户、还没有结果**。这份 `file_path.md` 里只有一行目录清单提及
（`db/3k/wild/{SampleID}.transfer.merge.chr.fasta  # 140+野生稻/近缘种染色体级
组装基因组`），仓库里没有任何地方写明这140+个样本具体是哪个物种/亚群——
需要用 `3kall_variety_map.tsv`（3K RGP官方样本元数据表，含
VARIETY_INDEX/NAME/IRIS_ID）交叉查询才能确认。**注意**：3K Rice Genome
Project官方主体是3024份**栽培稻**（indica/aus/aromatic/japonica等亚群），
不是野生稻——这批 `wild/` 目录既然单独命名和别的3K数据分开放，大概率是
额外补充的野生近缘种outgroup，但没有manifest佐证之前不要假设。**另见
`docs/RESEARCH_ROADMAP.md`第2节C——这140+个组装是否等同于Guo et al. 2025
pangenome论文提供的145个组装，是main分支也在追的同一个开放问题。**

**优先级下调**：这条本来是基于"数据库野生稻缺失"假说排的第一优先级，
现在7.1已经证明acc2taxid里Oryza属覆盖面本身就很宽，这条不再紧急，可以
按需做，不是阻塞项。

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

### 7.5（当前最优先）`--damage-window=5` 是否偏窄——需要实测损伤衰减曲线

见0/6.2的推理：差1惜败的read系统性比KEEP的read长，机制上最可能的解释是
固定5bp窗口对长读长不利。这个假说要成立，前提是这批样本的真实损伤信号确实
延伸超过5bp。**先问用户/查项目历史**有没有现成的损伤曲线分析（mapDamage2
或类似工具跑过的5'C→T / 3'G→A频率-vs-位置图，针对angkor这批样本或者项目
里任何同批次测序的样本）。如果没有，可以从besthit的smoke test BAM直接估算
一个粗糙版本（只需要 KEEP 结果里"确认是水稻"的read，避免非水稻read混进来
干扰损伤信号统计）：

```bash
# 思路：besthit_oryza.fastq.gz 里的read名字，去原始 by_sample BAM 里找到它们
# 的primary alignment，统计末端N个位置(比如N=15，覆盖比默认5bp更远的范围)
# 的错配是不是C->T(5'端)/G->A(3'端)富集，并且看这个富集程度随距离衰减到
# 背景水平大概在第几个bp——这个可以另写一个小脚本，不在
# oryza_besthit_damage_filter.py现有功能范围内，需要新写。
```
如果衰减距离明显超过5bp（比如10-15bp），建议把 `--damage-window` 调大，
重新跑 smoke test 看"差1惜败"的541条里有多少会因此翻盘成KEEP，再判断这
是不是6.2这个模式的完整解释，还是只能解释一部分。

**2026-08-07进展**：4个样本(LV6000619499/619917/620016/620032)已跑完全量
besthit，留存率(kept/input)分别为5.00%/5.00%/4.37%/4.92%，四个样本高度
一致，没有出现某个样本被5bp窗口明显"冤枉"的迹象——这是"5bp目前没出大问题"
的数据支撑，用户据此判断"5bp暂时可用"，倾向于不再深挖这个问题，转向
确认参考基因组体系（0.5/7.2节）和ecotype-pca-panel分支的标签来源问题。
详见`docs/RESEARCH_ROADMAP.md`。

## 8. 后续三条工作线（讨论过，还没定最终方案）

这三条是用户提出的下一步方向，讨论了可行性和优先级，但都**还没有开始实际
写代码/跑分析**，需要先决定做哪个/什么顺序。

### 8.1（优先级已下调，见7.1）Oryza 参考基因组数据库整理/扩充

**这条原本排第一优先级，前提是"数据库野生稻缺失"——这个前提已经被7.1的
实测数据证伪，所以不再是紧急项**，但参考基因组体系本身"还是很乱"（用户
原话）这个观察依然成立，值得找机会理一遍，只是不再挡着best-hit这条线：

- `db/16/`（资源组A，NCBI datasets 16基因组+Liftoff注释）和 `db/3k/`（资源组B，
  3K Rice Genome Project数据）是**两个不同来源**，关系没有理清。
- `asn720data/`（720份现代/近现代品种PLINK面板）跟16个angkor古代样本的关系
  "尚未最终确认"（是否包含、是否capture panel来源）——⚠️2026-08-07更新：
  `ecotype-pca-panel`分支发现`asn720data/asn720.pop.fam`的FID列是`OrA-OrF`
  群体标签的关键来源，已不再是"可以忽略"的旧数据，详见
  `docs/RESEARCH_ROADMAP.md`。
- `db/3k/wild/` 这140+个野生稻/近缘种组装身份不明（见7.2，优先级已下调）。
- 已经有一个独立的 `results/07.wild_rice_alignment/`（minimap2 asm10 预设，
  野生稻组装比对到IRGSP，用于paftools.js call变异，标注"进行中"）——这是
  另一条已经在跑的、跟competitive mapping平行的分析线，不要跟best-hit阶段
  混淆，但如果这些wild genome的minimap2比对已经产出了变异，可能可以复用来
  确认这些野生稻组装的物种/谱系身份。

**现在真正优先级最高的是0.5节的taxid修正**，其次是7.5（损伤窗口校准，
已有初步数据支撑"5bp可用"）。这条可以往后放。

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

## 9. 给下一个接手者的具体待办（按优先级，2026-08-08更新后已重排）

1. 【当前最优先】0.5节：在服务器跑`fix_asian_rice_panel_taxid.sh`
   （先dry-run再--apply），修正np7的taxid错标，然后重新统计第6.1节的
   物种contig数表，更新这份文档。
2. 【次优先】7.5节：损伤窗口诊断——已有4个样本的初步数据支撑"5bp可用"，
   如果用户认为证据已经够用，可以不再深挖，直接按当前参数继续跑剩余样本；
   如果想更严谨验证，见7.5节的诊断脚本思路。
3. 【诊断，优先级降低，非阻塞】7.2（`db/3k/wild/`物种身份，与main分支
   `RESEARCH_ROADMAP.md`的P0第2条是同一个问题）、7.3（taxid rank校验）、
   7.4（MD tag覆盖率完整log、`seff`内存实测）——都可以做，但不再挡着
   best-hit往前走。
4. 只有在7.5的损伤窗口问题排查/调整完之后，才回头评估要不要给
   `oryza_besthit_damage_filter.py` 加一个margin参数（当前硬编码`<=`，没有
   margin概念）——在没排除损伤窗口这个更根本的原因之前，先调margin是在
   错的层面上打补丁。
5. 对剩余13个还在mapping队列里的样本，随mapping进度用 `submit_oryza_besthit.sh
   submit all` 陆续跑besthit（这条不依赖1-4，可以随时做）。
6. 8.1（参考基因组体系整理）——`asn720data`标签发现后优先级有所回升，
   但仍不是阻塞项。
7. 8.2（selection scan）和8.3（57基因SV判生态型）都还没真正开始。8.3的
   覆盖度QC是个独立、低成本、能快速出结论的子任务，可以和上面几条并行推进；
   8.2依赖更完整的genotype pipeline，且同样会受损伤/参考质量问题影响，
   优先级排在0.5/7.5之后。
