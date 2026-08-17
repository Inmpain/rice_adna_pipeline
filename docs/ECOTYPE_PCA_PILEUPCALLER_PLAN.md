# Ecotype PCA: pileupCaller + unified MAF + plot annotations — 执行计划（草案）

状态：**计划草案，只写文档、不写执行代码**。本文件是下一轮实现的蓝图，
等待用户逐项确认后再落到脚本。

分支：`codex/ecotype-pca-pileupcaller`（从 `codex/ecotype-pca-panel` @ `a99efd5`
分出，未改任何代码/脚本，只新增本计划文档）。

---

## 0. 本轮已确认的决定

- 三个现代面板（3K / 720 / Civán）**本轮都做**；服务器上已有可复用的中间文件就复用。
- 位点只做 **ALL track**，不再分 TV（本轮 TV 没有用）。
- 统一 **MAF = 0.01**；`geno` 仍按 `config/ecotype_pca_v2.yaml`：
  - 3K = 0.05
  - 720 = 0.10
  - Civán = 0.05
- 古样本 calling **确定替换为 sequenceTools `pileupCaller --randomHaploid`**，
  替换：
  - 共享轴：`scripts/ecotype_pca_v2/10_call_ancient_fixed_markers.py`
  - 私有轴：`scripts/ecotype_pca/pseudo_haploid_call.py`
- 覆盖位点的统计漏斗顺序固定为
  **参考基因组覆盖 → panel 交集 → MAF → LD**；“panel”这里指 **UNK 剔除个体后、
  还没做 MAF/LD 的原始 panel**（位点全集）。
- 绘图加两项标注：**覆盖位点数**、**PC 解释度**；**去掉红星 highlight**。
- 伪单倍体调用流程参考用户提供的 `Snakefile.pseudohaploid.from_panel`（本文件
  末尾附其关键参数映射）。

---

## 0.1 术语：marker 就是“位点”

本计划里的 **marker / marker 集 / fixed marker list** 和“**位点 / SNP 位点集**”是
同一个东西：一个 marker 就是一个基因组 SNP 位点（chr + pos + REF/ALT），marker 集
就是一批这样的位点的列表（仓库里叫 `*.fixed.snplist`）。

- “共享 marker 集” = 共享轴实际要用的那一批位点。
- “MAF ∩ coverage 的交集当共享 marker 集用” = 把“过了 MAF 且被古样本覆盖到的
  位点”这一批，作为共享轴要用的位点集。

所以它和“有覆盖的位置”不是两个东西：先有一堆候选位点，用 MAF 和 coverage 两个
条件筛一遍，筛剩下的那批位点就叫 marker 集。

---

## 1. 先厘清：共享轴 vs 私有轴

用户的问题“3K 的 coverage-first 是给私有轴降规模用；那共享轴的位点不是这个有
coverage 的部分吗”对应如下：

| | 共享轴（辅助分析） | 私有轴（主分析） |
|---|---|---|
| marker 集 | 现代面板 **MAF → ancient-union coverage 交集 → LD（覆盖候选内剪枝）** | **MAF → 每样本自己覆盖的位点**（3K 靠这步把 29M 降到可控） |
| 坐标 | 16 个古样本投影到**同一套**冻结坐标轴 | 每样本独立子集、独立坐标轴，**样本间不可比** |
| 出图 | 全部样品一张图（`26_plot_pc_pairs.sh`） | 每样品一格（`plot_pca_projection.py --evec-glob`） |
| coverage 交集 | **是**，MAF 之后交 ancient union（16 样本任一覆盖） | **是**，MAF 之后交每样本自己的覆盖位点 |

结论：两个轴都要做 coverage 交集，统计漏斗统一是
**参考基因组覆盖 → panel（UNK 剔除后、未做 MAF/LD）交集 → MAF → LD**。
区别只在“交谁的 coverage”：共享轴交 ancient union，私有轴交每样本自己的覆盖位点。

---

## 1.1 Phase 0（前置）— panel REF/ALT vs 参考基因组校验

在进入 Phase A 之前，先对每个面板做一次 REF/ALT 与参考基因组的一致性校验，避免把
方向错误的 REF/ALT 带进 pileupCaller 和下游 0/2 编码。

工具：`scripts/ecotype_pca_v2/23_validate_snp_ref_against_fasta.py`

```text
--snp <panel 原始 .snp>
--fasta /home/scratch/yinmt202607/db/asian_rice_panel_index/irgsp.fa
--out <panel>.ref_vs_fasta.report.tsv
```

三种结果：

- **PASS**：干净匹配 REF。
- **干净整体反标**：全位点一致地 REF/ALT 互换，且无 `true_mismatch` / 越界 →
  记档，决定是否在 pileupCaller `.snp` / `--a2-allele` 层面对齐。
- **FATAL**：散乱不一致 → 先停，查数据。

当前状态：

- **Civán**：51 runner 的 Step 0 已跑过，确认是**系统性整体反标**（REF↔ALT 互换）。
  这个反标对**旧 pysam 伪单倍体**是“内部自洽、无害”（panel 矩阵和古样本都读同一份
  `.snp` 列）；但对 **pileupCaller 不无害**——见下方 3.4。
- **3K / 720**：都还没跑过，需要新跑；720 的 REF/ALT 方向在
  `ECOTYPE_PCA_PANEL_QC_DESIGN.md` 里被标为“最急、最不了解”。

---

## 2. Phase A — marker 准备（漏斗：参考覆盖 → panel 交集 → MAF → LD，ALL track）

目标：先看古样本在参考基因组上覆盖了多少位点，再和 **UNK 剔除个体后、未做 MAF/LD
的 panel** 取交集，然后依次看交集位点里过 MAF、过 LD 的数量，最终得到一个可控、
且一定被古样本覆盖到的位点集合，供 `pileupCaller` / 私有轴 / 共享轴使用。

### 2.1 三个面板统一 MAF-only（复用现成脚本，不改代码）

工具：`scripts/ecotype_pca_v2/07_make_fixed_markers.sh`

参数：

```text
--stage geno_maf_only
--track ALL
--sensitivity primary
--library-type pooled_mixed
--config config/ecotype_pca_v2.yaml
--bfile <panel PLINK bfile>
--keep <reference keep-list>
--label {3k|720|civan}
```

输入：

- `<panel PLINK bfile>`：由 `02_convert_eigenstrat_for_plink.sh` 从 panel 的
  `.snp + .ind + .eigenstratgeno/.geno` 转出（PLINK bed/bim/fam）。
- `<reference keep-list>`：由 `06_build_reference_sample_set.py` 从真实 `.fam`
  生成。

输出（每面板一份）：

```text
{label}.pooled_mixed.ALL.primary.geno_maf_filtered.{bed,bim,fam}
{label}.pooled_mixed.ALL.primary.geno_maf_manifest.tsv
```

说明：manifest 是 TSV 汇总（不是 PLINK 格式），内含 `after_site_missingness`
和 `after_MAF` 两个独立计数。

### 2.2 coverage 交集（共享轴和私有轴都做）

coverage 普查（对**未处理 panel 的位点全集**做）：`19_survey_ancient_coverage.py
--panel-snp <该面板原始 .snp> --bam SAMPLE=BAM ×16`
→ `ancient_union_sites.tsv`（每个位点被哪些古样本覆盖，samples_covered 列）。

- **Civán**：这份 `ancient_union_sites.tsv` 已有（Stage 50 产出），可复用。
- **3K / 720**：需要各自跑一次（Civán 那份是针对 Civán 自己位点算的，不能直接套）。

交集后的两条路：

- **共享轴**：`21_extract_fixed_snplist.py` + `25_intersect_snplists.py` 交
  “MAF 过的位点 ∩ ancient union”，得到共享候选；
  然后再 `27_ancient_coverage_first_ld_prune.py` 在覆盖候选内做 LD（本轮确认做）。
- **私有轴**：交 “MAF 过的位点 ∩ 每样本自己的覆盖位点”，得到 per-sample 子集
  （不做 LD）；3K 正是靠这一步把 29M 降到可控，再往下 per-sample 子集筛选。

**coverage 一次计算、多处复用**：

- 每个 panel 只跑一次 `19_survey_ancient_coverage.py`，扫一遍 16 个 BAM。
- 之后共享轴（ancient union）、私有轴（每样本覆盖子集）、MAF 交集（25）、LD（27）、
  以及 2.3 的漏斗统计，**全部复用这份 coverage 结果，不再重扫 BAM**。
- 确认做：让 `19` **顺带写出每样本的覆盖位点清单**（不重扫 BAM，只是把 survey
  已经算出的 `sample_covered` 落盘），这样私有轴也直接吃 survey，而不是从 calling
  后的 `.call_sites.tsv` 反推；后续迭代只改清单交集，计算量小、也好定位错误。

最后用 `29_convert_plink_to_eigenstrat.sh` 把过滤后的 PLINK 转回 EIGENSTRAT，
供 `pileupCaller` / 私有轴（v1 `run_sample_panel_pca.sh` 只认 EIGENSTRAT）使用。

### 2.3 Phase A 统计漏斗（每面板一份 coverage 统计）

Phase A 结束时，每个面板产出一张覆盖度漏斗表（建议文件名
`*.coverage_funnel.tsv`），口径如下：

| # | 指标 | 定义 | 来源 |
|---|---|---|---|
| 1 | 参考基因组覆盖位点 | 古样本 read 在 IRGSP 参考基因组上覆盖到的位置数（全基因组 breadth） | `summarize_panel_overlap.py` 每样本 `genome_positions_low` / `genome_positions_high` |
| 2 | panel 覆盖位点 | 参考覆盖 ∩ **未处理 panel**（UNK 剔除后、未做 MAF/LD）的 SNP 位点 | `19_survey_ancient_coverage.py` → `ancient_union_sites.tsv` 行数（union）；per-sample 见其 `per_sample_coverage_summary.tsv` |
| 3 | panel 覆盖但 MAF 未过 | 第 2 项中 MAF < 0.01 的位点 | = (2) − (4) |
| 4 | panel 覆盖且 MAF 过 | 第 2 项 ∩ MAF-pass 位点 | `25_intersect_snplists.py`（MAF 过 ∩ ancient union） |
| 5 | panel 覆盖 + MAF 过 + LD 过 | 第 4 项中再经 LD 剪枝 | `27_ancient_coverage_first_ld_prune.py`（若本轮做 LD）或 07 `fixed` 产物 |

说明：

- 第 1 项“在参考基因组上覆盖的”是**全基因组覆盖 breadth**。主表建议用 **16 样本
  union（一个数字）**，这样和后面第 2–5 项（union 口径）能连成一个漏斗；需要看
  每个样本自己的覆盖时，再附一张 per-sample 明细表（`summarize_panel_overlap.py`
  每样本已能出 `genome_positions_low/high`，`19` 的 `per_sample_coverage_summary.tsv`
  则给出每样本 panel 覆盖数）。
- 第 2–5 项基本都来自现有脚本产出（19 / 07 / 25 / 27），Phase A 只需要一个小汇总脚本
  把它们拼成一张 `*.coverage_funnel.tsv`，不做新的重计算。
- 这张表同时能看到共享轴（union）和私有轴（per-sample）两个口径，但主表按
  “每面板一行”给；若需要 16 样本明细，直接用 `19` 的 `per_sample_coverage_summary.tsv`。

### 2.4 位点数 gate（不设人为下限）

Phase A 产出 marker 集后，**只有 0 位点的样本/面板才天然跳过**，其余都往下做
`pileupCaller` / pseudo-haploid / PCA：

- **共享轴**：`MAF ∩ coverage → LD` 之后 marker 数为 0 → 标记 `NO_MARKERS`，不做
  calling/PCA；否则照常跑。
- **私有轴**：某样本 per-sample 覆盖子集为 0 → 跳过该样本；否则照常跑。

位点多少只记录在统计表里，不据此丢弃。

---

## 3. Phase B — pileupCaller 替换（风险最大）

### 3.1 前置：安装 sequenceTools（pileupCaller）环境（服务器侧，用户执行）

新环境只装 `sequencetools`；`samtools` 走 `module load`，`plink` 用 base：

```bash
mamba create -n sequencetools
mamba activate sequencetools
mamba install bioconda::sequencetools

module load samtools
which pileupCaller samtools plink
pileupCaller -h
```

SLURM 作业里同样按这个顺序：先 `mamba activate sequencetools`，再
`module load samtools`，并确认 `plink` 仍能解析（如果 activate 后 `plink` 不见了，
就记 base 里 `plink` 的绝对路径，或改用 `module load plink`）。

### 3.2 新增两个调用脚本

- 共享轴：替换 `10_call_ancient_fixed_markers.py`
- 私有轴：替换 `pseudo_haploid_call.py`

两者核心同一套流水线（单样本 × 单 marker 集）：

```text
samtools mpileup \
  -R -B \
  -q {MAPQ} -Q {BASEQ} \
  -l {sites.bed} \
  -f {irgsp.fa} \
  {sample.besthit_oryza.irgsp.bam} \
| pileupCaller \
    --randomHaploid \
    --seed {stable_sample_seed} \
    --sampleNames {SAMPLE} \
    --samplePopName Rice \
    -f {marker.snp} \
    -p {prefix}

plink --bfile {prefix} \
  --keep-allele-order \
  --a2-allele {ref_alleles.a2} 2 1 \
  --make-bed --out {prefix}.a2
```

本项目的 ancient BAM 已经是单一 pooled BAM（`*.besthit_oryza.irgsp.bam`），
所以 Snakefile 里的“多 assay BAM merge”步骤**跳过**，直接对单 BAM 调用。

### 3.3 输出格式对接（关键）

`pileupCaller` 原生输出是 PLINK（bed/bim/fam）。下游现有脚本仍认
`.calls.txt` / `.call_sites.tsv` / `.call_report.tsv`：

- 共享轴下游：`11_build_ancient_callability.py`、`13_merge_ancients_fixed_panel.py`、
  `22_classify_scientific_projection.py`、`24_extract_sample_covered_sites.py`
- 私有轴下游：`build_sample_panel_subset.py`、`merge_ancient_into_panel.py`

推荐方案（改动最小）：写一层 **`pileupCaller PLINK → .calls.txt + .call_sites.tsv
+ .call_report.tsv` 转换**，保持 0/2/9、与 fixed `.snp` 行序一致、REF=2/ALT=0，
这样下游 11/13/22/24 和 v1 两个脚本**都不动**。备选是让下游直接吃 PLINK，但改动面大，先不选。

### 3.4 必须验证的点

- **REF/ALT 方向（关键）**：pileupCaller 用 `samtools mpileup -f irgsp.fa`，参考
  基因组的碱基是锚点；若把反标的 REF/ALT 直接当 pileupCaller 的 `.snp` 喂进去，
  0/2 会整体反掉，和现代 EIGENSTRAT 的 `REF=2; ALT=0; MISSING=9` 对不上。
  - 先按 Phase 0 的 `23_validate_snp_ref_against_fasta.py` 结果，把 pileupCaller
    的 `.snp` **纠正到 FASTA 方向**（REF = irgsp.fa 里的碱基；Civán 要 REF↔ALT 换回）。
  - 再用 Snakefile 的 `plink --a2-allele <ref> 2 1 --make-bed` 把输出对齐到同一方向。
  - 最强验证：用一个 panel 里已知群体的现代样本做 leave-one-out（`pseudo_haploid_call.py`
    文档里写的那套），确认它投影回自己已知群体，而不是往相反群体漂——这比只看 FASTA
    匹配率更直接。
- **seed 语义**：`pileupCaller --randomHaploid --seed` 是“每样本一个稳定 seed”，
  与现 pysam 版“每站点一个 seed”不同；本轮 ALL-only、不再要求 TV/ALL 同抽，
  所以这个差异可接受，但要在 `.call_report.tsv` 里记清 seed 合同。
- **spike 先行**：先 1 个样本 × 1 个面板（建议用 Civán 的 1015-marker 集）跑通
  `mpileup | pileupCaller → 转换 → merge → smartpca`，确认结果说得通，再铺开到
  16 样本 × 3 面板，不要一上来批量换。

---

## 4. Phase C — 绘图（加位点数 + PC 解释度，去红星）

出图规格：

- **共享轴**：`26_plot_pc_pairs.sh`，16 样品同一坐标系、全部样品一张图，
  PC1–PC10（5 对）。
- **私有轴**：PC1–PC10 都要，且分两种图——
  1. **16 样品组合 grid**：每个样品一格拼成一张大图
     （`plot_pca_projection.py --evec-glob`），至少 PC1–PC2；如需其它 PC 对
     同样按 5 对循环各出一张 grid。
  2. **每样品单独 PC1–PC10**：每个样品跑 5 对（1/2、3/4、5/6、7/8、9/10）。

两个脚本（`plot_pca_projection.py`、`26_plot_pc_pairs.sh`）都做以下改动：

1. **PC 解释度**：从 smartpca 的 `.eval` 文件读取 eigenvalue，
   `解释度 = eigen_i / sum(eigen)`，写在轴标签或子图标题里。
2. **覆盖位点数**：从 `.smartpca.log` 解析（或直接数私有轴 `*.subset.snp` /
   共享轴 `*.fixed_reference.snp` 的行数），标注在图内。
3. **去掉红星**：删除 `plot_pca_projection.py` 里
   `marker="*", color="red", s=170` 的古样本星形 highlight，改回普通散点
   （保留样本名 annotation 即可，不再画星）。
4. PC1–10 五对图照旧各跑 5 次（1/2、3/4、5/6、7/8、9/10）换文件名，不改画图逻辑。

---

## 5. 参数表

| 参数 | 值 | 来源 |
|---|---|---|
| MAF | 0.01（三面板统一） | config |
| geno | 3K=0.05 / 720=0.10 / Civán=0.05 | config |
| track | ALL（本轮不做 TV） | 本轮决定 |
| sensitivity | primary | 计划 |
| ancient MAPQ | 30 | config `ancient.mapq` |
| ancient BaseQ | 30 | config `ancient.baseq` |
| mpileup flags | `-R -B -q -Q -l -f` | Snakefile |
| pileupCaller | `--randomHaploid --seed --sampleNames --samplePopName` | Snakefile |
| REF/ALT 对齐 | `plink --a2-allele <ref> 2 1 --make-bed` | Snakefile |

**待确认**：Snakefile 默认 `MAPQ=25 / BASEQ=25`，本仓库冻结值是 `30/30`。
建议沿用仓库 `30/30`；若要改成 25 需用户明确。

---

## 6. 可复用中间文件 vs 需新生成

> 下面服务器路径按仓库记录（config + file_path.md + handoff）整理；
> **实际存在性需用户在服务器上用 `find`/`ls` 确认**，因为本机不直接连服务器。

### 6.1 直接复用（已在服务器，按仓库记录）

| 文件 | 路径 | 用途 |
|---|---|---|
| 3K 原始/过滤后 | `/home/scratch/yinmt202607/db/29M_3k/NB_final_snp.{snp,filtered.ind,filtered.eigenstratgeno}` | Phase A 输入 |
| 720 原始/过滤后 | `/home/scratch/yinmt202607/db/6.7M_720/asn720.6m.{snp,filtered.ind,filtered.geno}` | Phase A 输入 |
| Civán 原始/过滤后 | `/home/scratch/yinmt202607/db/paper1/civan_snp.{snp,filtered.ind,filtered.eigenstratgeno}` | Phase A 输入 |
| Civán coverage 普查 | `<results_v2_root>/**/ancient_union_sites.tsv` | Civán 私有/共享候选（Stage 50 产出） |
| Civán reference keep | `<results_v2_root>/**/*.reference_samples.keep`（595） | Civán 参考集 |
| Civán REF/ALT 校验 | `<results_v2_root>/**/civan.ref_vs_fasta.report.tsv` | 方向核对 |
| 16 ancient BAM | `/home/scratch/yinmt202607/gene/results/ecotype_pca/bam_irgsp/*.besthit_oryza.irgsp.bam` | pileupCaller 输入 |
| IRGSP 参考 FASTA | `/home/scratch/yinmt202607/db/asian_rice_panel_index/irgsp.fa` | mpileup `-f` |
| v1 16×3 first-look | `/home/scratch/yinmt202607/gene/results/ecotype_pca/pca_runs/` | 仅作 before 对照，不直接复用为新结果 |

### 6.2 需要新生成

| 文件 | 工具 | 说明 |
|---|---|---|
| 三面板 `geno_maf_filtered` | 07 `--stage geno_maf_only` | Civán 的 ALL 份可能已由 51 跑过，需确认后复用 |
| 3K/720 coverage 普查 | 19 `--panel-snp <panel>.snp` | 新（Civán 复用 Stage 50 的） |
| 3K/720 coverage 交集 snplist | 25 交集 | 新（共享轴交 union；私有轴交 per-sample） |
| 3K/720 REF/ALT 校验 | 23 `--snp <panel>.snp --fasta irgsp.fa` | 新（Civán 复用 51 的） |
| 私有轴 EIGENSTRAT | 29 | 新 |
| pileupCaller 转换层输出 `.calls.txt/.call_sites.tsv/.call_report.tsv` | 新脚本 | Phase B |
| 改版后的 PC 图 | plot 脚本 | Phase C |

---

## 7. 建议执行顺序

`Phase 0 → A → C → B`，并对 B 做一处解耦：

- **Phase 0**：REF/ALT 校验，先确认方向，再进 A。
- **A**：机械、低风险、全复用现成脚本，且产出是 B 的输入。
- **C**：改动范围小、独立，先拿下。
- **B**：风险最大，但 **pileupCaller 的 spike 可以提前**——不依赖三面板 MAF，
  只要 `which pileupCaller` 有结果，就能拿 Civán 现成 1015-marker 集 + 1 个样本
  先验证工具/输出格式，把最大未知项最早暴露。

---

## 8. 服务器侧待确认 + 待拍板

1. `which pileupCaller` / `pileupCaller -h` 是否有输出（没装则先装 sequenceTools）。
2. 上面 6.1 列出的可复用中间文件，哪些在服务器上**真实存在**、实际路径是否一致。
3. MAPQ/BaseQ 用仓库冻结 `30/30`，还是改 Snakefile 的 `25/25`。
4. Phase B 输出对接：确认“写一层 PLINK→calls 转换、下游不动”这个方案。
5. Civán 的 `geno_maf_filtered`（ALL）如果 51 已产出，是否直接复用。
6. 共享轴 LD 已确认本轮做（`27`）；位点数下限**不设人为阈值**，有位点的都跑。

---

## 附：Snakefile 关键参数 → 本项目映射

| Snakefile 项 | 本项目对应 |
|---|---|
| `samtools mpileup -R -B -q -Q -l -f` | 保留，直接采用 |
| `pileupCaller --randomHaploid --seed --sampleNames --samplePopName -f -p` | 保留 |
| `plink --a2-allele ... 2 1 --make-bed` | 保留，用于 REF/ALT 对齐 |
| 多 assay BAM merge | 本项目跳过（单 pooled BAM） |
| `MAPQ=25 / BASEQ=25` | 建议改仓库冻结 `30/30`（待确认） |
| 输出 PLINK bed/bim/fam | 需转回 `.calls.txt/.call_sites.tsv/.call_report.tsv` 兼容下游 |
