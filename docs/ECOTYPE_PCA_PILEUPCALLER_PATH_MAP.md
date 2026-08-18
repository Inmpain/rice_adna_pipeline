# Ecotype PCA v2（pileupCaller 版）— 路径与状态交接文档

> 本文是 `codex/ecotype-pca-pileupcaller` 分支当前工作状态的**路径地图**，
> 按「仓库内脚本 → 仓库内文档 → 服务器输入 → 服务器中间/结果文件」四层记录。
> 新接手的人先读本文件，再按需翻 `ECOTYPE_PCA_PILEUPCALLER_PLAN.md`（计划）、
> `ECOTYPE_PCA_PILEUPCALLER_PHASE0_RESULTS.md` / `PHASEA_RESULTS.md`（各阶段结果）、
> `decisions_log.md`（决策原因）和 `file_path.md`（全项目 v5 路径总图）。

最后更新：2026-08-18

---

## 0. 一句话现状

Phase 0（REF/ALT 方向体检）与 Phase A（转 PLINK + 锁 A2=irgsp + MAF/geno）已跑通并
留档；pileupCaller 替换调用（`mpileup | pileupCaller --randomHaploid` → PLINK →
`.calls.txt`）与配套绘图脚本已落地到仓库。**末尾 smartpca `lsqproject` 投影失败
（insufficient data / 古代样本被剔出 `.evec`）已知、本轮不追**，只记为待办，
不影响本交接文档。详见第 6 节。

---

## 1. 仓库内脚本路径清单（`rice_adna_pipeline/` 下）

### 1.1 新脚本（本轮 pileupCaller 流水线，位于 `scripts/ecotype_pca_v2/`）

| 脚本 | 作用 |
|---|---|
| `scripts/ecotype_pca_v2/23_validate_snp_ref_against_fasta.py` | 校验 panel `.snp` 的 REF/ALT 与 `irgsp.fa` 方向；区分“干净整体反标”与“散乱不一致/真 mismatch” |
| `scripts/ecotype_pca_v2/extract_mismatches.py` | 把 mismatch 位点（FASTA 第三碱基/N）落成 TSV |
| `scripts/ecotype_pca_v2/make_flip_list.py` | 生成需翻链的 SNP ID 清单（720 混合方向用） |
| `scripts/ecotype_pca_v2/make_irgsp_ref_list.py` | 生成 `SNP_ID<TAB>irgsp_base` 参考等位基因清单，供 `plink --a2-allele` 锁方向 |
| `scripts/ecotype_pca_v2/tally_720_refalt.py` | 按染色体/等位基因对统计 720 的 match_ref/match_alt/mismatch |
| `scripts/ecotype_pca_v2/tally_720_mismatch_by_pop.py` | 对 14 个 mismatch 位点按群体统计 720 基因型 |
| `scripts/ecotype_pca_v2/27_ancient_coverage_first_ld_prune.py` | 在 ancient 覆盖候选内做 LD 剪枝（coverage-first） |
| `scripts/ecotype_pca_v2/28_maf_filter_eigenstrat_for_private_axis.py` | 私有轴 EIGENSTRAT MAF 过滤 |
| `scripts/ecotype_pca_v2/29_convert_plink_to_eigenstrat.sh` | PLINK → EIGENSTRAT（给私有轴 / pileupCaller） |
| `scripts/ecotype_pca_v2/coverage_funnel.py` | 逐样本 coverage 漏斗（raw → MAF → LD） |
| `scripts/ecotype_pca_v2/plot_panel_pca.py` | 三档 modern-only PCA 诊断图 |
| `scripts/ecotype_pca_v2/pileupcaller_shared_call.sh` | 单样本 `mpileup | pileupCaller --randomHaploid` 调用（共享轴） |
| `scripts/ecotype_pca_v2/pileupcaller_plink_to_calls.py` | pileupCaller 的 PLINK 输出 → `.calls.txt`（0/2/9，`.bim` 行序） |
| `scripts/ecotype_pca_v2/plot_smartpca_evec.py` | 画 smartpca `.evec/.eval`：群体着色、Ancient 黑点、PC 解释度、marker 数、无红星 |

### 1.2 workflow runners（`scripts/ecotype_pca_v2/workflow/runners/`）

| runner | 作用 |
|---|---|
| `51_civan_maf_ld_and_private_axis.sh` | Civán：MAF + LD + 私有轴（已按 coverage-first 顺序改造） |
| `61_panel_maf_shared_projection.sh` | 3K/720 的 MAF-only 共享投影通用 runner（`--panel A|B`） |

### 1.3 配置与冻结值

| 文件 | 内容 |
|---|---|
| `scripts/ecotype_pca_v2/config/ecotype_pca_v2.yaml` | 所有数值参数唯一来源（MAF/geno/LD/MAPQ/BaseQ/服务器输入路径等） |
| `scripts/ecotype_pca_v2/workflow/workflow.json` | 阶段顺序、tracked_files、manual gate |

关键冻结值（见 config）：

```text
MAF         = 0.01（三面板统一）
geno        = 3K 0.05 / 720 0.20 / Civán 0.05
track       = ALL（本轮不做 TV）
LD          = window 100kb / r2 0.20
ancient MAPQ= 30, BaseQ = 30
smartpca    = lsqproject: true, numoutlieriter: 0, num_pcs: 10
```

---

## 2. 仓库内文档路径（`docs/`）

| 文档 | 说明 |
|---|---|
| `docs/ECOTYPE_PCA_PILEUPCALLER_PLAN.md` | pileupCaller 替换 + 统一 MAF + 绘图标注的完整计划（554 行） |
| `docs/ECOTYPE_PCA_PILEUPCALLER_PHASE0_RESULTS.md` | Phase 0 方向体检结果与中间文件 |
| `docs/ECOTYPE_PCA_PILEUPCALLER_PHASEA_RESULTS.md` | Phase A 转 PLINK/锁 REF/MAF 数字与参数调整 |
| `docs/ECOTYPE_PCA_PILEUPCALLER_PATH_MAP.md` | 本文件（路径地图/交接） |
| `docs/decisions_log.md` | 关键决策（含 Civán 顺序纠正、720 建轴保持 all_modern 等） |
| `docs/file_path.md` | 全项目 v5 路径总图（更老但更全） |
| `docs/ECOTYPE_PCA_V2_SPEC.md` | v2 冻结规格 |

---

## 3. 服务器输入数据路径

| 数据 | 路径 |
|---|---|
| 3K 原始/过滤后 | `/home/scratch/yinmt202607/db/29M_3k/NB_final_snp.{snp,filtered.ind,filtered.eigenstratgeno}` |
| 720 原始/过滤后 | `/home/scratch/yinmt202607/db/6.7M_720/asn720.6m.{snp,filtered.ind,filtered.geno}`（注意 `.geno` 不是 `.eigenstratgeno`） |
| Civán 原始/过滤后 | `/home/scratch/yinmt202607/db/paper1/civan_snp.{snp,filtered.ind,filtered.eigenstratgeno}` |
| IRGSP 参考 FASTA | `/home/scratch/yinmt202607/db/asian_rice_panel_index/irgsp.fa` |
| 16 古样本 BAM | `/home/scratch/yinmt202607/gene/results/ecotype_pca/bam_irgsp/*.besthit_oryza.irgsp.bam` |
| pileupCaller 二进制 | `~/software/pileupCaller-linux`（v1.5.3.1，`PILEUP_CALLER` 环境变量可覆盖） |

`results_v2_root`：

```text
/home/scratch/yinmt202607/gene/results/ecotype_pca_v2
```

---

## 4. 服务器中间/结果文件路径

### 4.1 Phase 0（REF/ALT 方向体检）

| 文件 | 路径 |
|---|---|
| 3K 方向报告 | `/home/scratch/yinmt202607/gene/results/ecotype_pca_v2/phase0/3k.ref_vs_fasta.report.tsv` |
| 720 方向报告 | `/home/scratch/yinmt202607/gene/results/ecotype_pca_v2/phase0/720.ref_vs_fasta.report.tsv` |
| 720 翻链清单 | `/home/scratch/yinmt202607/gene/results/ecotype_pca_v2/phase0/720.flip.snplist` |
| 720 mismatch 明细 | `/home/scratch/yinmt202607/gene/results/ecotype_pca_v2/phase0/720.mismatch.tsv` |
| Civán 方向报告 | `<results_v2_root>/**/civan.ref_vs_fasta.report.tsv`（51 runner Step 0 产出） |

Phase 0 结论：

```text
3K    : systematic_ref_alt_swap（29,635,224 / 29,635,224 全部 match ALT）
Civán : systematic_ref_alt_swap（整体反标）
720   : inconsistent_requires_manual_review（match_ref=6,164,430 / match_alt=605,270 / mismatch=14）
```

### 4.2 Phase A（转 PLINK + 锁 A2=irgsp + MAF/geno）

| 文件 | 路径 |
|---|---|
| 3K 锁 REF 后 PLINK | `<results_v2_root>/phaseA/3k/plink/NB_final_snp.irgsp.{bed,bim,fam}` |
| 3K MAF 后 PLINK | `<results_v2_root>/phaseA/3k/maf_ld/3k.pooled_mixed.ALL.primary.geno_maf_filtered.*` |
| 3K ref 清单 | `<results_v2_root>/phaseA/3k/ref/NB_final_snp.irgsp_ref.txt` |
| 3K keep-list | `<results_v2_root>/phaseA/3k/reference_sets/3k.reference_samples.keep` |
| 720 锁 REF 后 PLINK | `<results_v2_root>/phaseA/720/plink/asn720.6m.irgsp.{bed,bim,fam}` |
| 720 ref 清单 | `<results_v2_root>/phaseA/720/ref/asn720.6m.irgsp_ref.txt` |
| 720 缺失率审计 | `<results_v2_root>/phaseA/720/audit/720.audit.{missingness,samples}.tsv` |

共享 marker 数（7.1.1）：

```text
3K = 2859，720 = 758，Civán = 1015
```

720 漏斗：`6,769,714 raw → geno0.20 ~1.26M → MAF0.01 1,208,247 → ∩ ancient coverage 834 → r² LD 758`。

### 4.3 其他可复用中间文件

| 文件 | 路径 |
|---|---|
| Civán coverage 普查 | `<results_v2_root>/**/ancient_union_sites.tsv`（Stage 50 产出） |
| Civán reference keep | `<results_v2_root>/**/*.reference_samples.keep`（595） |
| v1 16×3 first-look | `/home/scratch/yinmt202607/gene/results/ecotype_pca/pca_runs/`（仅作 before 对照） |

---

## 5. 阶段状态

| 阶段 | 状态 | 备注 |
|---|---|---|
| Phase 0 方向体检 | 完成 | 3K/Civán 整体反标，720 混合（605,270 翻链 + 14 mismatch） |
| Phase A marker 准备 | 3K 完成，720 待按 `geno 0.20` 重跑，Civán 复用 | 见 `PHASEA_RESULTS.md` |
| Phase B pileupCaller 替换 | 脚本已落地（spike 工具可用），批量铺开待做 | 见第 1.1 节新脚本 |
| Phase C 绘图 | 绘图脚本已落地 | 三档 modern-only + smartpca `.evec` |
| smartpca 投影 | 末尾失败，本轮不追 | 见第 6 节 |

---

## 6. 投影失败说明（不追，仅留档）

末尾 smartpca `lsqproject` 投影失败，表现为古代样本被 `insufficient data` 剔除出
`.evec`（`15_pca_qc.py` 的 FATAL 门槛会拦下）。根因与此前 Civán 共享轴一致：
古代覆盖位点稀疏、与 marker 集交集过小，导致可调用位点数不足。

本轮决定：**不继续追这条投影失败**。它不影响已完成的 Phase 0 / Phase A 留档，
也不影响仓库内 pileupCaller 脚本与绘图脚本的交付。后续若要恢复，优先排查：

1. 各样本 `coverage_funnel.tsv` 里 `ld_passed_covered` 是否为 0/极低；
2. 16 个古样本 BAM 与 marker `.snp` 的 contig 命名/排序是否一致（`chr01` vs 数字）；
3. `pileupcaller_shared_call.sh` 生成的 `.snp` 与 `.sites.bed` 是否与 mpileup `-l` 对齐；
4. 现代样本是否被错误放进 poplist 走投影路线（应为 Ancient 标签才允许缺失）。

---

## 7. 交接前必读（简要）

- 不要在任何步骤里硬编码数值参数，一律读 `config/ecotype_pca_v2.yaml`。
- 每次写新 PLINK bed 后都要重新 `--a2-allele <ref_list> 2 1 --keep-allele-order` 锁方向
  （尤其 07 的 plink2 MAF 之后）。
- 服务器路径以 `find`/`ls` 实际为准；本机不直连服务器，`<results_v2_root>` 一律指
  `/home/scratch/yinmt202607/gene/results/ecotype_pca_v2`。
- pileupCaller 二进制用 `~/software/pileupCaller-linux`（v1.5.3.1）；Bioconda 1.6.0.0
  在服务器 segfault，别用。
