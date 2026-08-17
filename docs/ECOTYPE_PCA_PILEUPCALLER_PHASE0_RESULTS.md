# Ecotype PCA pileupCaller — Phase 0 结果与中间文件

状态：**Phase 0 完成**（2026-08-17）。本文件记录“REF/ALT 方向体检”的结论、
产出的中间文件路径，以及这些结果后续怎么用。

---

## 1. 结论摘要

| Panel | pattern | 关键数字 |
|---|---|---|
| Civán | `systematic_ref_alt_swap` | 整体反标 |
| 3K | `systematic_ref_alt_swap` | 29,635,224 / 29,635,224 全部 match ALT |
| 720 | `inconsistent_requires_manual_review` | match_ref=6,164,430；match_alt=605,270；mismatch=14 |

处理方式：

- **3K / Civán**：干净整体反标 → 全局处理。
- **720**：混合方向 → 605,270 个位点逐位点翻 + 14 个 mismatch 单独标记缺失/剔除。
- **归一化机制**：给每 panel 生成 `SNP_ID + irgsp_base` 参考等位基因清单，用
  `plink --a2-allele <ref_list> 2 1 --keep-allele-order --make-bed` 在**每次写新 bed
  后**锁 A2 = irgsp 参考碱基（plink2 步骤用 `--ref-allele` 语义）。

---

## 2. 服务器中间文件清单

| 文件 | 路径 | 内容 |
|---|---|---|
| 3K 方向报告 | `/home/scratch/yinmt202607/gene/results/ecotype_pca_v2/phase0/3k.ref_vs_fasta.report.tsv` | 3K 整体反标，match_alt=29,635,224 |
| 720 方向报告 | `/home/scratch/yinmt202607/gene/results/ecotype_pca_v2/phase0/720.ref_vs_fasta.report.tsv` | 720 混合 + 14 个 mismatch 示例 |
| 720 翻链清单 | `/home/scratch/yinmt202607/gene/results/ecotype_pca_v2/phase0/720.flip.snplist` | 605,270 个 SNP ID（match_alt，需翻） |
| 720 mismatch 明细 | `/home/scratch/yinmt202607/gene/results/ecotype_pca_v2/phase0/720.mismatch.tsv` | 14 行：snp_id/chrom/pos/ref/alt/fasta_base |
| Civán 方向报告 | `<results_v2_root>/**/civan.ref_vs_fasta.report.tsv` | 51 runner Step 0 产出，整体反标 |

> 路径按仓库记录整理；`<results_v2_root>` =
> `/home/scratch/yinmt202607/gene/results/ecotype_pca_v2`。

---

## 3. 720 混合方向的详细分布

### 3.1 按染色体

| chrom | match_ref | match_alt | mismatch | alt_frac |
|---|---|---|---|---|
| chr01 | 701,676 | 74,099 | 0 | 0.0955 |
| chr02 | 612,410 | 61,989 | 0 | 0.0919 |
| chr03 | 621,820 | 65,849 | 0 | 0.0958 |
| chr04 | 539,602 | 49,207 | 10 | 0.0836 |
| chr05 | 534,625 | 47,026 | 0 | 0.0808 |
| chr06 | 534,481 | 52,823 | 0 | 0.0899 |
| chr07 | 469,675 | 46,937 | 0 | 0.0909 |
| chr08 | 500,158 | 44,938 | 0 | 0.0824 |
| chr09 | 375,484 | 35,072 | 0 | 0.0854 |
| chr10 | 392,969 | 43,886 | 2 | 0.1005 |
| chr11 | 477,594 | 45,606 | 2 | 0.0872 |
| chr12 | 403,936 | 37,838 | 0 | 0.0857 |

### 3.2 按等位基因对

| pair | match_ref | match_alt | mismatch | alt_frac |
|---|---|---|---|---|
| C/T | 2,235,459 | 223,873 | 9 | 0.0910 |
| A/G | 2,234,828 | 222,996 | 1 | 0.0907 |
| A/T | 526,276 | 45,812 | 2 | 0.0801 |
| A/C | 449,841 | 43,321 | 0 | 0.0878 |
| G/T | 449,259 | 42,896 | 1 | 0.0872 |
| C/G | 268,767 | 26,372 | 1 | 0.0894 |

结论：`match_alt` 占比在所有染色体和所有等位基因对上都在 8–10%，**无染色体/碱基
类型偏好，是均匀混合**。最可能是 720 面板的 REF/ALT 不是严格按 irgsp 参考定义
（例如按 major/minor 频率、或从其他参考 liftOver 而来）。

### 3.3 14 个 mismatch 明细

| SNP | chr:pos | declared REF/ALT | FASTA |
|---|---|---|---|
| 4np5072535 | chr04:5072535 | T/C | A |
| 4np5422723 | chr04:5422723 | T/C | A |
| 4np5422731 | chr04:5422731 | C/T | A |
| 4np19053798 | chr04:19053798 | C/T | A |
| 4np20939858 | chr04:20939858 | T/A | N |
| 4np22899953 | chr04:22899953 | T/C | G |
| 4np28353118 | chr04:28353118 | C/T | G |
| 4np28845702 | chr04:28845702 | A/T | G |
| 4np31370380 | chr04:31370380 | G/C | A |
| 4np32261676 | chr04:32261676 | G/A | C |
| 10np7475160 | chr10:7475160 | C/T | A |
| 10np21601148 | chr10:21601148 | C/T | G |
| 11np16411945 | chr11:16411945 | T/G | A |
| 11np20582932 | chr11:20582932 | C/T | G |

规律：10 个在 chr04、2 个 chr10、2 个 chr11；FASTA 上几乎都是第三个碱基（A/G/C），
1 个是 N（`chr04:20939858`）。属多等位/参考组装差异/低质量区，不是方向反。

---

## 4. 下游怎么用

- `720.flip.snplist` / `720.mismatch.tsv` 是**诊断记录**，用于核对和留档。
- 真正的方向归一化不直接靠翻链清单，而是：
  1. 生成 `SNP_ID + irgsp_base` 参考等位基因清单（`make_irgsp_ref_list.py`，Phase A 待写）；
  2. 在每次 `--make-bed` 后用 `--a2-allele` / `--ref-allele` 锁 A2 = irgsp；
  3. 用 leave-one-out 现代样本验证投影方向。
- 14 个 mismatch 位点在 marker 集里标记缺失或直接剔除。
