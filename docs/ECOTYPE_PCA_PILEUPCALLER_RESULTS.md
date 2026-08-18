# Ecotype PCA v2（pileupCaller 版）— 结果汇总

状态：截至 **2026-08-18** 的实际产出数字汇总。本文是结果入口；逐项明细在
`ECOTYPE_PCA_PILEUPCALLER_PHASE0_RESULTS.md`、`ECOTYPE_PCA_PILEUPCALLER_PHASEA_RESULTS.md`。

---

## 1. Phase 0：REF/ALT 方向体检

| Panel | pattern | 关键数字 |
|---|---|---|
| 3K | `systematic_ref_alt_swap` | 29,635,224 / 29,635,224 全部 match ALT |
| Civán | `systematic_ref_alt_swap` | 整体反标 |
| 720 | `inconsistent_requires_manual_review` | match_ref=6,164,430；match_alt=605,270；mismatch=14 |

720 混合方向无染色体/碱基偏好（各染色体 match_alt 占比 8–10%）；14 个 mismatch 中
10 个在 chr04、2 个 chr10、2 个 chr11，FASTA 上多为第三碱基（1 个 N）。

---

## 2. Phase A：转 PLINK + 锁 A2=irgsp + MAF/geno

| Panel | 转换规模 | 锁 A2=irgsp | ref keep-list | MAF（geno/maf） |
|---|---|---|---|---|
| 3K | 29,635,224 × 3000 | 29,635,224 assignments | 2733 参考样本 | 16,825,293 → 4,592,439 |
| 720 | 6,769,714 × 718 | 6,769,700 assignments（14 skip） | 718（all_modern） | 51,549 → 45,378（旧 geno 0.10） |
| Civán | 复用已有 | 待补 | 复用 595 | 复用 1015 marker |

720 缺失率审计结论：`geno 0.10` 只保留 0.76%，过严；**放宽 `geno 0.20`**（保留
~1.26M，18.6%）。

---

## 3. coverage 漏斗

### 720

```text
6,769,714 raw
  → geno 0.20         ~1.26M
  → MAF 0.01          1,208,247
  → ∩ ancient coverage  872
  → r² LD(100kb/0.2)     795
```

结论：古代覆盖度是真正瓶颈（16 个古 BAM 只覆盖 720 的 5,192 个位点），r² LD 本身
只温和减 77 个；古代覆盖位点成簇（872 个 MAF∩coverage 里一半相距 <5kb），故古代
投影保留 r² LD、不用 `--bp-space 5000` 物理抽稀。

### 3K（coverage_funnel，union 口径）

```text
raw_panel_covered        38,665
  → maf_passed_covered    4,415   ← MAF∩coverage（覆盖候选，未 LD）
  → ld_passed_covered     2,859   ← coverage-first r² LD 后
```

即 3K 的覆盖候选 = **4,415**，LD 后 = **2,859**（删 1,556）。与 720 的 872→795 相比，
3K 覆盖候选多得多，但 LD 剪枝也重得多（高密度 29M 面板，覆盖位点彼此 LD 密集）。

---

## 4. 共享 marker 数

```text
3K = 2859（覆盖候选 4415，LD 后）
720 = 795（覆盖候选 872，LD 后）
Civán = 1015
```

> 注意（2026-08-18 更正）：上面这行是「coverage-first」路线的共享 marker 数。
> 已跑通并验证的 720 共享投影（`720hybrid.v2.final.png`，20:14）用的是 **44,920
> marker** 的 pileupCaller 位点集（pileupCaller TotalSites=44920），不是 795——两者
> 是不同配方（44,920 的「5kb backbone + ancient covered」合成集 vs 795 的
> coverage-first LD）。3K 用哪个配方尚未定，见 PLAN「实际执行策略」。

---

## 5. pileupCaller 环境

| 版本 | 路径 | 状态 |
|---|---|---|
| v1.5.3.1 | `~/software/pileupCaller-linux` | ✅ 可用 |
| v1.6.0.0（Bioconda） | — | ❌ segfault（需 GLIBC_2.34，服务器 glibc 更老） |

配套：`plink` 在 `snakemake` 环境（`mamba install bioconda::plink`），`samtools`
走 `module load samtools`。

---

## 6. smartpca 投影（早期失败，已被 pileupCaller 版 supersede）

smartpca `lsqproject` 把古代样本判为 `insufficient data` 并从 `.evec` 剔除，触发
`15_pca_qc.py` FATAL。根因与早期 Civán 共享轴相同：古代覆盖位点与 marker 集交集过小。

**superseded（2026-08-18）**：改用 pileupCaller 的 44,920-marker 共享投影已跑通，
16 古样本全部进 `.evec`（call 21–1314）。本条只留作历史，别再当当前状态引用；
当前待办是「两张图 PC 解释度差异」的根因核对，见 `HANDOFF` 第 4 节 B7。
