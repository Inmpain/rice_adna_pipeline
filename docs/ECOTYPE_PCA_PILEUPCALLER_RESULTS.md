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

## 3. 720 coverage 漏斗

```text
6,769,714 raw
  → geno 0.20         ~1.26M
  → MAF 0.01          1,208,247
  → ∩ ancient coverage  834
  → r² LD(100kb/0.2)     758
```

结论：古代覆盖度是真正瓶颈（16 个古 BAM 只覆盖 720 的 5,192 个位点），r² LD 本身
只温和减 76 个；古代覆盖位点成簇（872 个 MAF∩coverage 里一半相距 <5kb），故古代
投影保留 r² LD、不用 `--bp-space 5000` 物理抽稀。

---

## 4. 共享 marker 数

```text
3K = 2859
720 = 758
Civán = 1015
```

---

## 5. pileupCaller 环境

| 版本 | 路径 | 状态 |
|---|---|---|
| v1.5.3.1 | `~/software/pileupCaller-linux` | ✅ 可用 |
| v1.6.0.0（Bioconda） | — | ❌ segfault（需 GLIBC_2.34，服务器 glibc 更老） |

配套：`plink` 在 `snakemake` 环境（`mamba install bioconda::plink`），`samtools`
走 `module load samtools`。

---

## 6. smartpca 投影（失败，本轮不追）

smartpca `lsqproject` 把古代样本判为 `insufficient data` 并从 `.evec` 剔除，触发
`15_pca_qc.py` FATAL。根因与早期 Civán 共享轴相同：古代覆盖位点与 marker 集交集过小。
本轮不再追，详见 `ECOTYPE_PCA_PILEUPCALLER_PATH_MAP.md` 第 6 节。
