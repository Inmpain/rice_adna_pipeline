# Ecotype PCA pileupCaller — Phase A 进度与结果

状态：**Phase A 进行中**（2026-08-18）。本文件记录“转 PLINK + 锁 REF + MAF”的
实际数字、720 的缺失率审计，以及据此做出的参数调整。

---

## 1. 三个面板 Phase A 进度

| 面板 | 02 转换 | 锁 A2=irgsp | 06 keep-list | 07 MAF(geno/maf) | 状态 |
|---|---|---|---|---|---|
| 3K | 29,635,224 × 3000 | 29,635,224 assignments | 2733 参考样本 | 16,825,293 → 4,592,439 | 完成 |
| 720 | 6,769,714 × 718 | 6,769,700 assignments（14 skip） | 718（all_modern，待改） | 51,549 → 45,378 | 待按新参数重跑 |
| Civán | 复用已有 | 待补 | 复用 595 | 复用 1015 marker | 待补归一化 |

### 1.1 3K 参考样本构成（06）

- reference/axis-builder N = 2733：
  `TEJ=319, IND=1743, TRJ=388, AUS=215, ARO=68`
- other（不建轴）：`JAPONICA_UNSPEC=132, INTERMEDIATE_TYPE=135`

### 1.2 720 参考样本构成（06，改动前）

- `axis_mode: all_modern` → 718 全建轴：
  `OrD=66, OrA=133, OrF=34, OrC=38, OrADM=41, OrE=24, OrB=64, RAY=9,
  TRJ=48, ARO=15, IND=148, JAPONICA_UNSPEC=18, INTERMEDIATE_TYPE=12,
  TEJ=30, AUS=38`
- WARNING：未用 `--smiss-file` 检查“全缺失样本”。

---

## 2. 720 缺失率审计（03）

位点缺失率分布：

| 缺失率区间 | 位点数 | 占比 |
|---|---|---|
| 0–0.01 | 132 | ~0% |
| 0.01–0.05 | 9,920 | 0.15% |
| 0.05–0.10 | 41,497 | 0.61% |
| 0.10–0.20 | 1,210,029 | 17.9% |
| >0.20 | 5,508,136 | 81.4% |

样本缺失率：野生稻（ERR/SRR）约 50–56%（例 `ERR068594=0.539`）。属**结构性高缺失**
（野生高、栽培低）。

结论：`geno 0.10` 只保留 0.76%，太严；**放宽为 `geno 0.20`**（保留 ~1.26M，18.6%）。

---

## 3. 本轮参数调整（已落入 config）

1. `panel_B_720.geno: 0.10 → 0.20`（已写入 `scripts/ecotype_pca_v2/config/ecotype_pca_v2.yaml`；此前文档已定、config 仍 0.10 是一处 drift）。
2. `panel_B_720.axis_mode`：**保持 `all_modern`**（野生稻参与建轴），不改成栽培锚点。

理由：720 来源论文用 `ngsCovar` 对基因型似然把栽培+野生一起建轴；先沿用这一思路，
用三档 PCA 图验证“硬基因型 + geno 0.20 + LD”下的结构再决定是否换轴。

---

## 4. 中间文件（Phase A）

| 文件 | 路径 |
|---|---|
| 3K 锁 REF 后 PLINK | `<results_v2_root>/phaseA/3k/plink/NB_final_snp.irgsp.{bed,bim,fam}` |
| 3K MAF 后 PLINK | `<results_v2_root>/phaseA/3k/maf_ld/3k.pooled_mixed.ALL.primary.geno_maf_filtered.*` |
| 3K ref 清单 | `<results_v2_root>/phaseA/3k/ref/NB_final_snp.irgsp_ref.txt` |
| 3K keep-list | `<results_v2_root>/phaseA/3k/reference_sets/3k.reference_samples.keep` |
| 720 锁 REF 后 PLINK | `<results_v2_root>/phaseA/720/plink/asn720.6m.irgsp.{bed,bim,fam}` |
| 720 ref 清单 | `<results_v2_root>/phaseA/720/ref/asn720.6m.irgsp_ref.txt` |
| 720 缺失率审计 | `<results_v2_root>/phaseA/720/audit/720.audit.{missingness,samples}.tsv` |

> `<results_v2_root>` = `/home/scratch/yinmt202607/gene/results/ecotype_pca_v2`

---

## 5. 下一步

- 720 重跑 07（geno 0.20 + 保持 `all_modern` 建轴；不要写成栽培锚点）。
- 三面板各跑 coverage 普查（19）→ MAF∩coverage（25）→ 共享轴 LD（27）。
- 三档 modern-only PCA 诊断图（raw / MAF-geno / LD），每档 PC1–PC10 + 解释度 + 位点数。
