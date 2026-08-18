# Ecotype PCA v2（pileupCaller 版）— 进展

状态：**进行中**（2026-08-18）。本文件记录当前进展：已完成什么、进行中什么、
阻塞什么、下一步做什么。数字口径见 `ECOTYPE_PCA_PILEUPCALLER_RESULTS.md`，
完整路径见 `ECOTYPE_PCA_PILEUPCALLER_PATH_MAP.md`，执行策略见
`ECOTYPE_PCA_PILEUPCALLER_PLAN.md` 顶部「实际执行策略」。

---

## 1. 阶段总览

| 阶段 | 内容 | 状态 |
|---|---|---|
| Phase 0 | REF/ALT 方向体检 | ✅ 完成 |
| Phase A | 转 PLINK + 锁 A2=irgsp + MAF/geno | 🔶 3K 完成；720 待按 geno 0.20 重跑；Civán 复用补归一化 |
| Phase B | pileupCaller 替换调用 | 🔶 脚本已落地，spike 工具可用；批量铺开待做 |
| Phase C | 三档现代诊断图 + 古代投影图 | 🔶 绘图脚本已落地；smartpca 投影末尾失败（不追） |

---

## 2. 已完成

- Phase 0 三面板方向体检，产出 720 翻链清单 + mismatch 明细。
- Phase A：3K 完整跑通（29.6M 位点 → MAF 后 4.59M）；720 完成转换 + 锁方向 +
  缺失率审计，结论是把 geno 从 0.10 放宽到 0.20。
- pileupCaller 环境确认：v1.5.3.1（`~/software/pileupCaller-linux`）可用，
  v1.6.0.0 segfault 已弃用。
- pileupCaller 调用、PLINK→`.calls.txt` 转换、coverage 漏斗、modern-only 与
  smartpca `.evec` 绘图脚本全部落地并提交到 `scripts/ecotype_pca_v2/`。

---

## 3. 进行中 / 待做

- 720 按 `geno 0.20` 重跑 07 MAF（当前旧数字 51,549 → 45,378 是基于 geno 0.10）。
- Civán 补 Phase A 归一化（02 转换 / 锁 A2 / 29 转回 EIGENSTRAT）。
- Phase B：先 spike（1 样本 × Civán 1015-marker），跑通
  `mpileup | pileupCaller --randomHaploid → .calls.txt → merge → smartpca`，
  确认结果合理后再铺开 16 样本 × 3 面板。

---

## 4. 阻塞 / 已知问题

- 末尾 smartpca `lsqproject` 投影失败（古代样本被 `insufficient data` 剔出
  `.evec`）。本轮决定**不追**，只留档在 `PATH_MAP` 第 6 节；恢复时优先查各样本
  `coverage_funnel.tsv` 的 `ld_passed_covered` 是否为 0/极低。

---

## 5. 下一步（建议顺序）

1. 720 重跑 07（geno 0.20）→ 锁 A2 → 更新 `PHASEA_RESULTS.md` 数字。
2. Civán 补 Phase A 归一化。
3. Phase B spike（1 样本 × Civán），验证 pileupCaller 输出与 merge/投影链路。
4. spike 通过后铺开 16 × 3，再进入 Phase C 出图。
