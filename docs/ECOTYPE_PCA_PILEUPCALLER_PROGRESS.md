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
| Phase B | pileupCaller 替换调用 | ✅ 720 共享投影已跑通验证（44,920 marker，16 样本 call 21–1314）；3K 待跑同样流程 |
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

- 两张 720 图 PC 解释度差异（21.71/19.70 vs 5.70/6.03）待查根因——已定位为
  图来源不同（`plot_panel_pca.py` modern-only plink2 诊断 vs `plot_smartpca_evec.py`
  smartpca 投影），不是 `lsqproject` 投影改解释度，见 `HANDOFF` 第 4 节 B7。

---

## 5. 下一步（建议顺序）

1. 查清两张 720 图 PC 解释度差异根因（先比两个 `.eval`，见 `HANDOFF` B7）。
2. 3K 面板照 `HANDOFF` 第 3 节跑 pileupCaller 共享投影（先 1 样本 spike，再 16 批量）。
3. 主分析私有轴：复用 v1 `run_sample_panel_pca.sh`（命令见 `PHASE1_COMMANDS.md` 第 7 节）。
4. 出最终图 + 人工核对现代结构、古样本投影位置；LV7008416294 标低置信。
5. Civán 补 Phase A 归一化（02 转换 / 锁 A2 / 29 转回 EIGENSTRAT）。
