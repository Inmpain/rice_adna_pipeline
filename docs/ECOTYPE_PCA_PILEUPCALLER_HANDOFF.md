# Ecotype PCA v2（pileupCaller 版）— 交接手册 / Runbook

> 本文是一份**可照着做的交接手册**，服务三件事：① 分析新数据；② 继续当前流程；
> ③ 出了问题时定位 bug。按顺序读即可。
>
> 配套文档：路径看 `ECOTYPE_PCA_PILEUPCALLER_PATH_MAP.md`，进度看
> `ECOTYPE_PCA_PILEUPCALLER_PROGRESS.md`，数字结果看
> `ECOTYPE_PCA_PILEUPCALLER_RESULTS.md`，策略/决策看
> `ECOTYPE_PCA_PILEUPCALLER_PLAN.md` 顶部的「实际执行策略」。

最后更新：2026-08-18

---

## 0. 一句话现状

720 面板的 **共享轴投影已经跑通并验证**：16 个古样本全部投影进 `.evec`（call 数 21–1314），
不是之前的全 9。剩余：3K 面板同样流程、私有轴、出图/收尾。

---

## 1. 环境与工具（服务器侧，已就绪）

```bash
source activate /home/usr/yinmt/.local/mamba/snakemake   # plink2 / python3 / matplotlib
module load samtools                                      # samtools mpileup
PILEUP_CALLER=~/software/pileupCaller-linux                # v1.5.3.1（v1.6.0.0 segfault，别用）
```

脚本最省事是克隆仓库（所有脚本都在）：

```bash
git clone -b codex/ecotype-pca-pileupcaller https://github.com/Inmpain/rice_adna_pipeline.git rice_adna_pipeline
```

或按需单独 curl（raw 地址格式）：

```text
https://raw.githubusercontent.com/Inmpain/rice_adna_pipeline/codex/ecotype-pca-pileupcaller/scripts/ecotype_pca_v2/<文件名>
```

---

## 2. 冻结参数（不要硬编码进脚本，读 config）

| 参数 | 值 |
|---|---|
| MAF | 0.01（三面板统一） |
| geno | 3K=0.05 / 720=0.20 / Civán=0.05 |
| track | ALL（本轮不做 TV） |
| LD | window 100kb / r2 0.20 |
| ancient MAPQ / BaseQ | 30 / 30 |
| smartpca | `lsqproject: YES`, `numoutlieriter: 0`, `num_pcs: 10` |
| 参考等位基因方向 | `plink --a2-allele <ref> 2 1 --keep-allele-order`（每次写新 bed 后重锁） |

---

## 3. 完整可执行流程（分析新数据照这个跑）

### 3.1 准备 marker bfile 与参考 EIGENSTRAT（Phase A 已产出则跳过）

```text
marker bfile        : <BFILE>.{bed,bim,fam}     # MAF/geno 过滤后、A2 已锁 irgsp
reference.snp       : 与 marker bfile 同一位点集的 EIGENSTRAT .snp（29 转出）
reference.eigenstratgeno / reference.ind / reference.poplistname : 现代参考三件套
```

### 3.2 生成 pileupCaller 输入（chr 归一化在这里发生）

`pileupcaller_shared_call.sh` 会从 `--bfile` 自动生成 `.snp` 和 `.sites.bed`，并做
`1 → chr01` 归一化。**不要自己另写 awk，直接用这份脚本。**

### 3.3 批量调用 16 个古样本

```bash
cd pileupcaller_work
PILEUP_CALLER=~/software/pileupCaller-linux
REF=/home/scratch/yinmt202607/db/asian_rice_panel_index/irgsp.fa
BAMDIR=/home/scratch/yinmt202607/gene/results/ecotype_pca/bam_irgsp
BFILE=<你的 marker bfile 前缀>

SAMPLES="LV6000619499 LV6000619917 LV6000620016 LV6000620032 LV6000620166 \
LV6000620172 LV6000654686 LV6000654698 LV7008416272 LV7008416280 LV7008416294 \
LV7008416329 LV7008416339 LV7008416349 LV7008416379 LV7008416407"

mkdir -p calls
for S in $SAMPLES; do
  ./pileupcaller_shared_call.sh \
    --bam "$BAMDIR/$S.besthit_oryza.irgsp.bam" \
    --sample "$S" --bfile "$BFILE" --ref-fasta "$REF" \
    --mapq 30 --baseq 30 --seed 0 \
    --out-dir calls --label "$S" || { echo "FAIL $S"; continue; }
  python3 pileupcaller_plink_to_calls.py --bfile "calls/$S" --out "calls/$S"
done

# 16 个样本 non-9 计数（不应有 0）
for S in $SAMPLES; do
  n=$(awk '$1!="9"{c++} END{print c+0}' "calls/$S.calls.txt"); echo "$S $n"
done
```

### 3.4 顺序核对（merge 前必做，否则会静默错位）

```bash
D=<reference.snp 所在目录>
diff <(cut -f2 calls/LV7008416379.bim) <(awk '{print $1}' "$D/reference.snp") \
  && echo ORDER_OK || echo ORDER_MISMATCH
```

> 注意：EIGENSOFT `.snp` 是**空格分隔**，用 `awk '{print $1}'`，别用 `cut -f1`。

### 3.5 merge → smartpca → 出图

```bash
python3 13_merge_ancients_fixed_panel.py \
  --reference-geno "$D/reference.eigenstratgeno" \
  --reference-ind  "$D/reference.ind" \
  --fixed-snp      "$D/reference.snp" \
  $(for S in $SAMPLES; do printf -- '--calls %s ' "$S=calls/$S.calls.txt"; done) \
  --ancient-poplabel Ancient --label <panel>.v2 --out-dir merge

./14_run_fixed_smartpca.sh --config ecotype_pca_v2.yaml \
  --geno merge/<panel>.v2.merged.eigenstratgeno --snp "$D/reference.snp" \
  --ind merge/<panel>.v2.merged.ind --poplist "$D/reference.poplistname" \
  --label <panel>.v2.pca --out-dir merge

python3 plot_smartpca_evec.py \
  --evec merge/<panel>.v2.pca.evec --eval merge/<panel>.v2.pca.eval \
  --ind merge/<panel>.v2.merged.ind --nmarkers <marker数> \
  --title "<panel> shared projection (pileupCaller, ALL, MAF)" \
  --out-prefix merge/<panel>.v2.final
```

---

## 4. 已知 bug / 陷阱对照表（找 bug 先查这里）

| # | 现象 | 根因 | 修复 / 检查 |
|---|---|---|---|
| B1 | 古代样本全 9、`.evec` 里消失 | 染色体命名不一致（`1` vs `chr01`），mpileup/pileupCaller 静默查不到 | 用 `pileupcaller_shared_call.sh`，它已把 `.snp` 第2列和 `.sites.bed` 第1列归一化为 `chr01` |
| B2 | pileupCaller 报 `NonMissingCalls>0`，但 `.calls.txt` 全是 9 | plink2 `--export A` 的 `.raw` 列名带等位基因后缀（`1np2833_T`），转换脚本精确匹配落空 | 用新版 `pileupcaller_plink_to_calls.py`（按最后一个 `_` 切开匹配 + 按 A2=REF 定向） |
| B3 | `diff`/`cut -f` 报整文件顺序不一致 | EIGENSOFT `.snp` 是空格分隔，`cut -f` 按 TAB 切会把整行当一个字段 | 空格分隔文件一律用 `awk '{print $N}'` |
| B4 | merge PASS 但结果错位 | `13_merge` 只查行数、不查 SNP ID 顺序 | merge 前跑 3.4 的顺序核对 |
| B5 | 文件变成 `.xxx`（带点）/ `No such file` | `$PREFIX`/`$BFILE` 为空、或 `cd` 失败后 pwd 错 | 显式赋值前缀；脚本已加 `[[ -s file ]]` 前置校验 |
| B6 | 某古样本 call 数极低、投影噪声大 | 覆盖本身少（如 `LV7008416294` 仅 21） | 不是 bug；单独标记为低置信，别和其他 15 个同等解读 |

### 4.1 投影失败排查顺序（遇到「古代样本消失」时）

1. 各样本 `calls/*.calls.txt` 的 `non9` 是否为 0。
2. `.snp` 第 2 列、`.sites.bed` 第 1 列是否都是 `chr01…chr12`，与 BAM `@SQ` 一致。
3. `samtools view -H BAM | grep '^@SQ'` 核对命名。
4. merge 前是否做了 SNP ID 顺序核对。

---

## 5. 服务器路径速查

```text
results_v2_root = /home/scratch/yinmt202607/gene/results/ecotype_pca_v2
phase0          = .../phase0/            # REF/ALT 方向体检 + 720 翻链/mismatch 清单
phaseA/720      = .../phaseA/720/        # 转换/锁方向/MAF 中间文件
720 参考 EIGENSTRAT = .../phaseA/720/eigenstrat/720hybrid.{snp,eigenstratgeno,ind,poplistname.txt}
16 ancient BAM  = /home/scratch/yinmt202607/gene/results/ecotype_pca/bam_irgsp/*.besthit_oryza.irgsp.bam
irgsp.fa        = /home/scratch/yinmt202607/db/asian_rice_panel_index/irgsp.fa
pileupCaller    = ~/software/pileupCaller-linux
```

完整四层路径见 `ECOTYPE_PCA_PILEUPCALLER_PATH_MAP.md`。

---

## 6. 当前进度 + 下一步

已完成：

- Phase 0（3K/Civán 整体反标，720 混合方向 + 翻链/mismatch 清单）。
- Phase A（3K 完整、720 转换/锁方向/缺失率审计，结论 geno 0.20）。
- pileupCaller 调用 + PLINK→calls 转换 + merge + smartpca 全链路，720 共享轴投影验证通过。

下一步（按序）：

1. **3K 面板**照第 3 节跑同样流程（先 1 样本 spike，再 16 样本批量）。
2. 主分析：复用 v1 `scripts/ecotype_pca/run_sample_panel_pca.sh` 做 per-sample 私有轴。
3. 用 `plot_smartpca_evec.py` 出最终图，人工核对现代结构 + 古样本投影位置。
4. 单独记录低覆盖样本（如 `LV7008416294`）为低置信。
5. 提交当前未提交的 docs（`PROGRESS`/`RESULTS`/`PLAN`/本文件）。

---

## 7. 一句话记忆点

- **投影不改变解释度**：`lsqproject` 轴和 eigenvalue 只由现代参考算，古样本被动投影。
- **方向**：REF/ALT 统一锁到 irgsp（A2=REF），pileupCaller 会按参考基因组自动对齐。
- **顺序**：`.calls.txt` 顺序必须等于 reference `.snp` 顺序，merge 前必查。
- **分隔符**：EIGENSOFT 文本文件是空格分隔，工具链里别混用 `cut -f` 和 `awk`。
