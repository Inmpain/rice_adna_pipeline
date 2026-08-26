# HANDOFF_ITP_BESTHIT — 转向「师兄 besthit 流程」工作交接

> 用途：把 angkor + nanzuo 的 reads 处理从「bowtie2 提取直投」切换为
> 「从 besthit 开始走师兄流程」，在**服务器上的 opencode 新会话**里继续。
> 断点：本机 Mac 的 opencode 会话结束；服务器 login01 上新建工作目录开新会话。
> 服务器：login01，`/home/scratch/yinmt202607/`。

## 0. 一句话现状（2026-08-26 交接时）

- **已完成**：376 angkor shotgun → 全 6.7M 面板 q25 调用 → merge → smartpca lsqproject 投影跑通；
  加 16 个 no-besthit 的子库 4 点（SG/C1/C2/BH）后共 **440 样、400 成功投影**（340 原始 + 60 子库/BH）。
  产物在 `gene/pileupcaller_work/ecotype_pca_v2_upload/full6M/merge/full6M.pca2.*`。
- **转向决策（用户已定）**：reads 处理**从 besthit 开始走师兄流程**，但：
  - 路径/目录用自己的（新建 ITP 工作目录，symlink 回 scratch）
  - **mapping 继续用我们的 bowtie2 -N1**（不换 bwa）
  - **去重用师兄的 `aedna_dedup.py`**，其余（竞争比对/besthit/选 Oryza reads/pileupCaller/合并 PLINK）用师兄的
  - **nanzuo + angkor pooled 到一起跑 besthit**（打标区分），跑完再按前缀拆
- **我们自己的 besthit 版本已在本服务器对 16 样跑过正式批**（留存率 13.72%）：
  `rice_adna_pipeline` `codex/oryza-competitive-mapping` 分支 `scripts/oryza_besthit/`。
  DB（wgs_eukaryota 131 库、asian_rice_panel、IRGSP、acc2taxid、taxonomy）都在 login01。

## 1. 师兄 besthit 流程（要走的链路）

```
逐库原始 reads
  → 竞争比对（bowtie2 对 131 wgs_eukaryota + asian_rice_panel + IRGSP）   [submit_oryza_competitive_mapping.sh]
  → besthit 过滤（oryza_besthit_damage_filter.py：damage-window 5bp，Oryza 属 4527 动态解析，KEEP/REJECT）
  → 选 Oryza reads（select_oryza_competitive_reads.py + extract_fastq_by_read_names.py）
  → IRGSP mapping（我们用 bowtie2 -N1）
  → 去重（师兄 aedna_dedup.py，consensus/unclipped-boundary）
  → pileupCaller 伪单倍体 → 合并 PLINK
```

参考脚本（ref_pipeline，本地 `/Volumes/SSD/claude_code/rice/ref_pipeline/`）：
- `zhe_pipeline/snakemake.filtered_bam_besthit.best10.optimized.260730.smk`（besthit）
- `zhe_pipeline/snakemake.extract_Oryza_reads.best10.260730.smk`（选 Oryza reads）
- `zhe_pipeline/snakemake.Oryza_reads.mapping.260731.smk`（IRGSP mapping + q25 + dedup）
- `zhe_pipeline/script/aedna_dedup.py`（去重，我们要用这个）
- `zhe_pipeline/Snakefile.pseudohaploid.from_panel.asn720_dense.260802.smk`（pileupCaller + 合并 PLINK）

我们自己的版本：`rice_adna_pipeline` `codex/oryza-competitive-mapping`：
- `scripts/oryza_besthit/submit_oryza_competitive_mapping.sh`（竞争比对，bowtie2，131 库）
- `scripts/oryza_besthit/oryza_besthit_damage_filter.py`（v2，Oryza 属动态，damage-window 5bp）
- `scripts/oryza_besthit/rebuild_all_wgs_asian_irgsp_acc2taxid.sh`（acc2taxid 重建，taxid 反标已修配方）
- `docs/ORYZA_BESTHIT_HANDOFF.md`（详细接手说明，含 16 样留存率表格）

## 2. 新工作目录（ITP → scratch）

用户打算：在 ITP 下新建文件夹，link 回 scratch。**具体路径待用户确认**
（「ITP 文件夹」指的是 `/itp/...` 还是新 scratch 子目录，服务器会话开工前先问清）。

建议结构（待确认后定稿）：
```
<workdir>/
├── data/                  # symlink 回 scratch 的原始数据/中间产物
│   ├── angkor_shotgun/    -> 5.angkor_shotgun_finished/data/reads/（360 + 16 no-besthit 的 shotgun 库）
│   ├── angkor_capture/    -> 3.angkor_capture_panel1、7_angor_capture_panel2（原始 fq）
│   ├── angkor_bt2_extract/-> tests/param_matrix_bt2_vs_bwa/00.extraction/bt2_new/fastq/（capture 提取产物）
│   ├── nanzuo_extract/    -> nanzuo/00.extract/{popgen,shotgun,function}/
│   └── metadata/          -> angkor_metadata/（angkor_final_metadata.tsv、plot_meta_440_full.tsv 等）
├── scripts/               # 拷贝 real_use/ + oryza_besthit/ + 师兄关键脚本
├── db/                    # symlink: asian_rice_panel、wgs_eukaryota、IRGSP、acc2taxid、taxonomy
└── results/               # 新流程产物
```

## 3. 各类型样本的 reads 来源（besthit 输入 = Oryza 候选 reads）

| 样本类型 | 候选 reads 来源 | 说明 |
|---|---|---|
| angkor shotgun（360 + 16 no-besthit 的 shotgun 库） | `5.angkor_shotgun_finished/.../{lib}.prefiltered.IRGSP1.mapped.fq` | 供应商已 IRGSP 预筛，**无 panel QC**，直接当候选喂竞争比对 |
| angkor capture1/2（16 no-besthit） | `00.extraction/bt2_new/fastq/{robot}_RicePanel{1,2}.bt2new.fastq.gz`（有 QC）或原始 `3./7_.angkor_capture_panel{1,2}/` | 有提取+QC 的用提取产物 |
| nanzuo（21 样） | `nanzuo/00.extract/{popgen,shotgun,function}/{s}.bt2.primary_mapped.fastq.gz`（有 QC） | 三子库；用户提到「nanzuo 两个 capture 中的某个」需开工时确认到底哪几个子库参与 |

> besthit 的竞争比对输入是**逐库 reads**（按 library 或 assay 分开），不是合并后的。
> pooled 打标：样本名加前缀 `angkor_` / `nanzuo_`（匹配师兄 BAM 正则
> `^(.+)\.(shotgun|capture_panel1|capture_panel2)\.name_sorted\.bam$`），besthit 跑完按前缀拆回。

## 4. 关键差异决策（已定，勿改）

- mapping：bowtie2 `-k 3 -L 22 -N 1 -i S,1,1.15 --mp 1,1 --rdg 0,1 --rfg 0,1 --score-min L,0,-0.1 --no-unal`
  （就是我们 `real_use/mapping/map_irgsp_single_fq.sh` 的参数）；竞争比对那步的 k 值按师兄的来
- 去重：师兄 `aedna_dedup.py`（不用我们 markdup 只标记）
- 其余（竞争比对/besthit/选 reads/pileupCaller/合并 PLINK）：师兄流程
- baseQ/seed 等调用参数：跟师兄 `Snakefile.pseudohaploid`（mapq25/baseq25/seed 派生）还是我们（baseq30/seed0），开工确认

## 5. 本地已备好可直接用的东西

- `real_use/`（`/Volumes/SSD/claude_code/rice/real_use/`）：我们跑通的全部脚本 + `README.md`（含与 ref_pipeline 的方法对比表）
- `ref_pipeline/`（`/Volumes/SSD/claude_code/rice/ref_pipeline/`）：师兄全流程脚本（zhe_pipeline/original_input/running_log/paleoclimate）
- 仓库：`rice_adna_pipeline`（codex/nanzuo-popgen + codex/oryza-competitive-mapping）、`panel-pca-pipeline`（main）

## 6. 服务器会话开工清单（新 opencode 会话第一件事）

1. 确认新工作目录路径（ITP?）与 symlink 方案
2. `curl` 拉取：`rice_adna_pipeline`（nanzuo-popgen + oryza-competitive-mapping 分支）、`panel-pca-pipeline`（main）、`real_use/`、`ref_pipeline/` 所需脚本
3. 确认 DB 可读：`wgs_eukaryota`（`/home/database/ref20250728/cph_euk`）、`asian_rice_panel`、`IRGSP`、`acc2taxid`、taxonomy（nodes.dmp/names.dmp）
4. 按第 3 节准备各类型样本的候选 reads（symlink）
5. 先小样本（如 1-2 个 angkor + 1 个 nanzuo）spike 跑通「竞争比对 → besthit → 选 reads → IRGSP bowtie2 → aedna_dedup」全链，再批量

## 7. 未决问题（开工先确认）

- 新工作目录确切路径（ITP?）
- nanzuo 到底哪几个子库参与（popgen/shotgun/function 全要，还是用户说的「两个 capture 之一」）
- 调用参数（mapq25/baseq25/seed派生 vs mapq25/baseq30/seed0）
- angkor 的 360 shotgun 是否全部要走 besthit（还是只对 16 个重点样走，其余维持现状）