# 子库路径地图（SUB_LIBRARY_PATH_MAP）— angkor 16 no-besthit + nanzuo

> 用途：为「16 个 no-besthit 样各拆 4 个投影点（shotgun / capture1 / capture2 / pooled）」
> 以及 nanzuo 同款处理，定位**子库 reads → 提取 → irgsp mapping → 下游可用 BAM** 的完整链。
> 所有路径均为 2026-08-21 在服务器 `/home/scratch/yinmt202607/` 实测核对。
> 服务器：login01，数据根 `/home/scratch/yinmt202607/`；本地 git 仓：
> `/Volumes/SSD/claude_code/rice_adna_pipeline`（分支 `codex/nanzuo-popgen`）。

---

## 0. 一句话现状

- 16 个 no-besthit 样的三份子库 reads **全部在服务器上**，且都已有「可直接下游」的形式：
  - shotgun：`5.angkor_shotgun_finished/`（已 irgsp 预筛）
  - capture1 / capture2：`tests/param_matrix_bt2_vs_bwa/00.extraction/bt2_new/`（bowtie2 -N1 已提取）
- **pooled（合并）BAM 已存在且已投影**：`map_irgsp_no_besthit/{robot}.dedup.bam`
- 缺的只是 **shotgun-only / capture1-only / capture2-only 三个独立 BAM**（需单独 irgsp mapping）
- nanzuo 三子库提取产物 + 合并 + 最终 BAM 全在（各 21）

---

## 1. angkor 16 个 no-besthit 样 — robot ↔ 库 ↔ 岩芯 映射

来源：`angkor_metadata/angkor_final_metadata.tsv`（`data_source` = `no-besthit(合并shotgun+capture)`）。

| robot_sample_id | shotgun 库(=miss16) | archive_sample_id | 岩芯(field) |
|---|---|---|---|
| LV6000619499 | LV7008879541 | LV3006121032 | CAM2509 |
| LV6000619917 | LV7008879245 | LV3006120633 | CAM2509 |
| LV6000620016 | LV7008890373 | LV3006120665 | CAM2509 |
| LV6000620032 | LV7008890353 | LV3006120680 | CAM2509 |
| LV6000620166 | LV7008890404 | LV3006120708 | CAM2509 |
| LV6000620172 | LV7008890359 | LV3006105274 | CAM2509 |
| LV6000654686 | LV7008879263 | LV3006121044 | CAM2509 |
| LV6000654698 | LV7008879193 | LV3006121074 | CAM2509 |
| LV7008416272 | LV7008961181 | LV3006006419 | CAM23-13 |
| LV7008416280 | LV7008960345 | LV3006006421 | CAM23-13 |
| LV7008416294 | LV7009026578 | LV3006006423 | CAM23-13 |
| LV7008416329 | LV7008957621 | LV3006007172 | CAM23-11 |
| LV7008416339 | LV7009024536 | LV3006007306 | CAM23-13 |
| LV7008416349 | LV7009022707 | LV3006007308 | CAM23-13 |
| LV7008416379 | LV7008961151 | LV3006007198 | CAM23-11 |
| LV7008416407 | LV7009024089 | LV3006007231 | CAM23-11 |

`/tmp/miss16.txt` = 上表「shotgun 库」列（16 个），已在 shotgun reads 目录 16/16 命中。

---

## 2. angkor 16 样 — 子库 reads（可直接下游的形式）

| 子库 | 路径 | 状态 |
|---|---|---|
| **shotgun** | `5.angkor_shotgun_finished/data/reads/{Library_ID}.prefiltered.IRGSP1.mapped.fq`（Library_ID = 上表 shotgun 库列） | 已 irgsp 预筛，**直接可 mapping** |
| **capture1** | `tests/param_matrix_bt2_vs_bwa/00.extraction/bt2_new/fastq/{robot}_RicePanel1.bt2new.fastq.gz` | bowtie2(-N1, vs asian_rice_panel.fa) 已提取，**需 irgsp mapping** |
| **capture2** | `tests/param_matrix_bt2_vs_bwa/00.extraction/bt2_new/fastq/{robot}_RicePanel2.bt2new.fastq.gz` | 同上 |

提取产物计数：`00.extraction/bt2_new/fastq/` = **32 个**（16 样 × 2 panel），命名 `{robot}_RicePanel{1,2}.bt2new.fastq.gz`。

每样合成 reads（param 测试产物，含 shotgun+panel1+panel2 三个文件 / robot）：
`tests/param_matrix_bt2_vs_bwa/01.reads_combined/bt2_new/`（48 个文件 = 16×3）：
- `{robot}.prefiltered.IRGSP1.mapped.fq`（shotgun，robot 命名版）
- `{robot}_RicePanel1.bt2new.fastq.gz`
- `{robot}_RicePanel2.bt2new.fastq.gz`

---

## 3. angkor 16 样 — 已存在的 mapped BAM

| BAM | 路径 | 说明 |
|---|---|---|
| **pooled（合并 shotgun+capture）** | `gene/results/ecotype_pca/map_irgsp_no_besthit/{robot}.dedup.bam` | **当前 376 投影用的就是它**，=「4 点」里的 pooled 点 |
| 原始 besthit BAM | `gene/results/ecotype_pca/bam_irgsp/{robot}.besthit_oryza.irgsp.bam` | oryza besthit 过滤路线，历史 |
| param 测试最终 BAM | `tests/param_matrix_bt2_vs_bwa/02.final_mapping/{combo}/final/` | 7 组合（bt2_new/bt2_old/bwa 提取 × bt2new/bt2old/bwa 比对），16 样 × 组合 |
| shotgun-only 生产 BAM | `gene/results/ecotype_pca/map_irgsp_shotgun_only/{lib}.dedup.bam`（427 = 360 angkor + 67 对照） | **不含 16 个 no-besthit**（它们的 shotgun 库 = miss16，未单独投） |

---

## 4. nanzuo 21 样 — 子库 reads / 中间产物

| 阶段 | 路径 | 计数 |
|---|---|---|
| 原始 popgen | `2.nanzuo_popgen_yancheng/*.bbduk.lowcomp_filtered.fq` | — |
| 原始 shotgun | `1.nanzuo_shotgun/*.taxa_cleaned.fq.gz` | — |
| 原始 function | `6.nanzuo_function_yancheng/*.bam` | — |
| **bt2(-N1) 提取 popgen** | `nanzuo/00.extract/popgen/{s}.bt2.primary_mapped.fastq.gz` | 21 |
| **bt2(-N1) 提取 shotgun** | `nanzuo/00.extract/shotgun/{s}.bt2.primary_mapped.fastq.gz` | 21 |
| **function bam→fq** | `nanzuo/00.extract/function/{s}.bam2fq.fastq.gz` | 21 |
| 合并 | `nanzuo/01.merge/{s}.combined.fastq.gz` | 21 |
| **最终 BAM** | `nanzuo/02.map_irgsp/{s}.dedup.bam` | 21 |

nanzuo 三子库（popgen/shotgun/function）各自已是提取后 fastq，**可直接单独 irgsp mapping**
（沿用 `scripts/nanzuo_popgen/map_irgsp_single.sh` 的 bowtie2 -N1 + markdup 流程）。

---

## 5. 参考 / 工具（mapping + projection 复用的固定件）

| 内容 | 路径 |
|---|---|
| irgsp 主参考 + bt2 索引 | `db/asian_rice_panel_index/irgsp.fa`、`irgsp_bt2idx` |
| 提取用多物种 panel + bt2 索引 | `db/asian_rice_panel_index/asian_rice_panel.fa`（`.bt2l` 大索引） |
| pileupCaller | `~/software/pileupCaller-linux`（v1.5.3.1，v1.6 禁用） |
| 全 6.7M A2 锁 bfile | `gene/results/ecotype_pca_v2/phaseA/720/plink/asn720.6m.irgsp.{bed,bim,fam}` |
| 全 6.7M shared 调用文件 | `gene/pileupcaller_work/ecotype_pca_v2_upload/full6M/calls/shared.snp` / `shared.sites.bed` |
| 全 6.7M 参考 EIGENSTRAT | `.../full6M/eigenstrat/full6M.ref.{snp,eigenstratgeno,ind,poplistname.txt}` |
| 16 样 q25 调用 | `.../full6M/calls/q25/{robot}.calls.txt`（已在 376 投影中） |

---

## 6. 下一步（4 点投影，待方案确认后执行）

1. **单独 mapping**：16 样 × 3 子库（shotgun / capture1 / capture2）= 48 个新 BAM。
   输入 reads 见第 2 节，参数沿用 `map_irgsp_single.sh`（bowtie2 -N1 → markdup 只标记）。
2. **q25 调用**：复用 `full6M/calls/shared.snp` + `sbatch_pileupcaller_array.sh`，样品 ID 用
   唯一后缀（如 `{robot}_SG / _C1 / _C2`）避免与 376 及彼此冲突。
3. **重新 merge + smartpca**（全 6.7M，~40 min/次）→ 新 `.evec`。
4. **出图**：每个 no-besthit 样 4 个点（shotgun/capture1/capture2/pooled），
   nanzuo 同款（popgen/shotgun/function/pooled）。

---

## 7. 临时清单（/tmp，重开窗口可能还在）

- `/tmp/angkor376_robot.txt` — 376 robot（投影样品清单）
- `/tmp/miss16.txt` — 16 个 no-besthit 的 shotgun 库（=第 1 节 shotgun 库列）
- `/tmp/angkor_calls_n.tsv` — 376 样品 non-9 调用数
- `/tmp/projected.txt` / `/tmp/all376.txt` — 投影出 340 / 全 376
