# DATA & PROJECTION PLAN — angkor + nanzuo 全流程（供服务器 opencode 使用）

> 服务器工作目录：`/home/scratch/yinmt_rice/`（数据在 `data/`，仓库已解包）。
> 本文件 = 数据全貌 + panel 情况 + metadata 说明 + 样品情况 + 投影方案 + 执行步骤。
> 服务器新会话开工前先读本文件 + `HANDOFF_ITP_BESTHIT.md`。

---

## 1. data/ 数据全貌（已 symlink 好）

| 条目 | 内容 | 是否 asian_rice_panel extract？ | 在流程中的角色 |
|---|---|---|---|
| `angkor_shotgun_reads/` | 443 个 `{Library_ID}.prefiltered.IRGSP1.mapped.fq` | ❌ **供应商 IRGSP1 预筛**（非 panel） | angkor shotgun 候选 reads（360 angkor + 16 no-besthit 的 shotgun 库 + 对照） |
| `angkor_capture1_raw/` | 16 个 `*.bbduk.lowcomp_filtered.fq` | ❌ capture 原始（仅 bbduk 低复杂度过滤） | 16 no-besthit 的 capture panel1 原始 reads |
| `angkor_capture2_raw/` | 16 个 `*.bbduk.lowcomp_filtered.fq` | ❌ 同上 | capture panel2 原始 reads |
| `angkor_bt2_extract_fastq/` | 32 个 `{robot}_RicePanel{1,2}.bt2new.fastq.gz` | ✅ **bowtie2 -N1 vs asian_rice_panel 提取**（有 QC） | capture 候选 reads 优先用这个 |
| `nanzuo_extract/popgen/` | 21 个 `YWL1-A3xxx.bt2.primary_mapped.fastq.gz` | ✅ **bowtie2 -N1 vs asian_rice_panel 提取** | 南佐 popgen 子库候选 |
| `nanzuo_extract/shotgun/` | 21 个 `YWL1-A3xxx.bt2.primary_mapped.fastq.gz` | ✅ **bowtie2 -N1 vs asian_rice_panel 提取** | 南佐 shotgun 子库候选 |
| `nanzuo_extract/function/` | 21 个 `YWL1-A3xxx.bam2fq.fastq.gz` | ❌ **原始 BAM→fq**（非 panel，类 angkor shotgun） | 南佐 function 子库候选 |
| `nanzuo_mapped_bams/` | 21 个 `YWL1-A3xxx.dedup.bam`(+bai+stats) | 已 IRGSP 比对的 pooled BAM（我们 bowtie2+markdup 做的） | 参考/对照；besthit 不吃它 |
| `angkor_metadata/` | 元数据表（见 §3） | — | 出图 / PCA / 损伤 |
| `asian_rice_panel_index/` | `irgsp.fa`、`asian_rice_panel.fa` + bt2 索引 | — | IRGSP 比对 + 竞争比对参考 |
| `wgs_eukaryota/`（待补） | `/home/database/ref20250728/cph_euk`（131 库） | — | 竞争比对判种 DB |

> ⚠️ `wgs_eukaryota` symlink 还没建（命令见 §6 步骤 0）。

**候选 reads 结论**：capture 用 `angkor_bt2_extract_fastq`；angkor shotgun 用 `angkor_shotgun_reads`（供应商 IRGSP1 预筛，直接当候选）；南佐 popgen/shotgun 用 `nanzuo_extract` 的 panel 提取产物，function 用 bam2fq 产物。

---

## 2. panel 情况（投影用现代参考）

| 内容 | 路径（/home/scratch/yinmt202607/） | 说明 |
|---|---|---|
| 720 现代面板 EIGENSTRAT | `db/6.7M_720/asn720.6m.{snp,ind,geno}` | 6,769,714 SNP × 718 样，IRGSP1 坐标 |
| 全 6.7M A2 锁 PLINK | `gene/results/ecotype_pca_v2/phaseA/720/plink/asn720.6m.irgsp.{bed,bim,fam}` | pileupCaller 位点集 |
| 全 6.7M 参考 EIGENSTRAT | `gene/pileupcaller_work/ecotype_pca_v2_upload/full6M/eigenstrat/full6M.ref.{snp,eigenstratgeno,ind,poplistname.txt}` | smartpca lsqproject 参考（718 现代） |
| 15 群体标签 | IND 148 / OrA 133 / OrD 66 / OrB 64 / TRJ 48 / OrADM 41 / OrC 38 / AUS 38 / OrF 34 / TEJ 30 / OrE 24 / JAPONICA_UNSPEC 18 / ARO 15 / INTERMEDIATE_TYPE 12 / RAY 9 | 建轴群体 |

投影方式：`13_merge`（古 calls 拼进 718 现代矩阵）→ `smartpca lsqproject YES`（轴只由现代算，古样被动投影）。

---

## 3. angkor metadata 内容与目的

`angkor_metadata/`（=/home/scratch/yinmt202607/angkor_metadata/）：

| 表 | 行数 | 内容 | 用途 |
|---|---|---|---|
| `angkor_final_metadata.tsv` | 377（含表头） | 376 angkor：`sample_id, robot_sample_id, library_id, archive_sample_id, field_sample_id(core), site, lat, lon, country, project, depth_cm, age(CE), age_unit(错,勿用), prep(单/双链), data_source(shotgun-only / no-besthit)` | **主元数据**，robot 主键 |
| `plot_meta_440_full.tsv` | 440 | 在 final 基础上加：`base_robot, libtype(SG/C1/C2/pooled_*/shotgun-only), besthit, prep` | 出图索引（投影点→岩芯/年代/子库/besthit） |
| `site_coords.tsv` | 4 | site 经纬度 + 样本数 | R 位置图 |
| `library_prep.tsv` | 377 | 376 样 建库方法（单链 113 / 双链 263） | 损伤（damage）计算用 |
| `age_depth_model.tsv` | 年代模型 | 只覆盖 CAM23-11 / CAM23-13 / CAM2509 | 年龄来源；age 列=CE |
| `sample_meta_data_20250922.tsv` | 10710 | CGG 大表（字段源，勿删） | 原始字段核对 |
| 已删 | — | `angkor_master_metadata`、`clean_*`、`library_list.txt` 等 | 清理过 |

> 年龄：**CE（公元）**，`age_unit: ka BP` 是位点级错误字段勿换算。CAM22-08/CAM2201 无年龄模型（age 空正常）。

---

## 4. 样品情况

### angkor（376）
- 组成：**360 shotgun-only + 16 no-besthit（合并 shotgun+capture）**
- 16 no-besthit：CAM2509 8 个 + CAM23-11/13 8 个；每样有 shotgun 库（miss16）+ capture1 + capture2 子库
- 岩芯分布（按 440 表统计）：CAM2509 239 / CAM23-13 74 / CAM23-11 71 / CAM22-08 49 / CAM2201 7
- 年龄范围：CAM23-11 1493–2021 CE、CAM23-13 971–1477、CAM2509 972–1472（⚠️CAM23-11 的 2021 异常待查）
- 已投影批次（本机已完成）：376 → 全 6.7M 投影，440 样 400 成功

### nanzuo（21）
- 样品：YWL1-A3483 … YWL1-A3503
- 每样 3 子库：popgen（`2.nanzuo_popgen_yancheng`）、shotgun（`1.nanzuo_shotgun`）、function（`6.nanzuo_function_yancheng`，bam→fq）
- pooled BAM 已有（`nanzuo_mapped_bams/`）

---

## 5. 投影方案（目标：子库各自投影 + pooled 也投影）

**每个分析单元都出一个投影点，全部落到同一个 6.7M 面板的共享轴上**：

| 样品集 | 投影单元（每单元 1 点） |
|---|---|
| angkor shotgun-only（360） | `{robot}`（pooled=shotgun 单库） |
| angkor no-besthit（16） | `{robot}`（pooled 合并）、`{robot}_SG`、`{robot}_C1`、`{robot}_C2`、`{robot}_BH`（besthit pooled） |
| nanzuo（21） | `{id}`（pooled）、`{id}_popgen`、`{id}_shotgun`、`{id}_function` |

即：**shotgun / capture(C1,C2) / popgen / function 各自投影，同时 pooled 也投**。总投影点数 ≈ 360 + 16×5 + 21×4 = 524。

> 子库点调用数会低，投影噪/被剔除属正常，按低置信标注。besthit 流程跑完后，每个投影单元得到一份 `.calls.txt`，再按第 6 节步骤 F merge + lsqproject。

---

## 6. 执行步骤（服务器会话按此推进）

**Step 0 — 补 link + 环境确认**
```bash
D=/home/scratch/yinmt_rice/data
ln -s /home/database/ref20250728/cph_euk "$D/wgs_eukaryota"
ls "$D/wgs_eukaryota" | head -3
ls -d /datasets/caeg_dataset/taxonomy/ncbi/20250530 2>/dev/null && echo TAX_OK || echo TAX_NEED_PATH
ls /home/scratch/yinmt202607/db/asian_rice_panel_index/*.acc2taxid
```

**Step A — 竞争比对**（逐库候选 reads → bowtie2 对 131 wgs + panel + IRGSP，k 值按师兄）
- 输入：§1 各候选 reads（capture 用 bt2_extract，shotgun/nanzuo-function 用现有 fq）
- 脚本：`rice_adna_pipeline/scripts/oryza_besthit/submit_oryza_competitive_mapping.sh`
- 命名：`{angkor|nanzuo}_{id}.{shotgun|capture_panel1|capture_panel2|popgen|function}.name_sorted.bam`

**Step B — besthit 过滤**
- `oryza_besthit_damage_filter.py`（damage-window 5bp，Oryza 属 4527 动态解析，KEEP/REJECT）
- 输出各 unit 的 KEEP reads（FASTQ）

**Step C — IRGSP mapping（我们 bowtie2 -N1）**
- `real_use/mapping/map_irgsp_single_fq.sh`（bowtie2 `-k3 -L22 -N1 ...` + markdup）

**Step D — 去重（师兄 aedna_dedup.py）**
- 用师兄的去重替代我们 markdup（consensus / unclipped-boundary）

**Step E — pileupCaller 伪单倍体**（每投影单元 1 份 calls）
- 复用 `real_use/pileupcaller/`（sbatch_pileupcaller_array.sh + shared_call + plink_to_calls），位点集 = 全 6.7M

**Step F — merge + lsqproject**
- `13_merge_ancients_fixed_panel.py` → `14_run_fixed_smartpca.sh`（718 现代参考，lsqproject YES）

**Step G — 出图**
- `plot_angkor_dashboard.py`（交互面板，meta 用扩展后的 plot_meta）

**先 spike**：1 个 angkor + 1 个 nanzuo，跑通 A→E 全链再批量。

---

## 7. 未决问题（开工确认）

- nanzuo 参与子库：popgen/shotgun/function 全要？（§5 按全要列）
- 调用参数：mapq25/baseq25/seed派生（师兄） vs mapq25/baseq30/seed0（我们）——跟师兄一致就全用师兄的
- angkor 360 shotgun 是否全部走 besthit（还是只对重点样）——§5 按全走列
- CAM23-11 2021 CE 异常年龄来源