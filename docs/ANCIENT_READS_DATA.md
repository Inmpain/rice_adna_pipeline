# 水稻古 DNA 数据组织（angkor + nanzuo）

最后更新：2026-08-19

## 1. 概览

两批古水稻样品，均比对到 IRGSP 参考基因组后做 ecotype PCA 投影：

| 批次 | 样品 | reads 类型 | 分析状态 |
|---|---|---|---|
| angkor（吴哥窟） | **16 个** | shotgun + capture1 + capture2 | ✅ 已跑通（besthit / no-besthit 两条线） |
| angkor（吴哥窟） | **几百个 shotgun-only** | 仅 shotgun | 🔶 待跑（本批次） |
| nanzuo（南佐） | **21 个**（20 有效） | shotgun + popgen + function（合并） | ✅ 已跑通（无 besthit） |

## 2. angkor 16 个分析样本

**样本 ID（robot_sample_id）**：`LV6000619499` … `LV7008416407`（16 个，`YWL` 前缀的才是南佐）。

**三个 reads 来源**（→ 合并 → bt2new 提取候选 reads → [besthit / 不 besthit] → map irgsp）：

| 来源 | 路径 | 命名 | 说明 |
|---|---|---|---|
| shotgun | `5.angkor_shotgun_finished/data/reads/` | `{Library_ID}.prefiltered.IRGSP1.mapped.fq` | 预筛选到 IRGSP1 的 shotgun reads |
| capture1 | `3.angkor_capture_panel1/data/reads/` | `{robot_sample_id}_RicePanel1.bbduk.lowcomp_filtered.fq` | 原始 cap1（未提取） |
| capture2 | `7_angor_capture_panel2/data/reads/` | `{robot_sample_id}_RicePanel2.bbduk.lowcomp_filtered.fq` | 原始 cap2（未提取） |

**提取后的候选 reads**（besthit 的输入，也是 no-besthit 直接 map 的输入）：
`tests/param_matrix_bt2_vs_bwa/01.reads_combined/{bwa,bt2_new}/`
- shotgun: `{robot_sample_id}.prefiltered.IRGSP1.mapped.fq`
- cap1: `{robot_sample_id}_RicePanel1.{bwa.primary_mapped|bt2new}.fastq.gz`
- cap2: `{robot_sample_id}_RicePanel2.{bwa.primary_mapped|bt2new}.fastq.gz`

**元数据**：`/home/scratch/yinmt202607/angkor_robot_library.txt`
- `robot_sample_id`（第 8 列）= 我们用的样本名（LV6000619499 等）
- `Library_ID`（第 2 列）= shotgun 文件命名用
- `Depth`（第 4 列，cm）/ `Age`（第 5 列，公历年代）

**注意**：`Library_ID` ≠ `robot_sample_id`（是两套编号，要经 metadata 桥接）。
MCP proxy 样本（`LV{数字}-LV{数字}-proxy`）不属于 angkor 16 样本，统计时显式排除。

## 3. angkor shotgun-only 样品（几百个，已定位）

- **位置**：`5.angkor_shotgun_finished/data/reads/*.prefiltered.IRGSP1.mapped.fq`
- **数量**：443 个 shotgun 库在 CGG 大表（`sample_meta_data_20250922.tsv`）里匹配到 metadata；
  排除 18 个 SmplNTC（阴性对照）+ 43 个空白 field_sample + 6 个 Fuglsø（丹麦，噪声）后，
  **保留 376 个 angkor shotgun 库**。
- **只有 shotgun reads**，无 capture → 覆盖比 16 个分析样本低。
- **流程**：shotgun reads → bowtie2 map irgsp（`map_one_sg.sh`，sbatch 并行）→ 覆盖普查 → 调用 → 投影。
- **元数据**：`angkor_metadata/angkor_shotgun_meta_clean.tsv`（library_id / robot_sample_id /
  archive_sample_id / field_sample_id / site_name / depth / age / prep）
  - 从 `sample_meta_data_20250922.tsv`（CGG 全球 eDNA 大表，10710 库）筛出 angkor 部分；
  - depth/age 从 `age_depth_model.tsv` 按 archive_sample_id join。

### shotgun 库的位点分布（4 处 Angkor 位置）

| 位点 | 库数 | depth/age |
|---|---|---|
| CAM23-11/13 北护城河（Northern Temple Moat, Angkor Wat） | 113 | ✅ 0.5–184 cm + age |
| CAM2509 北寺塘（Northern Temple pond, Angkor Wat） | 207 | ❌ 无（库缺 archive_sample_id） |
| CAM22-08 Angkor Thom 护城河 | 49 | ❌ |
| CAM2201 West Baray | 7 | ❌ |

**结论**：angkor shotgun 不是单一位置——至少来自王城北护城河、北寺塘、Angkor Thom 护城河、West Baray 四处。
- 只有 **CAM23-11/13（113 库）有 depth/age**（含已分析的 16 个）；
- 其余 263 库（CAM2509/22-08/2201）无 depth/age，标注 NA。
- **排除项**：SmplNTC（阴性对照 18）、空白（43）、Fuglsø（6，非 angkor）。

## 4. nanzuo（南佐）

- 21 个样本（`YWL1-A3483` … `YWL1-A3503`），`YWL1-A3495` 为坏文库（剔除）。
- 三条提取源（popgen/shotgun/function）合并 → `nanzuo/01.merge/{sample}.combined.fastq.gz`
- → bowtie2 map irgsp（无 besthit）→ `nanzuo/02.map_irgsp/{sample}.dedup.bam`
- 无深度/年代元数据。

## 5. 分析产物

| 产物 | 说明 |
|---|---|
| `master_table.tsv` / `master_table_nb.tsv` | besthit / no-besthit 全样本汇总表（reads/覆盖/调用/损伤/最近群体） |
| `merge/combined.groups.png(.html)` | besthit angkor + nanzuo 投影图 |
| `nobesthit/merge_nb/all_nobesthit.pca.evec` | no-besthit 新 panel 投影 |
| `metaDMG_summary.tsv` | 损伤汇总（besthit angkor + nanzuo） |
| `nobesthit_dmg.tsv` | no-besthit angkor 损伤 |
