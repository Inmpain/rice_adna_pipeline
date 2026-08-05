# Oryza competitive mapping / best-hit 接手说明

更新时间：2026-08-05

## 1. 当前目标

候选水稻 reads 已经由 shotgun、panel1、panel2 三份 FASTQ 合并完成。当前正在把每个样本分别比对到 WGS 真核数据库、亚洲水稻 panel 和 IRGSP。下一阶段需要根据所有比对结果执行 best-hit，排除更像非水稻物种的 reads。

## 2. 关键路径

### GitHub

仓库：

https://github.com/Inmpain/rice_adna_pipeline

当前分支：

`codex/oryza-competitive-mapping`

mapping 脚本：

`scripts/oryza_besthit/submit_oryza_competitive_mapping.sh`

测试脚本：

`scripts/oryza_besthit/test_submit_oryza_competitive_mapping.sh`

### 本机参考代码

原 best-hit 项目：

`/Users/inmpain/github/aeDNA_popgen`

WGS mapping/动态内存参考：

`/Users/inmpain/github/Pipeline_snakemake/new_single_multi`

本机 GitHub 工作副本：

`/Users/inmpain/Documents/angkor/rice_adna_pipeline_publish`

### 服务器

执行脚本：

`/home/scratch/yinmt202607/gene/scripts/submit_oryza_competitive_mapping.sh`

候选 FASTQ：

`/home/scratch/yinmt202607/gene/results/oryza_candidates_combined`

格式：

`<sample>.oryza_candidates.combined.fastq.gz`

mapping 输出根目录：

`/home/scratch/yinmt202607/gene/results/oryza_competitive_mapping`

各数据库 BAM：

`.../bam_by_database/<sample>/<sample>.<database>.bam`

每样本最终合并 BAM：

`.../by_sample/<sample>.competitive.name_sorted.bam`

日志：

`.../logs`

提交记录：

`.../submissions`

串行运行状态：

`.../series`

## 3. 数据库

WGS 真核数据库：

`/home/database/ref20250728/cph_euk/wgs_eukaryota.1.fas.gz`

至：

`/home/database/ref20250728/cph_euk/wgs_eukaryota.129.fas.gz`

亚洲水稻 panel：

`/home/scratch/yinmt202607/db/asian_rice_panel_index/asian_rice_panel.fa`

IRGSP：

`/home/scratch/yinmt202607/db/asian_rice_panel_index/irgsp_bt2idx`

合并 accession-taxid：

`/home/scratch/yinmt202607/db/asian_rice_panel_index/all_wgs_asian_irgsp.acc2taxid`

NCBI taxonomy：

`/home/database/ref20250728/taxonomy_CPH/ncbi/20250530`

目标 Oryza taxid：

- 4529：Oryza rufipogon
- 4530：Oryza sativa
- 4536：Oryza nivara

## 4. 当前 mapping 流程

每样本比对131个数据库：

- 129个 WGS shard
- asian_rice_panel
- IRGSP

Bowtie2 参数：

```bash
-k 100
-L 22
-i S,1,1.15
--mp 1,1
--rdg 0,1
--rfg 0,1
--score-min L,0,-0.1
--no-unal
