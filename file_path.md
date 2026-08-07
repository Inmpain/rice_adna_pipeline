# 完整路径清单 (v5 — 最终权威版，供工作交接使用)

> 本文档是整个项目的路径地图，每个目录都标注了**是什么、由哪个脚本产出、
> 什么场景下需要用它**。新接手的人应该先读这份文档，再去看
> `PROJECT_STATUS.md`(项目当前进度)和`decisions_log.md`(关键决策原因)。
>
> ⚠️ **本仓库还有3个未合并进main的工作分支**（`codex/oryza-screen-merge`、
> `codex/oryza-competitive-mapping`、`codex/ecotype-pca-panel`），只读
> main分支看不到完整现状。四分支横向汇总见 `docs/REPO_OVERVIEW_STATUS.md`。
> 科研目标/证据阶梯/实施阶段的项目级框架见 `docs/RESEARCH_ROADMAP.md`。

最后更新: 2026-08-07

---

## 一、服务器目录总览(`/home/scratch/yinmt202607/`)

```
/home/scratch/yinmt202607/
├── asn720data/                  # 720份现代/近现代品种PLINK面板(外部数据)
├── db/                          # 参考基因组、索引、注释、外部数据库
├── results/                     # 全部分析产出(按处理阶段编号)
├── script/                      # 早期独立脚本(已同步进git的server_originals)
├── 3.angkor_capture_panel1/     # 原始数据: capture panel1
├── 7_angor_capture_panel2/      # 原始数据: capture panel2
├── 5.angkor_shotgun_finished/   # 原始数据: shotgun
├── 4.mcp_reshotgun/             # 原始数据: MCP proxy样本(⚠️不属于angkor 16样本, 分析时需排除)
└── angkor_robot_library.txt     # 样本元数据(Library_ID/Depth地层深度/Age年代/robot_sample_id映射)
```

---

## 二、原始数据

| 路径 | 内容 | 来源 |
|---|---|---|
| `5.angkor_shotgun_finished/data/reads/` | shotgun原始reads(`{Library_ID}.prefiltered.IRGSP1.mapped.fq`) | 已经过预筛选比对到IRGSP1，不需要重新提取 |
| `3.angkor_capture_panel1/data/reads/` | capture panel1原始reads(`*.bbduk.lowcomp_filtered.fq`) | 杂交捕获测序，已过bbduk低复杂度过滤 |
| `7_angor_capture_panel2/data/reads/` | capture panel2原始reads | 同上 |
| `4.mcp_reshotgun/data/reads/` | MCP reshotgun proxy样本 | **不属于本项目16个angkor古稻样本**，文件名格式`LV{数字}-LV{数字}-proxy...`(短横线)，跟真实panel文件`LV{数字}_RicePanel{1|2}...`(下划线)不同，务必在任何统计脚本里显式排除 |
| `angkor_robot_library.txt` | 样本元数据总表 | 含`Library_ID`→`robot_sample_id`映射、`Depth`(考古地层深度,cm)、`Age`(公历年代)，跨度约公元1160-1981年 |
| `asn720data/asn720.pop.{bed,bim,fam}` | 720份现代品种PLINK基因型 | 94974个SNP位点，⚠️与16个angkor样本的关系尚未最终确认(是否包含、是否capture panel来源)——基因型数据密度太低不用于PCA，但`.fam`的FID列(`OrA-OrF`群体标签)是`ecotype-pca-panel`分支目前找到的关键标签来源，详见该分支`docs/ECOTYPE_PCA_PANEL.md` |

---

## 三、参考基因组与索引(`db/`)

```
db/asian_rice_panel_index/
├── irgsp.fa[.amb/.ann/.bwt/.pac/.sa/.fai]   # IRGSP1.0主参考(最终定量比对用)
├── irgsp_bt2idx.*.bt2                        # irgsp.fa的Bowtie2索引(9格矩阵测试时补建)
├── irgsp.asm10.mmi                           # irgsp.fa的minimap2索引(野生稻组装比对用)
├── asian_rice_panel.fa[.amb/.ann/.bwt/.pac/.sa]  # 多物种鉴定panel(提取阶段用，非irgsp.fa!)
├── asian_rice_panel.fa.*.bt2l                # 上述panel的Bowtie2索引(large index格式,约10.4GB)
├── blastdb/ref.*                             # BLAST索引(低复杂度区域交叉验证用)
└── lowcomplexity_qc/
    └── lowcomplexity.sorted.bed              # dustmasker标记的全基因组低复杂度区域

db/gene/
├── msu7.gff3                        # MSU7官方注释
├── msu7_pseudomolecule.fna          # MSU7拟分子fasta(Liftoff源参考)
├── flower_gene.txt                  # 57个开花基因清单(Level/Group/MSU_id/RAPDB_id/Name)
├── flower_gene_ids.tsv              # MSU_id→Name精简对照
├── msu7.genes.bed                   # 全部MSU7基因坐标
├── flower_gene.sorted.bed           # 57基因严格边界坐标(chr01格式)
└── flower_gene.flank1kb.bed         # 57基因±1kb扩展坐标

db/16/                               # 【资源组A】NCBI datasets 16基因组+Liftoff
├── download16.sh / reorganize_rice_genomes.sh / 16_3k.csv
├── asian_rice_panel_download/       # 原始下载(嵌套目录)
└── asian_rice_panel/                # 扁平结构, 每个子目录含genome.fna+liftoff_from_msu7.gff3
    (genome1_IRGSP, genome4_N22, genome5_AZ, genome6_IR64, genome7_ARC,
    genome8_LM, genome9_LX, genome10_KYG, genome11_LIMA, genome12_NABO,
    genome13_PR106, genome14_KN, genome15_CM, genome16_GS,
    genome27_MH63, genome28_ZS97)

db/3k/                               # 【资源组B】3K Rice Genome Project数据
├── 16_3k/                           # 官方16参考基因组面板zip(与资源组A不同来源，未解压)
├── NB_bialSNP_pseudo_canonical_ALL.vcf.gz    # ⚠️性质待确认: 可能是过滤精简版VCF
├── CNVnator_Q10_goodRD_noCN1-3.vcf.gz
├── NB_{DEL,DUP,INS,INV}_mergesam_clustered.tar.gz
├── Nipponbare_indel.{bed,bim,fam}.gz
├── tmp/                             # ⚠️核心: 原始最完整3200万位点矩阵
│   ├── 3kall_snpposition_map.tsv    # SNP_INDEX/CHROMOSOME/POSITION/REFCALL
│   ├── 3kall_variety_map.tsv        # VARIETY_INDEX/NAME/IRIS_ID/...
│   └── Universe_matrix_geno_NB      # 3024行(品种)x3200万字符(SNP)基因型矩阵
├── sv_extracted/
│   └── NB_DEL_mergesam_clustered.txt  # 解压后的DEL记录
└── wild/
    └── {SampleID}.transfer.merge.chr.fasta  # 140+野生稻/近缘种染色体级组装基因组
                                                # ⚠️物种身份未确认，是否与Guo et al. 2025
                                                # pangenome论文的145组装同源待查，见
                                                # ORYZA_BESTHIT_HANDOFF.md第7.2节 /
                                                # RESEARCH_ROADMAP.md第2节C/第6节P0第2条

db/29M_3k/                           # 【ecotype-pca-panel分支用】3K RG 29mio biallelic SNP
├── NB_final_snp.bed.gz / .bim.gz / .fam.gz   # 3024份材料，PLINK格式，与irgsp.fa坐标系一致

db/6.7M_720/                         # 【ecotype-pca-panel分支用】720份样本，野生稻为主
├── asn720.6m.geno / .ind / .snp     # EIGENSTRAT格式(smartpca原生输入格式)

db/asian_rice_panel_index/all_wgs_asian_irgsp.acc2taxid  # 【besthit分支用】accession→taxid映射
```

**5个早期近缘品种参考基因组**(Ensembl Plants下载，早于资源组A/B):
```
/home/scratch/yinmt202607/outputs/rice_test_panel/
├── nipponbare_irgsp10/  azucena_rs1/  ir64_osir64rs1/  mh63_rs2/  n22_osn22rs2/
```

---

## 四、主流程处理产出(`results/02.irgsp/`)

```
results/02.irgsp/
├── 00.reads_bwa/                      # BWA提取的shotgun+capture reads(当前主线输入)
│
├── 01.mapping_bwa/final/              # 【当前主线最终BAM】
│   ├── {robot}.dedup.bam[.bai]        #   全量, 重复只标记(flag 0x400)未删除
│   └── {robot}.dedup.q30.bam[.bai]    #   已排除重复+MAPQ≥30, 下游分析用这份
│
├── 01.mapping/final/                  # 【历史Bowtie2版最终BAM】(早期流程, markdup加了-r真删除)
│
├── 02.gene_hits/                      # 旧Bowtie2版57基因命中统计
│   └── gene_hits_with_metadata.tsv    # 42条命中, 含read序列/mapq/cigar/Depth/Age
├── 02.gene_hits_bwa/                  # 【当前主线】57基因命中统计
│   ├── sample_summary_bwa.tsv
│   ├── gene_by_sample_matrix.tsv
│   └── gene_by_sample_matrix_hitgenes_only.tsv  # 41/57个有命中的基因
│
├── 04.3krgp_flowergenes/               # 3K RG与57基因交叉分析
│   ├── flower_genes_3krgp_snps.tsv     # 12741个落在57基因内的3K RG已知SNP
│   ├── our_reads_matching_3krgp_snps.tsv  # 6条古稻reads精确匹配已知SNP位点
│   └── per_snp_columns/snp_{idx}.txt   # 逐位点提取的3024品种基因型
│
├── 05.flank1kb/                       # ±1kb侧翼区域探索(结论: 未捕捉到调控区)
│   ├── bam/{robot}.flank1kb.bam
│   └── fasta/{robot}_{gene}.fa         # consensus序列(多为空/近全N)
│
├── 06.simulated_reads/                # NGSNGS模拟数据(16样本x16锚点x3档reads量级)
├── 06.simulated_mapping/final/        # 模拟数据的比对结果
│
├── readlen_dist_bwa/{robot}.length_cdf.txt   # NGSNGS -lf 输入
├── mapdamage_out_bwa/{robot}/                # NGSNGS -mf 来源 + 古DNA损伤验证
│   ├── 5pCtoT_freq.txt / 3pGtoA_freq.txt
│   ├── misincorporation.txt
│   └── ngsngs_mf.txt                  # 转换后的NGSNGS损伤格式(200行)
├── mapdamage_summary_bwa.tsv           # 16样本损伤强度汇总
├── damage_curves_all_samples.png       # 损伤曲线对比图
│
├── gene_hits_lowcomplexity_check_v2.tsv  # 低复杂度过滤结果(35.7%命中存疑)
└── compare_bowtie2_vs_bwa_final.tsv       # 16/16样本BWA提升3.3-3.6倍的证据

results/asian_rice_compare/           # Bowtie2 vs BWA提取阶段历史对比(compare_rice_read_extractors.sh产出)
├── bowtie2/{bam,fastq}/               # ⚠️混有14个MCP proxy样本, 使用前需过滤
└── bwa/{bam,fastq}/                   # 同上

results/07.wild_rice_alignment/        # 野生稻minimap2比对(进行中)
├── bam/{sample}.bam[.bai]
├── bam/{sample}.paf                   # 供paftools.js call变异用
└── sbatch_logs/

results/igv_package/                   # IGV可视化打包(build_igv_package.sh产出)
├── bam_all/    # 全量BAM(16样本)
├── bam_q30/    # 过滤版BAM(16样本)
├── ref/        # irgsp.fa+索引
└── annotation/ # flower_gene.sorted.bed + flank1kb.bed
```

---

## 五、参数测试目录(`tests/param_matrix_bt2_vs_bwa/`)

独立于正式pipeline的测试区，验证"3种提取方法x3种定量比对工具"9格矩阵：

```
tests/param_matrix_bt2_vs_bwa/
├── scripts/           # 01-11号脚本(建目录→提取→建索引→比对→汇总→可视化)
├── 00.extraction/     # bt2_old(历史软链接) / bt2_new(新生成) / bwa(历史软链接)
├── 01.reads_combined/ # 三种提取方法各自的shotgun+panel整理结果
│                      # ⚠️其中 01.reads_combined/bwa/ 是 codex/oryza-screen-merge
│                      #   分支的输入来源
├── 02.final_mapping/  # 9个组合的最终BAM(2历史软链接+7新生成)
├── dup_recompute/     # 补算历史①②组合真实dup信息
└── summary/           # 最终结果表+可视化PDF
```

**核心结论**：提取阶段用BWA是决定性因素，定量比对阶段BWA/Bowtie2新参数皆可接受。
详见`docs/09_extraction_mapping_matrix_final.md`。

---

## 六、Git仓库结构(`github.com/Inmpain/rice_adna_pipeline`)

⚠️ **本仓库有4个分支，以下只是main分支的结构**。另外3个分支
（`codex/oryza-screen-merge`、`codex/oryza-competitive-mapping`、
`codex/ecotype-pca-panel`）各自有独立的`docs/`和`scripts/`内容，
不会出现在main分支里，需要单独切换查看。四分支汇总见
`docs/REPO_OVERVIEW_STATUS.md`；项目级科研框架(证据阶梯/五条工作线/
Phase 0-4实施计划)见`docs/RESEARCH_ROADMAP.md`。

```
rice_adna_pipeline/ (main分支)
├── PROJECT_STATUS.md          # ⚠️已知问题: 曾被意外截断，只剩最后一次追加内容，需修复
├── docs/
│   ├── file_paths.md                              # 本文档
│   ├── REPO_OVERVIEW_STATUS.md                    # 四分支跨线现状汇总
│   ├── RESEARCH_ROADMAP.md                        # 证据阶梯+五条工作线+实施阶段(新增)
│   ├── decisions_log.md                           # 关键决策记录
│   ├── research_goals.md                          # 四层研究目标拆解
│   ├── GIT_USAGE.md                                # Git操作指南
│   ├── flank1kb_msa_exploration.md                # 低复杂度QC发现全过程
│   ├── 3krgp_integration_and_simulation_prep.md   # 3K RG/SV/NGSNGS环境搭建
│   ├── 09_extraction_mapping_matrix_final.md      # 9格矩阵测试完整记录
│   └── IGV_visualization_guide.md                 # IGV可视化包使用说明
├── scripts/server_originals/   # 9个真实生产脚本(唯一权威版本)
├── tools/lowcomplexity_qc/     # 可复用低复杂度QC工具
├── tests/param_matrix_bt2_vs_bwa/  # 9格矩阵测试(脚本+结果表+图，不含BAM)
└── results/                    # 各阶段统计表格快照(不含BAM等大文件)
```

**⚠️git仓库当前已知问题(2026-08-05核实)**：
- `PROJECT_STATUS.md`被截断，只剩最新一次追加内容，原始核心内容丢失
- 其余文档(`decisions_log.md`等)是否完整、`IGV_visualization_guide.md`是否已提交，
  尚未逐一核实，建议接手人先跑一次`find . -type f`+`wc -l`核对，见下节

---

## 七、核对/维护建议

```bash
cd ~/rice_adna_pipeline

# 1. 看完整文件树
find . -type f -not -path "./.git/*" | sort

# 2. 核对每个关键文档的行数(判断是否被截断)
for f in PROJECT_STATUS.md docs/*.md; do
    echo "$f: $(wc -l < "$f") 行"
done

# 3. 看提交历史，确认哪些内容真的进去了
git log --oneline --stat

# 4. 确认没有未提交的改动堆积在本地
git status
```

**关于`PROJECT_STATUS.md`修复**：

```bash
git log --oneline -- PROJECT_STATUS.md
# 找到内容完整的那个commit hash(损坏之前的版本)
git show <commit_hash>:PROJECT_STATUS.md > /tmp/old_version.md
# 人工合并旧版本(项目总结/待办/文档索引)和新版本(最新进展)，重新写一份完整的
```

---

## 八、给新接手人/新对话的快速定位指南

| 我想知道... | 去看... |
|---|---|
| **我们到底要证明什么、现在能证明到哪一级**（先看这个） | `docs/RESEARCH_ROADMAP.md`(main分支，证据阶梯+五条工作线+实施阶段) |
| **四个分支现在各自做到哪一步** | `docs/REPO_OVERVIEW_STATUS.md`(main分支，跨分支横向汇总) |
| 项目现在做到哪一步、下一步该干什么 | `PROJECT_STATUS.md`(⚠️当前已损坏，需先修复) |
| 为什么选了BWA而不是Bowtie2 | `docs/decisions_log.md` + `docs/09_extraction_mapping_matrix_final.md` |
| 57基因命中数据为什么只有约2/3可信 | `docs/flank1kb_msa_exploration.md` |
| 3K RG数据怎么用、SV断点怎么查的 | `docs/3krgp_integration_and_simulation_prep.md` |
| 某个具体文件在服务器哪个路径 | 本文档(`docs/file_paths.md`) |
| 怎么用git提交新东西 | `docs/GIT_USAGE.md` |
| 怎么在IGV里看数据 | `docs/IGV_visualization_guide.md` |
| 某个脚本具体是干什么的 | `scripts/server_originals/`(真实脚本) 或
  `tests/param_matrix_bt2_vs_bwa/scripts/`(参数测试脚本) |
| Oryza best-hit过滤(竞争性比对+古DNA损伤校正)进展 | 切到`codex/oryza-competitive-mapping`分支，读`docs/ORYZA_BESTHIT_HANDOFF.md` |
| 旱稻/水稻生态型PCA判定进展 | 切到`codex/ecotype-pca-panel`分支，读`docs/ECOTYPE_PCA_PANEL.md` |
| 候选Oryza FASTQ怎么合并的 | 切到`codex/oryza-screen-merge`分支，读`oryza_screen_merge/README.md` |
