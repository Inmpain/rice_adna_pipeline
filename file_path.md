# 完整路径清单 (v5 — 最终权威版，供工作交接使用)

> 本文档是整个项目的路径地图，每个目录都标注了**是什么、由哪个脚本产出、
> 什么场景下需要用它**。新接手的人应该先读这份文档，再去看
> `PROJECT_STATUS.md`(项目当前进度)和`decisions_log.md`(关键决策原因)。

最后更新: 2026-08-16

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
| `4.mcp_reshotgun/data/reads/` | MCP reshotgun proxy样本 | **不属于本项目16个angkor古稻样本**，文件名格式`LV{数字}-LV{数字}-proxy...`(短横线)，跟真实panel文件`LV{数字}_RicePanel{1\|2}...`(下划线)不同，务必在任何统计脚本里显式排除 |
| `angkor_robot_library.txt` | 样本元数据总表 | 含`Library_ID`→`robot_sample_id`映射、`Depth`(考古地层深度,cm)、`Age`(公历年代)，跨度约公元1160-1981年 |
| `asn720data/asn720.pop.{bed,bim,fam}` | 720份现代品种PLINK基因型 | 94974个SNP位点，⚠️与16个angkor样本的关系尚未最终确认(是否包含、是否capture panel来源) |

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
├── 02.final_mapping/  # 9个组合的最终BAM(2历史软链接+7新生成)
├── dup_recompute/     # 补算历史①②组合真实dup信息
└── summary/           # 最终结果表+可视化PDF
```

**核心结论**：提取阶段用BWA是决定性因素，定量比对阶段BWA/Bowtie2新参数皆可接受。
详见`docs/09_extraction_mapping_matrix_final.md`。

---

## 六、Git仓库结构(`github.com/Inmpain/rice_adna_pipeline`)

```
rice_adna_pipeline/
├── PROJECT_STATUS.md          # ⚠️已知问题: 曾被意外截断，只剩最后一次追加内容，需修复
├── docs/
│   ├── file_paths.md                              # 本文档
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
| 项目现在做到哪一步、下一步该干什么 | `PROJECT_STATUS.md`(⚠️当前已损坏，需先修复) |
| 为什么选了BWA而不是Bowtie2 | `docs/decisions_log.md` + `docs/09_extraction_mapping_matrix_final.md` |
| 57基因命中数据为什么只有约2/3可信 | `docs/flank1kb_msa_exploration.md` |
| 3K RG数据怎么用、SV断点怎么查的 | `docs/3krgp_integration_and_simulation_prep.md` |
| 某个具体文件在服务器哪个路径 | 本文档(`docs/file_paths.md`) |
| 怎么用git提交新东西 | `docs/GIT_USAGE.md` |
| 怎么在IGV里看数据 | `docs/IGV_visualization_guide.md` |
| 某个脚本具体是干什么的 | `scripts/server_originals/`(真实脚本) 或
  `tests/param_matrix_bt2_vs_bwa/scripts/`(参数测试脚本) |

---

## 九、`codex/ecotype-pca-panel`分支路径(2026-08-14新增)

**跟上面八节描述的主流程(BWA提取/57基因命中)是平行的独立分析线**，
详见该分支`docs/ECOTYPE_PCA_PANEL.md`的📍handoff段落。这里只补充本文档
之前完全没覆盖的路径。

**⚠️没有软链接**：本节下面全部路径都是真实目录/文件，直接建在
`/home/scratch/yinmt202607/`下——**不经过`/itp`软链接层**(不同于本项目
其他部分"新产出放itp、软链接回scratch"的惯例，见长期记忆
`itp-storage-symlink-convention`)。这一点是根据本次会话里用户实际粘贴
的终端输出推断的，没有专门跑`ls -la`逐个核实每一条，如果发现某条其实
是软链接，以服务器真实情况为准，回来更新这一节。

```
/home/scratch/yinmt202607/
├── db/
│   ├── 29M_3k/                       # PCA-A：3K RGP栽培稻，3024份×~29.6M biallelic SNP
│   │   ├── NB_final_snp.{ind,snp,eigenstratgeno}       # 原始(convertf转换产出)
│   │   ├── NB_final_snp.labeled.ind                    # 群体标签(IND/AUS/ARO/TRJ/TEJ/...)
│   │   ├── NB_final_snp.label_report.tsv
│   │   ├── NB_final_snp.filtered.{ind,eigenstratgeno}  # UNK剔除后，下游一律用这份
│   │   └── references/rice_line_metadata_20141029.xlsx # 3K RGP官方元数据(来自main分支docs/references/)
│   ├── 6.7M_720/                     # PCA-B：野生稻为主，720份×~6.7M SNP
│   │   ├── asn720.6m.{ind,snp,geno,geno.gz}            # 原始(作者私发的加密版，来源/处理流程不透明)
│   │   ├── asn720.6m.labeled.ind / label_report.tsv    # 群体标签(来自asn720data的OrA-OrF)
│   │   └── asn720.6m.filtered.{ind,geno}               # UNK剔除后
│   ├── paper1/                       # PCA-C：Civáň 2019桥接面板，1056份×2,365,188 SNP
│   │   ├── civan_snp.{ind,snp,eigenstratgeno}          # 原始(VCF→convertf转换产出)
│   │   ├── Table_S1.csv / Table_S2.csv                 # 论文原始元数据(样本×祖源/SRA来源)
│   │   ├── civan_snp.labeled.ind / label_report.tsv
│   │   └── civan_snp.filtered.{ind,eigenstratgeno}     # UNK剔除后
│   └── wild_rice_pangenome_README.txt                  # 内容为空，未解决"OrA-OrF"定义来源问题
├── gene/
│   ├── scripts/                      # `scripts/ecotype_pca/`的部署位置(curl下载到这里执行)
│   └── results/ecotype_pca/
│       ├── panel_overlap/            # summarize_panel_overlap.py产出(位点重叠+QC统计)
│       ├── bam_irgsp/                # besthit_oryza→IRGSP全属比对BAM，mapping_summary.tsv
│       ├── loo_smoke/                # Civáň leave-one-out smoke test完整产出(唯一真正跑过的PCA)
│       └── pca_runs/                 # v1的16×3 first-look已产出；无MAF/LD、每样本独立坐标，不是final
```

**Git仓库(`codex/ecotype-pca-panel`分支)新增内容**：
```
scripts/ecotype_pca/
├── summarize_panel_overlap.py                  # 位点重叠+QC统计(pysam pileup)
├── build_{29m3k,720,civan}_population_labels.py # 三个panel的群体标签匹配
├── filter_panel_by_label.py                    # 按标签物理剔除样本(UNK剔除用这个)
├── build_sample_panel_subset.py                # 样本专属子集(⚠️设计缺陷，见QC_DESIGN第1节)
├── simulate_leaveoneout_projection.py           # leave-one-out正对照模拟
├── merge_ancient_into_panel.py                  # 古样本合并进panel
├── run_sample_panel_pca.sh                      # ①-④四步流水线封装+批量入口
├── summarize_projection_distances.py            # .evec最近现代群体排名
├── civan_domesticated_reference_labels.txt      # Civáň panel axis-building标签清单(poplist bug修复配套)
├── pseudo_haploid_call.py / map_besthit_to_irgsp.sh  # 更早期已有脚本
└── par.MERGE / run_convert_merge.sh             # ⚠️已废弃，仅存档mergeit调试历史，不要用

docs/
├── ECOTYPE_PCA_PANEL.md              # 本分支的"为什么"+📍handoff入口
├── ECOTYPE_PCA_EXECUTION_PLAN.md     # 10步执行计划+完成状态
├── ECOTYPE_PCA_PHASE0_COMMANDS.md    # Phase 0(位点重叠census)服务器命令
├── ECOTYPE_PCA_PHASE1_COMMANDS.md    # Phase 1(标签匹配→子集PCA→LOO)服务器命令
├── ECOTYPE_PCA_PANEL_QC_DESIGN.md    # MAF/LD pruning设计+GPT review记录+待办顺序
└── references/3k_rice_genomes_project/  # 3K RGP官方文档(从main分支同步)
```

对应的完整"为什么/怎么用"说明见`docs/ECOTYPE_PCA_PANEL.md`📍段落，不在
本文档重复。

### 9.1 PCA v2顺序执行与debug状态目录（2026-08-16新增）

**仓库内权威路径**：

```text
scripts/ecotype_pca_v2/
├── bootstrap_ecotype_pca_v2.sh       # 按40位commit下载不可变仓库快照
├── config/ecotype_pca_v2.yaml        # 全阶段统一参数源；变更会使receipt失效
└── workflow/
    ├── workflow.json                 # 唯一阶段顺序/门禁定义
    ├── ecotype_pca_workflow.py       # 状态机控制器
    ├── collect_server_evidence.py    # 720/718与BAM单次flagstat只读证据
    ├── runners/                      # 当前已开放阶段的精确runner
    └── tests/                        # 顺序、digest、stale receipt回归测试

docs/ECOTYPE_PCA_V2_WORKFLOW.md       # 功能/运维说明
docs/ECOTYPE_PCA_V2_SPEC.md           # 参数/统计细节说明（与功能说明分层）
```

**建议服务器路径**：

```text
/home/scratch/yinmt202607/gene/workflow_sources/
└── rice_adna_pipeline-<FULL_COMMIT>/ # 每次debug修复安装一个新commit，不覆盖旧版

/home/scratch/yinmt202607/gene/results/ecotype_pca_v2/workflow_state/
├── receipts/                         # 成功阶段的内容寻址receipt
├── attempts/<stage>/<timestamp>/     # 每次独立尝试；失败也保留
└── debug_bundles/                    # 发回本地debug的小型tar.gz
```

源码版本可以更换，`workflow_state`应复用。控制器会根据config、阶段定义、
tracked scripts和上游receipt哈希判断旧结果仍有效还是`STALE`，不得手改receipt。
stage 00可在login节点运行；stage 10/20必须由SLURM执行。
