# Ecotype PCA v2：冻结统计设计与执行规范

## 最高优先级指令

```text
STATISTICAL DESIGN IS FROZEN.

Your role is implementation only.

You may write scripts, configuration files, validation checks, logs,
summaries, smoke tests, and documentation.

You MUST NOT redesign the analysis, change parameters, substitute
methods, redefine reference populations, select parameters based on
the appearance of PCA plots, or make undocumented analytical decisions.

All numerical parameters, sample roles, marker definitions, analysis
tracks, execution phases, validation rules, and confidence thresholds
specified below are immutable.

If any required input, metadata, software capability, sample label,
bait BED, population mapping, or file relationship is missing,
ambiguous, inconsistent, or technically impossible:

STOP.

Report the exact problem and wait for user instruction.

Do not choose a fallback.
Do not infer a replacement.
Do not silently skip the affected step.
Do not change any parameter.
```

---

# 一、分析目标

建立由现代参考样本定义的固定 PCA 坐标系，然后将古代样本投影到固定坐标系中。

核心原则：

```text
现代参考样本
→ 定义 marker universe
→ 完成 missingness、MAF、LD/thinning
→ 冻结 fixed marker set
→ 建立 fixed PCA axes
→ 投影 wild、ADM 和 ancient
```

古样本不得参与：

* marker selection；
* site missingness；
* MAF；
* LD pruning；
* PCA axes 构建。

古样本之间不需要覆盖相同 SNP。

禁止计算：

```text
Ancient A ∩ Ancient B ∩ Ancient C
```

并把交集作为正式 PCA marker set。

固定的是：

```text
reference marker universe
+
PCA loadings
```

不是古样本的 non-missing SNP。

古样本未覆盖的 fixed markers 必须编码为：

```text
9 = missing
```

增加第二个或更多古样本时，不得重新计算 marker set 或 PCA axes。

---

# 二、旧方法的地位

旧的 sample-specific 方法：

```text
ancient callable SNP
→ modern panel 裁成相同 SNP
→ 重新计算 PCA
```

可以保留，但必须命名为：

```text
MATCHED_MARKER_DIAGNOSTIC_PCA
```

它只能回答：

> 在该古样本实际覆盖的 marker subset 上，现代结构和该古样本的 affinity 是什么？

严禁将其用于：

* 跨古样本比较 PC 坐标；
* 时间轨迹分析；
* production population assignment；
* 替代 fixed-reference PCA。

---

# 三、Shotgun、Capture、TV、ALL

## 3.1 Shotgun和Capture必须分开

Shotgun marker universe：

```text
panel SNP
→ modern QC
→ MAF
→ LD/thinning
→ shotgun fixed markers
```

Capture marker universe：

```text
panel SNP
∩ rice capture bait BED
→ modern QC
→ MAF
→ LD/thinning
→ capture fixed markers
```

不得在与 bait BED 求交之前进行 capture reference marker selection。

每个 panel 原则上建立四套独立坐标系：

```text
SHOT
```
你现在负责执行 `Ecotype PCA v2`。以下统计设计和全部参数已经冻结。

# 一、权限边界

你的职责仅限于：

* 检查现有仓库、输入文件和软件环境；
* 按本指令编写模块化脚本和配置；
* 运行输入审计、smoke test 和规定的正式分析；
* 生成日志、manifest、QC、图表和汇总结果；
* 按规定顺序提交 git commit；
* 如实报告错误、缺失输入、软件不兼容和资源限制。

你无权：

* 修改任何 MAF、missingness、LD、MAPQ、BaseQ、fold、置信度或位点数量阈值；
* 更换 reference population；
* 改变哪些群体定义 PCA axes；
* 根据 PCA 图是否“好看”选择参数；
* 为增加 ancient SNP 数取消过滤；
* 求所有 ancient 的共同 SNP intersection；
* 为每个 ancient 重新建立正式 PCA axes；
* 混合 shotgun 和 capture 坐标系；
* 把 TV 降为次要分析，或因 ALL 位点更多而自动提升 ALL；
* 把 diagnostic 分析冒充 production 分析；
* 根据历史预期挑选结果；
* 擅自增加 sensitivity 参数；
* 遇到失败后换用未授权方法；
* 自行解释 OrA–OrF 的地理或祖先含义；
* 未经明确授权提交大规模 cluster jobs。

如果严格执行遇到障碍：

1. 立即停止受影响的步骤；
2. 保留日志和中间文件；
3. 明确报告失败位置、错误信息和影响；
4. 提出问题等待用户决定；
5. 不得自行更改参数或统计方法绕过问题。

不要重新讨论、优化或 redesign 本方案。你的任务是执行，不是制定统计策略。

旧目录必须保持不变：

```text
results/ecotype_pca/
```

所有新结果写入：

```text
results/ecotype_pca_v2/
```

# 二、核心统计原则

正式分析必须先用现代参考样本建立并冻结：

```text
modern reference samples
→ reference marker universe
→ missingness filtering
→ MAF filtering
→ LD pruning / fixed thinning
→ fixed SNP list
→ fixed PCA axes and loadings
```

然后再处理 ancient：

```text
ancient BAM
→ call whatever fixed markers are covered
→ uncovered fixed markers coded as 9
→ project onto the fixed modern PCA axes
```

必须固定的是：

```text
reference marker universe + PCA loadings
```

不是：

```text
ancient non-missing SNP set
```

不同 ancient 可以覆盖完全不同的 SNP。禁止计算：

```text
Ancient A SNP ∩ Ancient B SNP ∩ Ancient C SNP
```

并将其作为正式 PCA marker set。

新增第二个或更多 ancient 时，不得重新计算 reference PCA，不得改变 fixed marker list 或 PCA axes。

# 三、shotgun 与 capture 分开

Shotgun marker universe：

```text
panel SNP
→ modern QC
→ TV/ALL
→ MAF
→ LD/thinning
→ fixed marker set
```

Capture marker universe：

```text
panel SNP
∩ rice capture bait BED
→ modern QC
→ TV/ALL
→ MAF
→ LD/thinning
→ fixed marker set
```

每个 panel 原则上建立四套独立坐标系：

```text
SHOTGUN.TV
SHOTGUN.ALL
CAPTURE.TV
CAPTURE.ALL
```

不同坐标系的 PC 数值不得直接比较。

# 四、TV与ALL

TV 是 primary ancient-DNA track，仅保留：

```text
A/C
A/T
C/G
G/T
```

排除 transitions：

```text
A/G
C/T
```

ALL 是 secondary track，保留全部 SNP。

不得因为 ALL 位点数更多而宣布 ALL 更可靠。ALL 能否提升为主要证据，只能由用户结合 UDG 信息和 damage profile 决定。

# 五、三个panel的固定用途

## Panel A：29M_3K

用途：

```text
cultivated group assignment
```

PCA axis builders 只能是：

```text
IND
AUS
ARO
TRJ
TEJ
```

以下只能 projection：

```text
ADM
Ancient
```

UNK 不恢复。

MAF、missingness 和 LD 只能在五个 axis-builder groups 中计算。Ancient 和 ADM 不得参与。

Primary 参数：

```text
biallelic only
chromosomes 1–12 only
unique position/allele identifier
site missingness <= 0.05
MAF >= 0.01
LD window = 100kb
LD r² = 0.20
```

PLINK2 使用：

```bash
--indep-pairwise 100kb 0.2
```

不得使用：

```text
100kb 10 0.2
```

Sensitivity 只能有：

```text
S1: 50kb,  r²=0.20
S2: 200kb, r²=0.20
S3: 100kb, r²=0.10
S4: primary LD + MAF=0.05
```

不得增加其他组合。

## Panel B：6.7M_720

用途：

```text
wild / feral / cultivated-related structure
```

该文件必须标记为：

```text
AUTHOR-PROVIDED HIGH-DENSITY VERSION
FILTER PROVENANCE UNKNOWN
NOT THE EXACT WANG2017 PCA MATRIX
```

不得假设作者已经完成 MAF 或 LD 过滤。

先执行 audit，不删数据、不跑 PCA。输出：

```text
720.audit.samples.tsv
720.audit.maf.tsv
720.audit.missingness.tsv
720.audit.spacing.tsv
720.audit.summary.txt
```

Audit 至少包括：

```text
N samples
N SNP
per-pop sample counts

MAF bins:
0
0–0.001
0.001–0.005
0.005–0.01
0.01–0.05
0.05–0.10
>0.10

site missingness bins:
0–0.01
0.01–0.05
0.05–0.10
0.10–0.20
>0.20

sample missingness

adjacent SNP spacing:
median
P10
P25
P75
P90
fraction <1kb
fraction <5kb
fraction <10kb
```

Panel B 固定运行两条路线。

### B-primary

PCA axis builders：

```text
all technically valid 720 modern samples
```

只能因 malformed sample 或全部 genotype missing 这类技术失败排除样本。不得因群体位置奇怪而排除。

参数：

```text
site missingness <= 0.10
MAF >= 0.01
LD window = 100kb
LD r² = 0.20
```

Sensitivity 只能有：

```text
50kb / 0.20
200kb / 0.20
100kb / 0.10
```

### B-paperlike

先应用：

```text
site missingness <= 0.10
MAF >= 0.01
```

再按 non-overlapping 5000-bp window 随机选择最多一个 SNP：

```text
window = 5000 bp
seed = 20260814
```

名称必须是：

```text
paperlike_5kb
```

禁止命名为：

```text
Wang2017_exact
```

如果 B-primary 和 paperlike_5kb 拓扑明显不同，只报告：

```text
PANEL_B_MARKER_SENSITIVITY = TRUE
```

不得自行决定哪张图更正确。

OrA–OrF 只能作为 genetic cluster labels。metadata 未核实前，不得赋予地理、祖先或生态含义。

## Panel C：Civáň

这是第一个 v2 prototype。

PCA axes 只能由595个 domesticated accessions 建立：

```text
indica
aus
aromatic
japonica
japonica_(temperate)
japonica_(tropical)
```

以下全部只能 projection，不得进入 `poplistname`：

```text
O._rufipogon
O._barthii
O._glaberrima
O._longistaminata
O._meridionalis
Ancient
```

必须核实 domesticated axis builders 总数是否为595。若不是，停止并报告，不得自行修补类别。

MAF、missingness 和 LD 只能在595个 domesticated samples 中计算。

Primary 参数：

```text
site missingness <= 0.05
MAF >= 0.01
LD window = 100kb
LD r² = 0.20
```

Sensitivity 只能有：

```text
50kb / 0.20
200kb / 0.20
100kb / 0.10
```

另外保留：

```text
CIVAN_PAPERLIKE_UNPRUNED
```

其规则是：

```text
使用原有 Civáň SNP matrix
595 domesticated define axes
wild + ancient projected
不新增 MAF/LD pruning
```

它只属于：

```text
diagnostic / literature-comparison
```

正式 ancient assignment 始终以 LD-pruned primary 为主。不得因为 unpruned 图更符合预期而选择 unpruned 作为正式答案。

# 六、fixed marker冻结

Ancient calling 之前，每个分析坐标系必须先生成：

```text
*.fixed.snplist
*.marker_manifest.tsv
```

Manifest 至少包含：

```text
panel
library_type
track
reference_samples_n
raw_snps
bait_overlap_snps
after_TV_ALL
after_site_missingness
after_MAF
after_LD_or_thinning
parameters
md5
```

一旦 `fixed.snplist` 生成：

* Ancient 不得增加位点；
* Ancient 不得删除位点；
* Ancient 不得触发重新选点；
* 未覆盖位置必须保留并编码为9。

Capture 必须先生成：

```text
panel.capture_compatible.snp
```

并报告：

```text
panel total SNP
SNP inside bait BED
capture-compatible ALL SNP
capture-compatible TV SNP
```

# 七、ancient calling固定参数

实现：

```text
10_call_ancient_fixed_markers.py
```

输入：

```text
BAM
fixed SNP list
track TV/ALL
sample
```

固定参数：

```text
MAPQ >= 30
BaseQ >= 30
exclude secondary alignments
exclude supplementary alignments
exclude QC-fail reads
exclude PCR duplicates
```

Pseudo-haploid calling：

```text
randomly choose one surviving read allele
```

随机种子必须稳定、可重复。先由下列内容确定 run seed：

```text
sample + panel
```

再由下列内容确定每个位点的 site seed：

```text
run seed + contig + position
```

`track` **不得**进入随机种子。TV 与 ALL 在共同的 transversion 位点必须
从同一批合格 reads 中抽到同一个 allele；否则 TV/ALL 差异会混入随机抽样
差异，不能解释为 transition 信号。panel SNP 行顺序和其他位点是否被访问也
不得改变该位点的抽样结果。

Genotype 编码：

```text
REF = 2
ALT = 0
no usable allele = 9
```

这里遵循 EIGENSTRAT 官方定义：genotype 数字是 reference allele copy 数。
因此 pseudo-haploid REF 为 2、ALT 为 0；不得按 VCF ALT dosage 习惯写反。
权威格式定义见DReichLab/EIG的
[`CONVERTF/README`](https://github.com/DReichLab/EIG/blob/master/CONVERTF/README#L62-L72)。

第三等位或不匹配等位：

```text
allele_mismatch
→ code as missing
→ report count
```

# 八、smartpca固定参数

所有 production smartpca 使用：

```text
lsqproject: YES
numoutevec: 10
numoutlieriter: 0
numchrom: 12
numthreads: 2
```

`poplistname` 必须严格对应：

Panel A：

```text
IND
AUS
ARO
TRJ
TEJ
```

Panel B：

```text
全部技术上合法的 modern population labels
```

Panel C：

```text
indica
aus
aromatic
japonica
japonica_(temperate)
japonica_(tropical)
```

Ancient 永远不得进入 `poplistname`。

# 九、callability与overlap报告

必须生成：

```text
ancient_callability.tsv
```

字段：

```text
sample
age
depth
panel
library_type
track
fixed_marker_n
callable_n
callable_fraction
information_flag
```

固定 flags：

```text
callable_n < 200       → VERY_LOW
200–499                → LOW
500–1999               → MODERATE
>=2000                 → HIGHER
```

这些是 workflow information flags，不得解释为普适科学阈值。

如果 LV7008416379 的 Civáň TV 仍为147 SNP，必须标记：

```text
VERY_LOW
```

还必须生成：

```text
callable_count_matrix.tsv
jaccard_matrix.tsv
```

Overlap matrix 只作描述，不得据此改做 ancient-intersection PCA。

# 十、PCA输出与assignment

不能只输出 PC1–PC2。必须输出：

```text
PC1–PC2
PC1–PC3
PC2–PC3
```

还要在以下维度计算距离：

```text
PC1–3
PC1–5
PC1–10
```

每个 ancient 至少计算：

```text
population centroid distance
nearest 20 modern individuals
nearest 50 modern individuals
```

Panel C 至少比较：

```text
INDICA
AUS
AROMATIC
JAPONICA_COMBINED
```

图上可以分别展示：

```text
japonica
japonica_(temperate)
japonica_(tropical)
```

但 broad classification 必须同时提供：

```text
JAPONICA_COMBINED
```

# 十一、exact-mask validation

每个 ancient 都必须进行 exact-mask validation。

例如 ancient 只有147个 callable TV SNP，就用完全相同的147个位点测试已知 modern samples：

```text
known modern sample
→ mask to exact ancient loci
→ pseudo-haploidize
→ project
→ classify
```

正式 validation 必须采用：

```text
stratified 5-fold
```

每个 fold：

```text
80% modern
→ define PCA and loadings

20% held-out modern
→ mask to exact ancient loci
→ pseudo-haploidize
→ project
→ classify
```

禁止将参与建立 PCA 的同一个 modern individual mask 后投回自己的 PCA，避免 data leakage。

Panel C broad classes：

```text
INDICA
AUS
AROMATIC
JAPONICA
```

Panel A classes：

```text
IND
AUS
ARO
TRJ
TEJ
```

如果任何 class：

```text
n < 10
```

标记：

```text
INSUFFICIENT_REFERENCE_N
```

停止该 class 的5-fold metric，不得临时换用其他验证方法。

# 十二、最终confidence规则

Ancient assignment 必须首先满足：

```text
nearest centroid in PC1–5
```

与：

```text
top20 neighbour majority
```

指向同一 broad group。

然后结合 exact-mask validation。

固定规则：

```text
HIGH:
precision >= 0.80
recall >= 0.80
and two assignment methods agree

MODERATE:
precision >= 0.60
recall >= 0.60
and two assignment methods agree

LOW / UNRESOLVED:
all other cases
```

例如 ancient 最近群体为 ARO，但 exact-mask validation 显示 aromatic 与 japonica 无法可靠区分，正式结果只能写：

```text
AROMATIC/JAPONICA-SIDE
UNRESOLVED
```

不得正式标记为 aromatic。

Sensitivity 汇总只允许报告：

```text
stable_3_of_3
stable_2_of_3
unstable
```

不得从多个参数结果中挑选最符合预期的结果。

# 十三、旧方法的地位

保留原有 sample-specific 方法：

```text
ancient callable SNP
→ subset modern panel to the same SNP
→ recompute PCA
```

但必须改名：

```text
MATCHED_MARKER_DIAGNOSTIC_PCA
```

它只能回答：

```text
在该 ancient 实际覆盖的 marker subset 上，
modern structure 与 ancient affinity 如何？
```

禁止用它：

* 跨 ancient 比较 PC coordinates；
* 绘制时间 trajectory；
* 作为 production population assignment；
* 替代 fixed-reference projection。

# 十四、唯一配置来源

所有参数只能来自一个版本控制下的 YAML config。

不得把参数散落硬编码到不同脚本中。脚本读取 config 后，必须把实际参数写入日志和 manifest。

配置至少包括：

```yaml
version: ecotype_pca_v2

ancient:
  mapq: 30
  baseq: 30
  primary_track: TV
  secondary_track: ALL

pca:
  num_pcs: 10
  numoutlieriter: 0
  numchrom: 12
  numthreads: 2
  lsqproject: true

panel_A_3k:
  axis_labels:
    - IND
    - AUS
    - ARO
    - TRJ
    - TEJ
  project_labels:
    - ADM
  geno: 0.05
  maf: 0.01
  ld_window_kb: 100
  ld_r2: 0.20

panel_B_720:
  axis_mode: all_modern
  geno: 0.10
  maf: 0.01
  ld_primary:
    window_kb: 100
    r2: 0.20
  paperlike_5kb:
    window_bp: 5000
    seed: 20260814

panel_C_civan:
  axis_labels:
    - indica
    - aus
    - aromatic
    - japonica
    - japonica_(temperate)
    - japonica_(tropical)
  project_all_other_modern: true
  geno: 0.05
  maf: 0.01
  ld_window_kb: 100
  ld_r2: 0.20

sensitivity:
  - window_kb: 50
    r2: 0.20
  - window_kb: 200
    r2: 0.20
  - window_kb: 100
    r2: 0.10

validation:
  folds: 5
  confidence_high:
    precision: 0.80
    recall: 0.80
  confidence_moderate:
    precision: 0.60
    recall: 0.60

information_flags:
  very_low: 200
  low: 500
  moderate: 2000
```

Panel A 的 `MAF=0.05` sensitivity 也必须明确写入 config，不得在脚本内硬编码。

# 十五、脚本结构

建立：

```text
scripts/ecotype_pca_v2/
```

脚本按下列职责拆分：

```text
00_validate_inputs.py
01_make_panel_manifest.py
02_convert_eigenstrat_for_plink.sh
03_audit_panel.py
04_audit_720_ld.py
05_intersect_panel_baits.py
06_build_reference_sample_set.py
07_make_fixed_markers.sh
08_make_5kb_thinned_markers.py
09_export_fixed_reference_eigenstrat.py
10_call_ancient_fixed_markers.py
11_build_ancient_callability.py
12_build_ancient_overlap_matrix.py
13_merge_ancients_fixed_panel.py
14_run_fixed_smartpca.sh
15_pca_qc.py
16_projection_summary.py
17_exact_mask_validation.py
18_validation_metrics.py
19_survey_ancient_coverage.py
20_filter_coverage_sites_to_transversions.py
21_extract_fixed_snplist.py
22_classify_scientific_projection.py
```

不得用一个大型单体 shell script 同时承担全部逻辑。

每个脚本必须：

* 支持 `--help`；
* 验证输入；
* 失败时返回非零状态；
* 记录参数和软件版本；
* 不静默跳过错误；
* 使用稳定排序和可重复随机种子；
* 不覆盖已有结果，除非显式指定安全的 overwrite 参数。

# 十六、执行顺序

必须严格按阶段执行。

## Phase 1：基础设施与audit

完成：

```text
config
input validation
format checks
panel manifests
MAF audit
missingness audit
720 spacing/LD audit
bait overlap audit
small smoke tests
```

这一阶段不运行正式 ancient PCA。

提交：

```text
1. ecotype-pca-v2 manifests and audit
2. fixed reference marker framework
3. capture bait intersection
```

## Phase 2：Civáň prototype

只使用：

```text
LV7008416379
```

运行：

```text
Civáň panel
→ verify 595 domesticated axis builders
→ define capture/shotgun marker universe
→ split TV/ALL
→ reference-only missingness and MAF
→ LD pruning
→ freeze fixed marker list
→ build fixed modern PCA
→ project wild
→ call and project LV7008416379
→ exact-mask validation
→ sensitivity summary
```

提交：

```text
4. Civan fixed reference PCA
5. fixed-marker ancient projection
6. exact-mask validation
```

只有 prototype 验收全部通过后才能进入 Phase 3。

## Phase 3：3K

构建一次 fixed reference，先只测试：

```text
LV7008416379
```

提交：

```text
7. 3K fixed reference PCA
```

## Phase 4：720

执行：

```text
LD-primary
paperlike_5kb
topology comparison
marker sensitivity report
```

提交：

```text
8. 720 audit and two marker routes
```

## Phase 5：全部ancient

只有三个 reference frameworks 全部冻结后，才批量 projection：

```text
1160
1239
1245
1265
1281
1314
1334
1372
1392
1414
1434
1439
1467
1555
1709
1981
```

以及已经测试的：

```text
LV7008416379
```

不得因加入这些样本重新选择 reference markers 或重算 PCA axes。

提交：

```text
9. multi-ancient batch projection
10. cross-panel summary
```

# 十七、Civáň prototype验收条件

以下项目必须全部通过：

```text
[ ] exactly 595 domesticated samples define axes
[ ] 461 wild samples are projection only
[ ] ancient is projection only
[ ] MAF is calculated from 595 domesticated samples only
[ ] missingness is calculated from 595 domesticated samples only
[ ] LD is calculated from 595 domesticated samples only
[ ] capture panel is intersected with bait BED first
[ ] fixed marker list exists before ancient calling
[ ] fixed marker manifest and md5 exist
[ ] uncovered ancient loci remain coded as 9
[ ] adding another ancient does not redefine PCA axes
[ ] sample-specific PCA is diagnostic only
[ ] PC1–PC2, PC1–PC3 and PC2–PC3 are generated
[ ] PC1–3, PC1–5 and PC1–10 distances are generated
[ ] LV7008416379 callable count is reported
[ ] 147 TV calls, if unchanged, are flagged VERY_LOW
[ ] exact-mask stratified 5-fold validation is completed
[ ] no automatic aromatic assignment is made without validation support
[ ] primary and sensitivity stability is reported
```

任何一项失败，都不得继续后续 panel 或全部 ancient 批处理。

# 十八、每阶段汇报格式

每个阶段结束后，只汇报事实：

```text
1. completed steps
2. files created
3. exact commands executed
4. input and output counts
5. parameters read from config
6. checks passed
7. warnings and failures
8. git commit hash
9. whether the next phase is authorized by the acceptance gate
```

不得在没有验证支持时给 ancient 强行分类。

现在开始执行 Phase 1。先检查仓库、输入、依赖和现有目录，随后严格按照上述规范实施。

再次强调：

```text
STATISTICAL DESIGN IS FROZEN.

DO NOT CHANGE ANY PARAMETER.

DO NOT REDESIGN THE ANALYSIS.

DO NOT INTERSECT CALLABLE SNPS ACROSS ANCIENT SAMPLES.

DO NOT RECOMPUTE PRODUCTION PCA AXES PER ANCIENT SAMPLE.

IF A REQUIRED INPUT, LABEL, TOOL OR ASSUMPTION FAILS VALIDATION,
STOP AND REPORT IT. DO NOT SUBSTITUTE ANOTHER METHOD.

IMPLEMENT, EXECUTE, VALIDATE AND REPORT ONLY.
```
