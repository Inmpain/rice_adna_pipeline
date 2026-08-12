# 生态型 PCA 执行计划（GPT 复核整理，2026-08-12）

> 本文档整理自用户转发的两轮 GPT 策略讨论，经 Claude Code 对照仓库
> `codex/ecotype-pca-panel` 分支当前 commit 的四个实际脚本
> (`check_ref.py` / `map_besthit_to_irgsp.sh` / `pseudo_haploid_call.py` /
> `merge_ancient_into_panel.py`) 和 `main` 分支主流程
> (`scripts/server_originals/mapping.sh`) 核实后写成。**这是一份可执行
> 清单，不是历史记录**——历史讨论、三级祖源框架的完整推导过程仍在
> `docs/ECOTYPE_PCA_PANEL.md`，本文档只写"接下来具体做什么、按什么顺序、
> 每一步的验收标准是什么"。两份文档配合读：本文档回答"下一步做什么"，
> `ECOTYPE_PCA_PANEL.md` 回答"为什么是这个设计"。

---

## 0. 结论先行

GPT 这两轮意见的整体方向是对的，尤其是"先做 IRGSP 覆盖普查，再决定
PCA 用哪个位点集合"这条主线，采纳。但 GPT 对当前代码的两处具体批评
需要先核实再执行——核实结果见第 2 节：**一处批评部分不成立（当前代码
其实没有"根据 check_ref.py 结果全局触发 --swap-ref-alt"这回事，问题
出在文档/注释里给出的判断逻辑，不是代码本身的自动化行为），另一处
批评完全成立且已确认是真实 bug（`merge_ancient_into_panel.py` 的
IndexError 崩溃路径、`map_besthit_to_irgsp.sh` 确实跳过了主流程的
`fixmate`+`-F 0x904` 步骤）**。第 2 节三处代码修复已经在本次提交里
完成，不是待办。

**不新开 git 分支**——包括本文档新增的 IRGSP 覆盖普查/样本专属子集
工作，以及第 10 节的旱稻/水稻生态型 panel (Lyu 2014) 工作，理由和
适用条件见第 10 节末尾。

---

## 1. 下一步执行顺序（10 步，合并两轮 GPT 建议 + 本次核实结果）

1. **【已完成，见第 2 节】** 修复三处代码问题：REF/ALT 方向判定文档
   纠正、比对去重与主流程对齐、合并脚本硬检查+原子写入。
2. 去 `codex/oryza-competitive-mapping` 分支核实 `target_aa_complex`
   /`oryza_genus`/`other_oryza` 分级输出是否已经存在
   (`ECOTYPE_PCA_PANEL.md` 3.2 节待办 0b)。**这一步必须先做**，否则
   Phase 0 覆盖普查的输入范围本身就是错的。
3. 运行 **Phase 0：IRGSP 覆盖普查**（第 3 节）——16 个样本全部跑
   `map_besthit_to_irgsp.sh`（已修复版）+ 新增的
   `summarize_irgsp_coverage.sh`，**不运行 PCA**。
4. 输出 `sample_panel_overlap.tsv` 总表：每个样本 × 每个 panel
   (Civáň/404K/4.8M/29M/720) 的真实可调用位点数（见第 3 节表结构）。
5. 根据第 4 步实测结果，选 2 个代表样本（覆盖最好 + 覆盖最差）。
6. 继续推进三个 panel 的现代群体标签工作（`ECOTYPE_PCA_PANEL.md`
   3.2 节待办 1/2/2c，顺序 2→1→2c），与第 3-5 步并行，互不阻塞。
7. 开发 `build_sample_panel_subset.py`（第 5 节）+
   `simulate_leaveoneout_projection.py`（第 6 节），先在 Civáň 面板上
   用第 5 步选出的 2 个样本冒烟。
8. 用已知现代样本做 leave-one-out 模拟正对照（第 6 节），**同时这一步
   也是验证 REF/ALT 方向是否正确的权威检验**，见 2.1 节。
9. 标签匹配 + leave-one-out 正对照都通过后，Civáň、3K、720 三套标签
   匹配同时推进，冒烟扩大到全部 16 个样本 × 3 个 panel。
10. 根据 16 个样本在 404K/4.8M/2.36M(Civáň)/6.7M 中的实测可调用位点数，
    决定每个样本最终用哪个密度的 panel（第 4 节的优先级表），不预设
    "29M 最密就该用 29M"。

---

## 2. P0：三处代码修复（本次提交已完成）

### 2.1 REF/ALT 方向判定——不能用 check_ref.py 的 FASTA 匹配率来决定 --swap-ref-alt

**核实结果**：当前 `pseudo_haploid_call.py` 的 `--swap-ref-alt` 本来
就是一个需要手动传入的 CLI flag，代码里**没有**"读取 check_ref.py
输出、自动设置 swap"这种全局联动逻辑——GPT 描述的"当前逻辑：
`if args.swap_ref_alt: ref, alt = alt, ref` 不应由 check_ref.py 结果
触发"这句话里，代码本身从来没有被 check_ref.py 结果触发过。**真正的
问题在文档/脚本 docstring 里给使用者的建议**：原来的说法是"先用
check_ref.py 的 snp 模式核对该 panel 的 REF/ALT 方向，如果方向反了就
加 --swap-ref-alt"——这一步推理本身是 GPT 指出的那个混淆：

> check_ref.py 检查的是"`.snp` 哪一列的碱基与 IRGSP FASTA 相同"；
> 真正需要确认的是"`.eigenstratgeno` 的 0/2 分别对应 `.snp` 的哪个
> 等位基因"。

这两者不是一回事：一个 panel 的原始数据完全可能把 REF/ALT 相对基因组
标反了（FASTA 匹配率低），但只要这个"反标"在 convertf 转换全程是
**内部自洽**的（`.snp` 列与已有 `.eigenstratgeno` 矩阵的编码彼此一致），
`pseudo_haploid_call.py` 只要用同样的列约定去解读，新算出来的古代样本
基因型列就依然能正确拼接进现代矩阵——FASTA 匹配率低不代表编码不一致，
匹配率高也不能反过来证明编码一致。

**已落地的修复**：`pseudo_haploid_call.py` docstring 第 4 点已重写
（本次提交），不再建议"看 check_ref.py 的 FASTA 匹配率决定要不要加
--swap-ref-alt"。权威检验改为第 6 节的 **leave-one-out 模拟**——用
panel 里一个已知标签的现代样本，遮蔽成古代样本的覆盖模式，跑一遍跟
古代样本完全相同的 pseudo-haplotype 调用+投影逻辑，如果 0/2 编码方向
真的和该 panel 已有的现代基因型矩阵不一致，这个模拟出来的"伪古代
样本"会系统性地投影偏离它真实的现代群体（往往偏向某个基因型"相反"的
群体，而不只是噪声变大）——这比任何 FASTA 抽查都更直接。**这个模拟
脚本本身还没写，是第 6 节的待办，不是本次修复范围**；本次修复只是
先把错误的判断逻辑从文档里去掉，避免继续按错误依据设置 --swap-ref-alt。

`check_ref.py` 本身不需要改代码——它回答的问题（.snp 列是否匹配
FASTA）仍然是有效信息，只是不能单独拿来决定 swap，仍然值得跑、结果
仍然值得记录（例如 720 panel 91.5% 匹配率这类信息，可以提示"这个
panel 的原始数据来源可能标注习惯不同"，但不能直接当作调用参数依据）。

### 2.2 比对去重步骤——已与主流程对齐

对照 `scripts/server_originals/mapping.sh`（main 分支）核实：

| 项目 | 主流程 `mapping.sh` | 修复前 `map_besthit_to_irgsp.sh` | 修复后 |
|---|---|---|---|
| `bwa aln` 参数 | `-l 1024 -n 0.01 -o 2` | 完全相同 | 不变（确认一致，不再是"未核对"状态） |
| samse 后过滤 | `samtools view -bh -F 0x904` | **没有**，未过滤直接排序 | **已加**，位置与主流程一致 |
| 去重链路 | `collate -O \| fixmate -m - - \| sort \| markdup` | 直接 `samtools markdup` 在坐标排序 BAM 上，fixmate 是"报错再补"的注释 | **已改成与主流程相同的链路**，fixmate 是流程的一部分，不是补丁 |
| markdup 是否 `-r` | `-r`（物理删除重复） | 不删除，只打标记 | **保留不删除**（见下方"刻意保留的差异"） |
| "mapped reads" 计数 | 隐含已被 `-F 0x904` 过滤，`-c` 直接可信 | `samtools view -c`（未过滤，可能把 unmapped 算进去） | **改为 `-c -F 4`**，并新增 MAPQ≥30/MAPQ≥20 两档计数 |

**刻意保留的一处差异**：不加 `-r`。理由：`pseudo_haploid_call.py` 已经
在 pileup 阶段用 `aln.is_duplicate` 过滤重复，物理删除是多余操作；
Phase 0 覆盖普查（第 3 节）还需要重复率本身作为 QC 指标，删掉重复会
丢失这个信息。如果未来某个下游用途确实需要物理去重的 BAM，在调用处
加 `-r`，不改这里的默认值。

MAPQ 30 用于正式 SNP calling（与主流程一致），MAPQ 20 保留作敏感性
分析对照——这与用户原话"正式 SNP 调用建议与原主流程一致使用 MAPQ 30；
MAPQ20可以保留为敏感性分析"完全对应。

### 2.3 合并脚本硬检查——已修复真实 bug

核实确认 GPT 指出的问题都是真实存在的：

- **IndexError 崩溃路径确认存在**：原代码在 `for i, line in
  enumerate(fin): extra = "".join(s[i] for s in call_strings)`
  这一行，如果 panel `.eigenstratgeno` 行数多于 call 文件长度，会在
  `i` 超出 `call_strings` 长度时直接抛出未捕获的 `IndexError`，程序
  在此崩溃退出，**永远走不到后面的行数校验和 warning**——GPT 的描述
  准确。
- **行数不足时静默写出错误文件确认存在**：如果 panel 行数少于 call
  长度，原代码循环正常跑完（只处理到 panel 的行数），随后打印一条
  warning，**但仍然把 `.geno`/`.ind` 写到了最终输出路径**——如果调用者
  没有认真看 stderr，这份不完整的文件会被后续 smartpca 直接吃进去。
- **没有字符校验、没有重复 ID 检查、非原子写入**：确认原代码都没做。

**已落地的修复**（本次提交，见脚本 docstring "SAFETY" 段落）：
1. 行数不匹配现在是**立即硬退出**——多的情况在循环内一旦发现
   `i >= n_snps` 就立刻 `sys.exit` 并清理临时文件，不会走到
   IndexError；少的情况在循环结束后立刻硬退出，同样清理临时文件，
   **最终输出路径不会出现任何不完整文件**。
2. 每个 call 文件读入后立刻校验字符集合是 `{0,1,2,9}` 的子集，否则
   硬退出并报告非法字符和涉及的样本 ID。
3. 新增重复 ID 检查：`--calls` 参数内部重复、以及新样本 ID 与 panel
   `.ind` 里已有现代样本 ID 冲突，两种情况都硬退出。
4. 输出改为"写临时文件（与目标文件同目录，保证 `os.replace` 是同
   文件系统内的原子操作）+ 全部校验通过后原子改名"，中途任何失败路径
   都会清理临时文件，绝不会在目标路径留下部分写入的文件。

---

## 3. Phase 0：IRGSP 覆盖普查（先做，不建 PCA）

**目的**：在决定用哪个 panel、要不要建样本专属子集之前，先搞清楚每个
样本比对到 IRGSP 之后真实覆盖了哪些位置——这一步产出的是事实，不是
建模选择，所以放在所有 PCA 相关工作之前。

**链路**：

```
besthit reads (确认已是 target_aa_complex，见第1节步骤2)
    ↓
map_besthit_to_irgsp.sh（已修复版，第2.2节）
    ↓ 产出：<sample>.besthit_oryza.irgsp.bam（含重复标记，未删除）
summarize_irgsp_coverage.sh（新脚本，待写）
    ↓
与三个 panel（Civáň / 3K的404K-4.8M-29M三档 / 720）分别求交
    ↓
sample_panel_overlap.tsv（总表）
```

**`summarize_irgsp_coverage.sh`（待写，未实现）每个样本至少输出**：

| 类别 | 字段 |
|---|---|
| 比对情况 | 输入 reads、mapped(`-F 4`)、primary mapped、MAPQ≥30、MAPQ≥20、duplicates_flagged、最终保留 reads |
| 全基因组覆盖 | ≥1×覆盖碱基数、覆盖比例、平均深度、≥2×/≥3×覆盖比例 |
| 染色体分布 | chr01–chr12 各自 reads 数、覆盖碱基、平均深度 |
| 窗口分布 | 每 100kb/1Mb 窗口的 reads 数与覆盖度（用于判断是否集中在少数区域） |
| panel 交集 | 每个 panel 各自：覆盖的 panel SNP 数、通过 MAPQ/BQ 数、颠换位点数、最终可调用数 |
| 位点质量 | REF 匹配、ALT 匹配、第三等位基因（对应 `pseudo_haploid_call.py` 新拆分的 `allele_mismatch` 字段）、重复区域比例、低 MAPQ 比例 |

**`sample_panel_overlap.tsv` 表结构（第 4 步产出）**：

```
sample_id  panel        total_panel_snps  covered_snps  tv_only_covered  callable_after_qc
LV...      civan_2.36M   2365188           ...           ...              ...
LV...      3k_404k       404xxx            ...           ...              ...
LV...      3k_4.8m       4.8Mxxx           ...           ...              ...
LV...      3k_29m        29635224          ...           ...              ...
LV...      720_6.7m      6.7Mxxx           ...           ...              ...
```

不做"16 个样本共同覆盖交集"——古 DNA 覆盖本来就稀疏，16 个样本再求
公共交集大概率趋近于零，这条路径直接不采用（GPT 原话已指出，采纳）。

---

## 4. Panel 选择策略——不默认拿 29M 直接建轴

29M 面板适合"提高古样本命中现代变异位点的概率"，但不等于适合直接
拿全部 2900 万位点建 PCA 轴：SNP 高度连锁、局部高密度区域被过度加权、
稀有位点很多、纯文本 EIGENSTRAT 全量大小约 90GB、`smartpca` 计算量大。
EIGENSOFT 自带 `ldregress`/`killr2` 控制 LD、并提醒检查 PC 与缺失率
的相关性。

**优先级表（按第 3 节实测覆盖数决定用哪个，不预设越密越好）**：

| 面板 | 目的 |
|---|---|
| Civáň 2.36M | 先验证栽培/野生大群体归属（PCA-C，桥接层） |
| 3K 404K CoreSNP | 已做过两轮 LD pruning 的规范化 PCA 集合 |
| 3K 4.8M filtered | 404K 覆盖不足时的中间方案 |
| 3K 29M | 高密度探索/群体频率似然用途，不默认直接建轴 |
| 720 6.7M | 野生群体精细分析（PCA-B） |

先看第 3 节实测的 overlap 数，再决定是否需要从 404K 扩到 4.8M 或
29M——不是每个样本都必须用同一档。

---

## 5. 样本专属 panel 子集设计

**逻辑**（每个"样本 × panel"各生成一份）：

```
panel 全部 SNP
    ↓
只保留该古样本真正覆盖的坐标（第3节覆盖普查的产出）
    ↓
MAPQ≥30、BQ≥20
    ↓
排除 duplicate/secondary/supplementary
    ↓
主分析保留 transversion（TV 轨，第7节）
    ↓
碱基必须等于 panel 两个等位基因之一（否则计入 allele_mismatch，第2.3/7节）
    ↓
从现代 panel 的 .geno 中同步提取相同行（按行号，不按 SNP ID 事后拼接）
    ↓
`merge_ancient_into_panel.py`（已加固版）拼入该古样本的 pseudo-haploid 列
```

**必须同时流式读取 `.snp` 和 `.eigenstratgeno`，按同一行号提取**——这一点
`pseudo_haploid_call.py` 现在的实现已经是这样做的（输出严格保持与
panel `.snp` 相同的行数和顺序），不是待改项，写在这里是为了在
`build_sample_panel_subset.py`（待写）里延续同一个不变量，不要在
"子集化"这一步引入按 SNP ID 事后 join 的写法。

**代价（必须在报告里体现，不能只报 PC 坐标）**：每个样本使用的 SNP
集合不同，PC 轴逐样本不可比——16 个古样本不能画在同一张 PCA 图上直接
比较坐标。因此样本专属 PCA 的输出格式应该是：

```
Sample01:
  最近现代群体: TEJ 76%
  次近: TRJ 18%, ARO 4%, 其他 2%
  第一名与第二名距离差: <数值>
  bootstrap 重复中的分类比例: <数值>
  已知现代样本模拟后的误分类率: <数值>（见第6节）
```

而不是简单地说"Sample01 落在 TEJ 附近"。

**辅助分析（不是主分析，不要求先做）**：用现代样本的 404K CoreSNP 等
LD 剪枝集合建一套固定坐标轴，古样本按各自覆盖情况以大量缺失投影
（`lsqproject` 本来就是为高缺失投影设计的）。如果固定轴投影和样本
专属子集分类结论一致，说服力最强。

---

## 6. 不确定性量化——组合方案，不能只靠单一抽样

如果多数古 DNA 位点深度只有 1，单次或重复 100-200 次随机抽 read 会
在深度为 1 的位点上每次都抽到同一条 read，得到完全相同的基因型，点云
会显得虚假地稳定。**必须组合以下几种，不能只用其中一种**：

1. **多次随机抽 read 重复调用**（现有 `pseudo_haploid_call.py`
   `--seed` 参数已预留接口，循环调用逻辑还没写——这是唯一"重复抽样"
   本身能贡献的不确定性来源，且在深度=1的位点上贡献为零，必须配合
   下面几项）。
2. **位点/染色体区块 bootstrap**——对覆盖到的位点集合按染色体区块
   重采样，而不是对每个位点独立重采样。
3. **现代样本按古样本实际位点和缺失模式降采样**——这正是第 6 节
   leave-one-out 模拟脚本 `simulate_leaveoneout_projection.py`
   （待写）要做的事：取已知 IND/AUS/ARO/TEJ/TRJ 现代样本，只保留该
   古样本覆盖的同一批位点，模拟 pseudo-haploid，投影，看是否仍能
   分回真实群体，计算混淆矩阵。**这一步同时服务两个目的**：(a) 判断
   该古样本的位点集合是否有分辨能力（第 4 节"位点够不够"的判断
   依据，不预设"1000个SNP就够"这种固定阈值）；(b) 验证 2.1 节的
   REF/ALT 方向是否正确。
4. **现代 reads 走同一套比对流程**——用现代样本的原始 reads（而不是
   已知基因型直接降采样）重新走一遍 `map_besthit_to_irgsp.sh` +
   `pseudo_haploid_call.py` 全流程，可以同时反映 IRGSP reference bias，
   比直接对已知基因型做统计学降采样更真实，但成本更高，作为第 3 项
   验证通过之后的加强项，不是冒烟阶段必须。

**判断"位点数够不够"的方法**：不预设固定阈值，而是通过现代样本模拟
确定——如果已知 TEJ 样本在某古样本的位点集合上模拟后仍能稳定分回
TEJ，说明这个位点集合有分辨力；如果现代 TEJ/TRJ 自己都分不开，古样本
也不能细分到 TEJ/TRJ。除了位点总数，还要检查：是否覆盖 12 条染色体、
是否集中在少数窗口、是否集中在高拷贝/重复区域、LD独立位点大约多少、
群体高 FST 位点有多少。

---

## 7. TV/ALL 双轨输出 + 报告字段拆分

默认只调用颠换位点（TV）是保守合理的主结果，特别是在末端损伤修剪
尚未定型的现状下（见 `pseudo_haploid_call.py` docstring 第 3 点，
本次已补充双轨运行的操作说明）。**每个样本对每个 panel 应产生两版**：

- **TV**：`pseudo_haploid_call.py` 默认参数（transversion-only），主结果
- **ALL**：加 `--no-transversions-only`，敏感性结果——**注意**：ALL
  不等于"损伤已修剪"，只是"不排除 transition 位点"，因为本项目目前
  没有末端损伤 trimming/rescaling 步骤（见第 10 节 DROT1 的注意事项，
  这是同一个限制在功能位点分析里的具体表现）。

如果 TV 和 ALL 群体结论一致，可信度明显提高；不一致时以 TV 为准，
ALL 的分歧本身是一个需要报告的信号，不是取舍谁对谁错的问题。

**报告字段拆分**（`pseudo_haploid_call.py` 本次已实现，见 `--report`
输出）：不再把"无覆盖"和"碱基不匹配 panel 任一等位基因"混在一个
`n_uncovered` 里，现在拆成：

```
total_panel_snps
transition_skipped
no_allele_info
no_coverage        # 真正没有 read 覆盖
allele_mismatch    # 有 read 覆盖，但碱基不是 panel 两个等位基因中任一个
uncovered_total    # = no_coverage + allele_mismatch，向后兼容口径
called
missing
```

`build_sample_panel_subset.py`（待写）应该在样本专属子集报告里同样
拆分这两个字段，不要重新合并。

---

## 8. 标签方案——标签不是冒烟测试的硬阻塞

最终分析必须有真实标签，但不妨碍机制测试。冒烟阶段：

```
现代样本 → 统一标为 Modern
古样本   → 统一标为 Ancient
poplistname.smoke.txt → 只写一行 "Modern"
```

用于验证 `merge → smartpca → ancient projection` 整条链路是否跑通，
不代表最终结果。真实标签（IND/AUS/ARO/TEJ/TRJ/OrA-OrF/Civáň 的
栽培-野生-亚群分类）与冒烟测试**并行推进**，不是冒烟通过之后才开始
——具体待办见 `ECOTYPE_PCA_PANEL.md` 3.2 节待办 1/2/2c 与本文档第 1
节步骤 6。

---

## 9. `par.PROJECT.template` 需要补充的参数

现有模板故意留空 `poplistname`（等标签工作完成）。除此之外，模板还
应包含（如果当前版本没有，视为待补）：

```
numchrom:       12
numthreads:     8
numoutlieriter: 0
lsqproject:     YES
```

`numoutlieriter: 0` 是因为古代投影样本不参与建轴、也不应该被 smartpca
默认的异常值剔除迭代误伤；`numchrom: 12` 对应水稻染色体数，避免
smartpca 默认的人类 22 条染色体假设漏掉或错配位点。

---

## 10. 后续工作线：旱稻/水稻生态型 panel（Lyu et al. 2014 等）—— 是否新开分支

### 10.1 现状判断

目前没有一个像 3K 亚群那样现成、统一的"亚洲稻生态型 EIGENSTRAT
panel"，需要自己整理。第一优先级是 **Lyu et al. 2014**
（166 份材料，84 upland + 82 irrigated，全基因组重测序，`SRA066116`，
2,623 个生态分化 SNP，74 个 EDR，154 个候选基因）——这是唯一真正按
upland/irrigated 选样并做全基因组比较的数据集。Zhao 2018（3K+202份
upland landrace+446份野生稻，根系表型）、DROT1 2022（271份japonica，
59 upland/212 lowland，抗旱表型）、AfricaRice 2018（330份非洲材料）
作为补充/验证，不作为主 panel。

### 10.2 设计要点（采纳 GPT 的三点核心判断）

1. **不能直接把全部 upland 对全部 irrigated 做一次 PCA**——Lyu 论文
   本身证明 indica-japonica 分化（(F_{ST}=0.13)，japonica内部）远大于
   生态型分化（全体 upland/irrigated (F_{ST}=0.06)），全体 PCA 的 PC1
   会先分开 indica/japonica 而不是旱稻/水稻。正确设计是**两层**：
   先用本分支已有的三级祖源框架（Civáň→3K/720）判定古样本更接近
   IND/AUS/ARO/TEJ/TRJ 或野生稻，再在相同遗传背景内部套用生态型轴
   （比如判定为 japonica 后，只用 `upland japonica vs irrigated
   japonica` 建轴投影）。
2. **旧论文坐标需要转换/核对**——Lyu 2014 用 `IRGSP/RAP build 5`，
   Zhao 2018 用 RGAP7，都不能直接和当前 `irgsp.fa` (IRGSP-1.0) 混用，
   需要先转换坐标或用 `check_ref.py` 核对。
3. **DROT1 的关键变异是 C/T transition**——如果按第 7 节默认的
   TV-only 主结果，这个位点会被直接排除，只能出现在 ALL 敏感性轨，
   不能指望它出现在主分析里。

### 10.3 建议建立三个互补 panel（记录设计，尚未开始）

- `ECO_GLOBAL_LYU`：Lyu 166 份全体，方法复现/全局参考
- `ECO_JAPONICA`：Lyu 的 japonica 子集 + DROT1 的 271 份 japonica，
  主要生态型投影用
- `ECO_FUNCTIONAL`：2,623 个 EDS + Zhao 2018 的 5,779 个根系候选基因
  SNP + DROT1 等已知位点，用于基因型计分，不用于建全基因组 PCA

**与低覆盖古 DNA 匹配的运行方式**：不是给每个古样本永久建一个不同的
现代 panel，而是复用第 5 节"样本专属子集"的同一套逻辑——
`panel SNP ∩ 该样本可调用位点 ∩ transversion`，在交集上重新做 LD
剪枝，现代样本建轴，古样本投影，并用现代样本做同位点数/同缺失模式
的遮盖模拟验证（第 6 节 leave-one-out 逻辑直接复用，不需要另写一套）。

### 10.4 是否新开分支——不需要

理由，与 `ECOTYPE_PCA_PANEL.md` 第 5.4 节已经确立的先例一致（当时
针对 qSH1/qSD1-2/SPS1 功能位点分析线做过同样的判断）：

- **目标相同**：Lyu/DROT1 生态型 panel 最终还是回答"这个古样本更接近
  旱稻还是水稻"，是本分支既有目标的延伸，不是新目标。
- **基础设施完全复用**：第 5 节样本专属子集脚本、第 6 节 leave-one-out
  模拟脚本、`map_besthit_to_irgsp.sh` 产出的 BAM，生态型 panel 全部
  直接复用，不需要另起一套输入链路。
- **拆分支会增加协作成本而无对应收益**——按仓库既有的
  `github-repo-protocol`，每个活跃分支需要维护自己的 handoff 文档、
  与其他分支同步；这项工作目前还处于"整理外部数据+坐标转换"的原型
  验证阶段，尚不确定 Lyu panel 最终能提供多少额外分辨力（取决于第 6
  节模拟结果），此时拆分支等于在验证价值之前先承担维护成本。

**建议做法**：在本分支新建子目录 `scripts/ecotype_panel_lyu/` 做数据
整理和坐标转换原型，复用第 5/6 节的子集化与模拟脚本验证是否真的有
区分力。**只有在验证后确认这是一条需要长期独立投入、且提交历史会
显著与本分支核心工作（覆盖普查/PCA投影）混杂到难以追踪的独立数据集
整理工作时，才考虑拆分支**——参照 `codex/oryza-screen-merge` 当初从
`codex/oryza-competitive-mapping` 分离出去的先例（那是因为它变成了
一个真正独立、上游的小工具，而不是因为主题不同就拆分）。
