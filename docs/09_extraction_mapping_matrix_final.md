# 提取方法 × 定量比对工具 9格矩阵测试 —— 完整记录 (最终版)

记录范围：从最初的方案设计，到全部9个组合(3种提取×3种定量比对)测试完成、
可视化、"单样本最优组合"细化分析为止。

最后更新: 2026-07-29

---

## 一、测试目的

系统性回答两个独立问题：
1. 哪种提取方法(Bowtie2旧参数 / Bowtie2新参数`-N1` / BWA aln)从capture panel
   原始reads里捞出的疑似水稻reads最多？
2. 不同提取方法产出的reads，再用不同定量比对工具(BWA aln / Bowtie2旧参数 /
   Bowtie2新参数)比对回irgsp.fa，最终有效数据量(q30 reads/覆盖度/基因命中数)
   分别是多少？

## 二、核心设计前提

- **shotgun部分固定不变**，不参与任何变量测试(原始shotgun数据不在本地，
  无法重新提取)，所有9个组合共用同一份shotgun reads
- **capture panel1/2**才是提取方法这个变量真正作用的对象
- 阶段①(提取)参考基因组: `asian_rice_panel.fa`(多物种鉴定panel)
- 阶段②(定量比对)参考基因组: `irgsp.fa`(IRGSP1.0)
- **⚠️必须显式排除MCP reshotgun proxy样本**：历史提取产出目录里混杂了14个
  MCP proxy文件(命名格式`LV{数字}-LV{数字}-proxy...`，用短横线，跟真实
  panel文件`LV{数字}_RicePanel{1|2}...`用下划线的格式不同)，不排除会污染
  统计结果

## 三、9格矩阵与结果概览

| 提取方法 \ 定量比对工具 | BWA aln | Bowtie2旧参数 | Bowtie2新参数(-N1) |
|---|---|---|---|
| **Bowtie2旧参数** | 历史数据 | 新增 | 新增 |
| **Bowtie2新参数(-N1)** | 新增 | 新增 | 新增 |
| **BWA aln** | 历史数据(现有主线) | 新增 | 新增 |

**总q30 reads热力图**(16样本加总):

```
                bt2old_map   bwa_map    bt2new_map
bt2_old提取       约?        7708         ~
bt2_new提取       ~          ~            ~
bwa提取           ~          ~            10651最高
```
(实际数字见 `summary/final_mapping_summary.tsv`，具体图见
`summary/heatmap_q30_total.pdf`)

## 四、脚本清单(以服务器实际文件为准，标注最终权威版本)

```
scripts/
├── 01_setup_dirs_and_symlinks.sh          # 建目录+软链接历史①②组合
├── 02_extract_bt2_new.sh                  # [已弃用] 顺序执行版, 保留仅作参考
├── 02a_extract_bt2_new_single.sh          # 单文件Bowtie2新参数提取逻辑
├── 02b_submit_extract_bt2_new.sh          # sbatch批量提交提取作业(实际使用)
├── 03_build_bt2_index_irgsp.sh            # 给irgsp.fa补建Bowtie2索引
├── 04_prepare_reads_combined.sh           # [旧, 有bug] 不识别fastq子目录+未过滤MCP
├── 04_prepare_reads_combined2.sh          # [旧, 部分修复] 识别fastq子目录, 未过滤MCP
├── 04_prepare_reads_combined3.sh          # ✅【最终权威版】fastq子目录兼容+MCP过滤都有
├── 05_run_final_mapping_single.sh         # 单(组合,样本)比对+merge+去重+q30逻辑
├── 06_submit_final_mapping_matrix.sh      # sbatch批量提交7个新组合x16样本=112作业
├── 07_summarize_extraction.sh             # [需要与04同批修复] 提取阶段双口径统计
├── 08_summarize_final_mapping.sh          # ✅定量比对阶段9组合统计(已修复软链接遍历bug)
├── 09_visualize.py                        # 中文版可视化(需要CJK字体配置)
├── 09_visualize_en.py                     # ✅【推荐使用】英文版, 无字体依赖
├── 10_recompute_dup_for_historical_combos.sh  # 补算①②历史组合的真实dup_count
└── 11.best_combo_per_sample.py            # 单样本最优组合分析(避免跨样本误读)
```

### ⚠️建议在git提交前做的清理

```bash
cd /home/scratch/yinmt202607/tests/param_matrix_bt2_vs_bwa/scripts
rm -f 04_prepare_reads_combined.sh 04_prepare_reads_combined2.sh
mv 04_prepare_reads_combined3.sh 04_prepare_reads_combined.sh
rm -f 02_extract_bt2_new.sh   # 已被02a+02b的sbatch版本取代
# 09_visualize.py(中文版)可以保留作为备选，也可以删除只留英文版，看个人偏好
mv 11.best_combo_per_sample.py 11_best_combo_per_sample.py   # 统一命名风格(下划线)
```

清理后重新生成一份`README.md`放在scripts目录下，注明"编号即建议运行顺序，
04已合并为单一权威版本"，避免以后自己或别人对着一堆历史版本发懵。

## 五、完整执行步骤(以清理后的最终脚本为准)

```bash
cd /home/scratch/yinmt202607/tests/param_matrix_bt2_vs_bwa

# 1. 建目录+软链接历史数据
bash scripts/01_setup_dirs_and_symlinks.sh

# 2. 给irgsp.fa建Bowtie2索引(仅需一次)
bash scripts/03_build_bt2_index_irgsp.sh

# 3. sbatch批量提交Bowtie2新参数(-N1)提取
bash scripts/02b_submit_extract_bt2_new.sh
squeue -u $USER   # 等待跑完

# 4. 整理三种提取方法的reads_combined(含fastq子目录兼容+MCP过滤)
bash scripts/04_prepare_reads_combined.sh   # (清理后即最终权威版)

# 5. sbatch批量提交7个新组合x16样本=112个定量比对作业
bash scripts/06_submit_final_mapping_matrix.sh
squeue -u $USER   # 等待跑完

# 6. 汇总统计
bash scripts/07_summarize_extraction.sh
bash scripts/08_summarize_final_mapping.sh

# 7. 补算①②历史组合的真实去重信息(原流程用markdup -r直接删除, 无法直接统计)
bash scripts/10_recompute_dup_for_historical_combos.sh
# (随后需要用python脚本把补算结果拼回final_mapping_summary.tsv, 见decisions_log)

# 8. 可视化(英文版, 无字体依赖)
python3 scripts/09_visualize_en.py

# 9. 单样本最优组合分析(避免跨样本误读)
python3 scripts/11_best_combo_per_sample.py
```

## 六、踩过的坑与解决方案(重要，避免以后重复踩)

### 坑1：`asian_rice_panel.fa`从未建过Bowtie2索引
提取阶段第一次用Bowtie2新参数时报`Could not locate a Bowtie index`。
解决：`bowtie2-build asian_rice_panel.fa asian_rice_panel.fa`(注意实际索引
前缀带`.fa`，之前误以为不带，浪费了一次重复建索引)。

### 坑2：Bowtie2内存不足OOM
`asian_rice_panel.fa`索引总大小约10.4GB(`.bt2l`大索引格式)，最初给
`--mem 4G`导致部分作业被OOM杀掉。解决：调到`--mem 16G`以上，重新提交
(脚本有断点续跑保护，已成功的样本不受影响)。

### 坑3：`sbatch --wrap`里用`source activate`找不到环境
Slurm作业执行环境不会自动加载`~/.bashrc`/conda初始化脚本。解决：改用
`export PATH=<conda环境bin目录>:$PATH`直接把工具路径加入PATH，绕开
conda activate机制。

### 坑4：`04_prepare_reads_combined`脚本不识别fastq子目录
`bt2_old`/`bwa`是软链接历史数据(文件直接在该层)，而`bt2_new`是新生成数据
(文件在`fastq/`子目录下)，两种目录结构不一致导致脚本遍历时漏掉`bt2_new`
的文件。解决：同时遍历`$src`和`$src/fastq`两层。

### 坑5：MCP reshotgun proxy样本污染提取统计
历史提取产出目录(`asian_rice_compare/bowtie2,bwa/fastq`)混杂了14个MCP
proxy样本(因为当初`compare_rice_read_extractors.sh`默认扫描范围包含了
`4.mcp_reshotgun`目录)。解决：用正则`^LV[0-9]+_RicePanel[0-9]`精确匹配
真实angkor panel文件名，显式排除MCP格式(`LV{数字}-LV{数字}-proxy...`)。

### 坑6：脚本08用`*/`通配符遍历软链接目录时行为不一致
bash glob对"以`/`结尾的模式匹配指向目录的软链接"，在不同shell选项设置下
表现不完全稳定，导致①②两个历史软链接组合被漏统计。解决：改用
`find -L "$BASE/02.final_mapping" -mindepth 1 -maxdepth 1 -type d`
显式follow软链接遍历，不依赖glob默认行为。同时①②组合的目录结构本身
(`final`层直接就是数据层，没有再套一层`final`子目录)与新生成的7个组合
不同，脚本里加了自动判断兼容两种结构。

### 坑7：①②历史组合的`markdup`加了`-r`，重复reads被直接删除
导致这两个组合的`dup_count`/`dup_rate_pct`统计出来恒为0(不是真的没有
重复，是物理删除后无法事后统计)。解决：回到markdup之前的中间产物
(`.sorted.bam`)，重新跑一次不加`-r`的markdup单独补算，再拼回主汇总表
(见`10_recompute_dup_for_historical_combos.sh`)。

### 坑8：matplotlib无法正确渲染中文标题
系统里`Noto Sans CJK`字体的`.ttc`文件包含SC/TC/JP/KR/HK五个语言变体名，
但matplotlib的字体管理器只注册其中一个(本例中是"Noto Sans CJK JP")，
写`font.sans-serif = ['Noto Sans CJK SC']`匹配不到会导致中文显示为方框。
最终决定：**改用全英文版可视化脚本(`09_visualize_en.py`)**，彻底避免
字体环境依赖问题，服务器换环境也不需要重新配置字体。

## 七、最终结论

### 7.1 提取阶段(阶段①)：BWA全面领先

3种方法在capture panel1+2上提取出的reads总量(`panel_total`列)排名：
**BWA > Bowtie2新参数(-N1) > Bowtie2旧参数**。Bowtie2加`-N1`比旧参数
确实有提升，但离BWA还有明显差距。

### 7.2 定量比对阶段总量：`bwa_extract__bwa_map`总q30最高

9个组合里，`bwa_extract__bwa_map`(现有主线方案)总q30 reads最高，比
排名最低的组合(`bt2_old_extract__bwa_map`，即最早的历史流程)高出约1.67倍。

### 7.3 【关键细化】单样本最优组合分析：真正起决定作用的是提取方法，不是定量比对工具

用`11_best_combo_per_sample.py`逐样本比较(而非只看总量，避免跨样本
对比造成的误读)，发现：

- q30维度：`bwa_extract__bwa_map`9/16样本第一，`bwa_extract__bt2new_map`
  7/16样本第一
- gene_hit维度：`bwa_extract__bwa_map`8/16，`bwa_extract__bt2new_map`4/16，
  `bt2_new_extract__bt2new_map`3/16，`bwa_extract__bt2old_map`1/16

**所有拿到"单样本最优"的组合，提取阶段清一色都是`bwa`**——没有任何一次
是`bt2_old`或纯`bt2_new`提取拿到冠军(除了gene_hit维度上有3次是
`bt2_new_extract__bt2new_map`)。**这说明真正决定性的变量是提取方法(阶段①)，
不是定量比对工具(阶段②)**：提取阶段用BWA之后，后续接BWA还是Bowtie2新参数
做最终比对，两者互有胜负(9:7，比例接近，不是压倒性优势)；但提取阶段一旦
用Bowtie2(不管新旧参数)，基本没有机会拿到单样本最优。

### 7.4 附加发现：q30最优与gene_hit最优只有6/16样本一致

说明"reads数量最多"不完全等价于"基因命中数最多"，这与之前发现的
"低复杂度序列会使reads数量虚高但命中质量存疑"的结论互相呼应，值得在
后续分析中留意这个区分。

## 八、Git仓库组织建议

```
~/rice_adna_pipeline/
└── tests/
    └── param_matrix_bt2_vs_bwa/
        ├── README.md                          # 本文档精简版/链接
        ├── scripts/                           # 清理后的最终脚本(见第四节)
        └── summary/                           # 只提交这几张小体积结果表+图,
            ├── extraction_summary.tsv           # 不提交BAM等大文件
            ├── final_mapping_summary.tsv
            ├── final_mapping_summary_fixed.tsv   # (含补算dup信息的最终版)
            ├── best_combo_per_sample.tsv
            ├── heatmap_q30_total.pdf
            └── barplot_16samples_9combos.pdf
```

```bash
cd ~/rice_adna_pipeline
mkdir -p tests/param_matrix_bt2_vs_bwa/{scripts,summary}

# 拷贝清理后的脚本和结果表(不含BAM等大文件)
cp /home/scratch/yinmt202607/tests/param_matrix_bt2_vs_bwa/scripts/*.sh \
   /home/scratch/yinmt202607/tests/param_matrix_bt2_vs_bwa/scripts/*.py \
   tests/param_matrix_bt2_vs_bwa/scripts/
cp /home/scratch/yinmt202607/tests/param_matrix_bt2_vs_bwa/summary/*.tsv \
   /home/scratch/yinmt202607/tests/param_matrix_bt2_vs_bwa/summary/*.pdf \
   tests/param_matrix_bt2_vs_bwa/summary/

# 把本文档也放进去
cp 09_extraction_mapping_matrix_final.md docs/

git add tests/param_matrix_bt2_vs_bwa/ docs/09_extraction_mapping_matrix_final.md
git commit -m "Add 9-combo extraction x mapping matrix test: scripts, results, conclusions (BWA extraction is the decisive variable)"
git push
```

## 九、待办/后续建议

- 把本次结论(提取方法是决定性变量)补充进主线`docs/decisions_log.md`，
  作为"为什么坚持用BWA"这个决策的最终、最完整证据
- `final_mapping_summary_fixed.tsv`(含补算dup信息)应作为正式对外的版本，
  原始`final_mapping_summary.tsv`(①②两组dup恒为0)建议标注为中间产物
- 如果后续还要测试更多参数组合，建议复用本次的目录结构和脚本框架，
  避免重新踩一遍第六节列出的8个坑
