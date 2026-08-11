
# Oryza competitive mapping / best-hit 接手说明

更新时间：2026-08-11（**第八次更新：16个样本mapping全部跑完了。开始实测
v2脚本时卡在一个环境问题——服务器上跑`bash submit_oryza_besthit.sh check`
打出来的还是v1那种老格式的一行输出(`ORYZA_TAXIDS=4529 4530 4536
DAMAGE_WINDOW=5 TOP_N=10`)，不是v2新增的`Oryza scope: whole genus,
ORYZA_GENUS_TAXID=4527`那种提示——说明服务器上`scripts/`目录里的
`submit_oryza_besthit.sh`大概率还是8号推送v2之前的旧版本，`curl`重新下载
这一步没有真正生效（或者根本没跑）。**已经让用户在服务器上跑
`grep -c "genus" submit_oryza_besthit.sh oryza_besthit_damage_filter.py`
确认是不是0，同时给了带`-o`强制覆盖文件名的重新下载命令，**还没收到反馈，
下一个接手的session/窗口第一件事就是等/要这个grep结果，不要假设已经解决**。
第七次更新：besthit主脚本v2重写完成——Oryza范围
从硬编码3个种(rufipogon/sativa/nivara)改成动态解析整个Oryza属，取代
`--oryza-taxids`默认值；v1脚本已归档为`oryza_besthit_damage_filter_v1.py`
（`git mv`保留历史）。本地合成数据验证通过，但没有服务器真实BAM/真实数据
的验证，详见新增的5.1b节，含"综合师兄脚本能不能采纳"的逐条取舍表**。第六
次更新（当晚收工前）：服务器真实dry-run已跑过，结果符合预期方向，但脚本
又发现并修复了一个新问题（irgsp.acc2taxid有表头行，会污染合并文件），
--apply还没跑，是明天开工第一件事，详见0.6节。第五次更新：用户截图确认
`asian_rice_panel.acc2taxid`是**整份文件**
4529/4530系统性反标——不是"np7等12个基因组各自标错"，而是"凡标4529的都
应为4530、凡标4530的都应为4529"，全文件二元互换。修正脚本已从"按基因组
名字分别赋值"改成"整列做4529↔4530互换"，逻辑更简单也更可靠（不依赖基因组
身份判断）。第四次更新：找到了`all_wgs_asian_irgsp.acc2taxid`
真实的三源合并配方，把0.5节的修正方案从"补丁式替换"改成"完整重建"，脚本
换成了`rebuild_all_wgs_asian_irgsp_acc2taxid.sh`，旧的`fix_asian_rice_panel_
taxid.sh`已从仓库删除。第三次更新：发现并诊断了`asian_rice_panel.acc2taxid`
里`np7`(IRGSP/Nipponbare)的taxid错标问题（当时以为只是这一个基因组的孤立
问题，后来发现是整文件系统性反标的冰山一角）。第二次更新记录：acc2taxid
Oryza覆盖度诊断结果出来后，证伪了最早"数据库野生稻缺失"的假设，改用"损伤
窗口固定长度 vs 读长"这个新证据更充分的假说，相应重排了第7/8/9节的优先级）

> 读这份文档前提：你没有服务器直接执行权限，所有服务器端命令的输出都需要
> 人类协作者手动跑了贴回来。不要假设你能直接看到服务器文件系统。

## 0. 现状一句话总结

Best-hit 脚本（`oryza_besthit_damage_filter.py` + `submit_oryza_besthit.sh`）
已经在服务器上实际跑通（3个样本里的1个做过1000-read smoke test，其余2个mapping
已完成但besthit还没跑全量）。脚本本身没有已知 bug 了（过程中发现并修复了3个，
见5.6）。**当前卡点不是代码**。

**方法论上走过一轮弯路，已经纠正，写在这里避免下一个人重复踩**：最早怀疑是
"acc2taxid 数据库里 Oryza 参考基因组覆盖不够、大概率只有 sativa"——**这个
假设已经用实测数据证伪**（见7.1）：Oryza 属18个已知种在数据库里普遍有几十到
几百条 contig，*rufipogon*(717条) 反而比 *sativa*(96+15=111条) 覆盖更好。
真正吻合观测数据的假说是**末端损伤校正窗口（`--damage-window`默认5bp）相对
读长是固定的，读长越长、两端5bp之外"没资格被校正"的中段就越大，累积一个
额外错配（无论是真实SNP、非末端损伤、还是测序错误）的概率也越高**——
smoke test里"差1个编辑距离惜败"的541条read，均值/中位数长度（~68bp）明显
比"打平守住KEEP"的103条read（~62bp）更长，方向和这个假说完全吻合（见6.2）。
**这个假说还差最后一步验证**：这批样本有没有做过标准的损伤衰减曲线分析
（末端N个位置的C→T/G→A频率随距离衰减的图，mapDamage风格），衰减到背景水平
大概在多少bp——如果实际损伤延伸远超5bp，说明`--damage-window`设得太窄，
这是下一个人接手后第一件要确认的事（见7.5）。**2026-08-07进展**：4个样本
的全量besthit留存率高度一致（约4.4-5.0%），用户据此判断"5bp暂时可用"，
优先级让位给0.5节的taxid修正（见7.5末尾）。

## 0.5 ⚠️`asian_rice_panel.acc2taxid` 整份文件 4529/4530 系统性反标（2026-08-08，
### 已找到真实三源合并配方，修正方案改为完整重建 + 全文件二元互换）

**问题发现过程分两步，第二步推翻了第一步的范围判断，接手时以第二步为准**：

**第一步（最初发现，范围判断有误）**：用户和GPT的独立讨论（完整过程见
`docs/asian_rice_panel_reference_design_conversation.md`）里，用染色体长度+
逐条染色体MD5比对，确认了`asian_rice_panel.fa`里的`np7.Chr1-12`与`irgsp.fa`
的`chr01-12`**完全一致(IDENTICAL)**——即`np7` = Nipponbare/IRGSP-1.0 =
*O. sativa*。但`asian_rice_panel.acc2taxid`里`np7.*`全部被标成了4529
(*O. rufipogon*)，应为4530(*O. sativa*)。**当时误判为"只是np7这一个基因组
的孤立错标"，据此设计了按基因组名字逐个赋值的修正脚本**。

**第二步（2026-08-08，用户截图纠正范围）**：用户直接贴出`asian_rice_panel.
acc2taxid`文件截图，指出**这是整份文件的系统性反标，不是个别基因组的问题**：
凡是应该标4530(sativa)的行，文件里全标成了4529；凡是应该标4529(rufipogon)
的行，文件里全标成了4530。截图证据（均为文件里的原始误标，箭头是应该改成
什么）：
- `G25_ruf_W2064.*`、`G25_ruf_W3037.*`（野生rufipogon组，本应4529）——
  文件里标成了 **4530**
- `IRGSP-1.0_genome.chr01-12`（sativa标准参考，本应4530）——文件里标成了
  **4529**
- `X24_kas.*`（sativa代表品种，本应4530）——文件里标成了 **4529**

也就是说，真正的规律是**整个4529/4530二元标签对调了**，不需要（也不应该）
按基因组名字分别判断"这个该是哪个taxid"——不管这一行原来是什么名字，只要
原始taxid是4529就该改成4530，原始是4530就该改成4529。原来的"12个基因组
按名字分别赋值"这套逻辑虽然对np7这一个例子给出了正确答案，但只是恰好蒙对
了方向，不是正确的诊断——已经废弃。

**对已有分析的影响**（这部分结论不变，但影响幅度比最初以为的大）：
- **不影响besthit的KEEP/REJECT二元判定**——`--oryza-taxids`默认
  `"4529 4530 4536"`，4529和4530都在白名单里，无论标签方向如何，panel里
  的reads仍会被判定为"目标Oryza"，不会被误判成非Oryza。
- **但影响第6.1节的"各物种contig数"统计表，且影响范围比最初估计（只有
  np7的16行）大得多**——既然是整份`asian_rice_panel.acc2taxid`的标签对调，
  意味着panel里**全部**标记为sativa代表品种的基因组（np7/mh63/X24_kas/
  azu/arc/liuxu）此前很可能都被计入了rufipogon桶，全部标记为野生rufipogon
  的7个`G25_ruf_W*`基因组则可能被计入了sativa桶——不是16行的量级，而是
  这份panel文件的**全部行数**。6.1节"sativa 111条 / rufipogon 717条"这张
  表在修正前完全不可信，必须等服务器重新统计后才能用。

**⚠️`all_wgs_asian_irgsp.acc2taxid`的真实三源合并配方已经在服务器上找到**
（原以为合并脚本没有留存，只能"补丁式"局部替换——这个判断是错的，用户在
服务器`taxonomy_CPH`目录下确认了三个源文件，修正方案改成了**完整重建**，
不再是局部补丁）：

1. **WGS真核库taxid**（独立大文件，406M，跟panel目录不在一起）：
   `/home/database/ref20250728/taxonomy_CPH/wgs_eukaryota.acc2taxid`
2. **亚洲水稻panel taxid**：`<panel目录>/asian_rice_panel.acc2taxid`
   （整文件4529/4530对调，本次修正对象）
3. **IRGSP taxid**：`<panel目录>/irgsp.acc2taxid`（独立文件，是否也有同样
   错标问题**还没验证过**，新脚本只统计报告，不自动改，见下方）

**修正方案**：`scripts/oryza_besthit/rebuild_all_wgs_asian_irgsp_acc2taxid.sh`
——本地已用还原截图里那三种反标模式（rufipogon组标成4530、IRGSP-1.0_genome
和X24_kas标成4529）的合成数据测试过dry-run和--apply两种模式，验证互换方向
完全正确、三源拼接行数吻合、备份文件齐全，逻辑是：
1. 对`asian_rice_panel.acc2taxid`第3列做**整列二元互换**：原来是4529的行
   改成4530，原来是4530的行改成4529，**不再依赖基因组名字判断**，写出
   修正后的副本
2. 打印`irgsp.acc2taxid`前5行和第3列(taxid)的分布统计，**不自动修改**——
   需要人工看这份输出，结合这次的教训（"看起来不对的时候，真实范围往往
   比第一眼判断的更大"）确认它是不是也有同样的整体反标问题，再决定要不要
   另外处理
3. `--apply`模式下：备份原始`asian_rice_panel.acc2taxid`和
   `all_wgs_asian_irgsp.acc2taxid`，用修正后的panel文件替换原文件，然后
   **完整重建**`all_wgs_asian_irgsp.acc2taxid` = `cat wgs_eukaryota.acc2taxid
   asian_rice_panel.acc2taxid(修正后) irgsp.acc2taxid`（原样拼接，不改
   WGS和IRGSP各自的内容）

**已从仓库删除**：`scripts/oryza_besthit/fix_asian_rice_panel_taxid.sh`
（最早的"补丁式局部替换"方案，基于"合并脚本已丢失"这个错误前提）和第一版
`rebuild_all_wgs_asian_irgsp_acc2taxid.sh`（基于"只有12个基因组按名字错标"
这个范围误判，commit `3edf2dbc`）都已被当前版本取代（commit `8ba66b4`）。

⚠️基因组**身份**验证（哪个accession实际是哪个品种/组）和taxid**方向**修正
是两件独立的事——身份验证方面仍然**只有np7做过逐染色体MD5的独立验证**，
其余基因组(mh63/X24_kas/azu/arc/liuxu/7个G25_ruf_W*，以及截图里出现的
独立`IRGSP-1.0_genome`条目和`np7`是否指向同一份序列)的具体身份仍待确认，
不受这次taxid方向修正影响——但taxid方向的修正本身现在不再依赖身份判断，
已经不是这个开放问题的阻塞项了。

**待办**：
```bash
# 1. 先 dry-run，确认改动符合预期(不写文件)——重点看 irgsp.acc2taxid 的
#    taxid分布输出，贴回来确认它有没有同样的整体反标问题
curl -fsSL -o rebuild_all_wgs_asian_irgsp_acc2taxid.sh \
  https://raw.githubusercontent.com/Inmpain/rice_adna_pipeline/codex/oryza-competitive-mapping/scripts/oryza_besthit/rebuild_all_wgs_asian_irgsp_acc2taxid.sh
bash rebuild_all_wgs_asian_irgsp_acc2taxid.sh \
  --wgs /home/database/ref20250728/taxonomy_CPH/wgs_eukaryota.acc2taxid \
  --panel-dir /home/scratch/yinmt202607/db/asian_rice_panel_index

# 2. 确认无误后加 --apply 真正落盘(会自动备份原文件)
bash rebuild_all_wgs_asian_irgsp_acc2taxid.sh \
  --wgs /home/database/ref20250728/taxonomy_CPH/wgs_eukaryota.acc2taxid \
  --panel-dir /home/scratch/yinmt202607/db/asian_rice_panel_index \
  --apply

# 3. 重新统计第6.1节的物种contig数表，更新这份文档
```

## 0.6 ⚠️2026-08-08当晚：真实dry-run已跑，结果符合预期，但发现新问题——
### `--apply`还没跑，是明天开工第一件事

**服务器dry-run实测结果**（用户跑的是当时的版本，commit `8ba66b4`，还没有
下面提到的表头过滤修复）：

```
== Step 1 ==
修正前: 705行=4529, 84行=4530  (changed_lines=789，即整份文件都改了)
修正后: 84行=4529, 705行=4530
抽样: G25_ruf_W0169.Chr1  4530 -> 4529   (对，野生rufipogon本应4529，
                                          之前被错标成4530，方向正确)

== Step 2 ==
irgsp.acc2taxid 前5行第一行是表头: accession/accession.version/taxid
真实12行数据全部是 4530 (chr01-chr12)

== Step 3 ==
wgs_eukaryota.acc2taxid:              11,038,401 行
asian_rice_panel.acc2taxid(修正后):          789 行
irgsp.acc2taxid:                              13 行(12数据+1表头)
现有all_wgs_asian_irgsp.acc2taxid:    11,039,202 行
```

**两个关键结论**：
1. **互换方向确认完全正确**——不是靠合成数据推断，是服务器真实文件上
   验证过的：789行里705行原本标4529、84行原本标4530，互换后709行变
   84/705变705，抽样对照的`G25_ruf_W0169`(野生rufipogon)从错误的4530
   变成正确的4529。
2. **`irgsp.acc2taxid`本身没有taxid反标问题**——12行真实数据(chr01-12)
   全部正确标成4530(sativa)，不需要修正taxid方向。**但**它有一行字面
   表头(`accession	accession.version	taxid`)，如果被原样`cat`进
   1100万行的合并文件，第3列会变成字符串`"taxid"`而不是数字，besthit
   脚本按taxid解析这一列时会出问题。

**已修复**（commit `2965348`，当晚新增的Step 2b + 过滤逻辑）：脚本现在
会对三个源文件都检查"第3列是不是纯数字"，`--apply`阶段自动丢弃任何
非数字taxid的行（不只是irgsp，wgs和panel文件如果也有类似问题也会被
一起处理）。**这个修复只在本地用合成数据测试过（还原了"文件里有一行
字面表头"这个场景），没有在服务器上重新跑过dry-run确认**。

**明天开工的第一步，是重新下载这个修复后的脚本、重新跑一次dry-run确认
Step 2b显示`irgsp.acc2taxid 非数字taxid行数 = 1`（这次应该被正确识别
并即将被过滤），确认没有新问题后，再跑`--apply`**：

```bash
cd /home/scratch/yinmt202607/gene/scripts
curl -fsSL -o rebuild_all_wgs_asian_irgsp_acc2taxid.sh \
  https://raw.githubusercontent.com/Inmpain/rice_adna_pipeline/codex/oryza-competitive-mapping/scripts/oryza_besthit/rebuild_all_wgs_asian_irgsp_acc2taxid.sh

# 重新dry-run一次，确认Step 2b正确识别irgsp.acc2taxid的表头行
bash rebuild_all_wgs_asian_irgsp_acc2taxid.sh \
  --wgs /home/database/ref20250728/taxonomy_CPH/wgs_eukaryota.acc2taxid \
  --panel-dir /home/scratch/yinmt202607/db/asian_rice_panel_index

# 确认无误后--apply（会自动备份原文件，预期新的all_wgs_asian_irgsp.acc2taxid
# 行数应该是 11,038,401 + 789 + 12 = 11,039,202（跟旧文件行数一样，因为
# 只是替换了12个原有基因组自己内部的标签方向，加上去掉了irgsp那1行表头，
# 净变化应该是"减1行"——如果看到的不是11,039,201，先别急着继续，回来核对）
bash rebuild_all_wgs_asian_irgsp_acc2taxid.sh \
  --wgs /home/database/ref20250728/taxonomy_CPH/wgs_eukaryota.acc2taxid \
  --panel-dir /home/scratch/yinmt202607/db/asian_rice_panel_index \
  --apply
```

**--apply之后要做的事**（还没做）：
1. 重新统计第6.1节的物种contig数表，替换掉现在标着"⚠️待修正后重新统计"
   的旧数字
2. 重新评估"sativa覆盖度偏低"这个结论在真实数字下是否仍然成立（第0/6.1/
   7.1节现在都还是基于旧的、已知有整体反标问题的数字，方向可能会变）
3. 更新本文档第5.6节bug表格第4行的commit hash（目前写的是"待补"/旧commit，
   最终版本是`2965348`）

## 1. 项目背景

古 DNA / 环境 DNA 水稻项目（angkor，16个考古样本，年代跨度约公元1160-1981年，
元数据见服务器 `/home/scratch/yinmt202607/angkor_robot_library.txt`）。候选水稻
reads 已经由 shotgun、capture panel1、panel2 三份 FASTQ 合并完成。当前正在把
每个样本分别比对到 WGS 真核数据库、亚洲水稻 panel 和 IRGSP（competitive
mapping），然后执行 best-hit，排除更像非水稻物种的 reads，产出"确认是 Oryza"
的干净 FASTQ 供下游使用。

## 2. 关键路径

### GitHub

仓库：https://github.com/Inmpain/rice_adna_pipeline

当前分支：`codex/oryza-competitive-mapping`

本机（Mac）GitHub 工作副本：`/Users/inmpain/Documents/angkor/rice_adna_pipeline_publish`

本机参考代码（**不在这个 git 仓库里，是另外两个独立项目，只用来抄逻辑，不是
本项目依赖**）：
- 原 best-hit 项目（ngsLCA 风格通用分类器，最终**没有**直接复用，见第5节）：
  `/Users/inmpain/github/aeDNA_popgen`
- WGS mapping/动态内存参考：`/Users/inmpain/github/Pipeline_snakemake/new_single_multi`

### 服务器（`/home/scratch/yinmt202607/gene/scripts/`）

这个目录下当前实际存在的文件（用户在服务器上 `ls` 看到的，注意有历史版本
文件，不要跑错）：

```
merge_oryza_fastq.sh
oryza_besthit_damage_filter.py          # <- besthit 主脚本，本次工作的产出
oryza_screen_merge/
submit_oryza_besthit.sh                 # <- besthit 提交脚本，本次工作的产出
submit_oryza_competitive_mapping_old.sh
submit_oryza_competitive_mapping.sh     # <- 当前在用的 mapping 提交脚本
submit_oryza_competitive_mapping_v2.sh
submit_oryza_competitive_mapping_v3.sh
test_submit_oryza_competitive_mapping.sh
test_submit_oryza_competitive_mapping_v2.sh
test_submit_oryza_competitive_mapping_v3.sh
```

**未确认**：`submit_oryza_competitive_mapping.sh`（无版本号后缀）是否就是
git 仓库 `scripts/oryza_besthit/submit_oryza_competitive_mapping.sh` 的当前
版本，还是本地又手改过、和 git 不同步——git仓库里只有一份（无 v2/v3），服务器
上有多个版本号文件，说明 mapping 阶段的脚本在服务器上有过 git 之外的迭代。
接手时先 `diff` 一下确认，不要默认服务器和 git 一致。

服务器上下载 git 脚本的方式（这台机器 login node 能连外网 github.com，已验证
可用）：
```bash
cd /home/scratch/yinmt202607/gene/scripts
curl -fsSL -o <文件名> \
  https://raw.githubusercontent.com/Inmpain/rice_adna_pipeline/codex/oryza-competitive-mapping/scripts/oryza_besthit/<文件名>
```

服务器 python 环境：conda/mamba 环境 `/home/usr/yinmt/.local/mamba/snakemake`
（提示符显示 `(base)` 或 `(.../snakemake)`），`python3` 直接在 PATH 上，已确认
装了 `pysam`。**没有 `python/` modulefile**（`module load python/` 会报错，
脚本里已经用 `|| true` 吞掉了，见第4.1节）。

候选 FASTQ：`/home/scratch/yinmt202607/gene/results/oryza_candidates_combined/<sample>.oryza_candidates.combined.fastq.gz`

mapping 输出根目录：`/home/scratch/yinmt202607/gene/results/oryza_competitive_mapping/`
- 各数据库 BAM：`bam_by_database/<sample>/<sample>.<database>.bam`
- 每样本合并 BAM（besthit 的输入）：`by_sample/<sample>.competitive.name_sorted.bam`（+ `.finished` 标记）
- mapping 日志/提交记录/串行状态：`logs/` `submissions/` `series/`

besthit 输出根目录：`/home/scratch/yinmt202607/gene/results/oryza_competitive_mapping/besthit/`
（脚本会自动建，见第5节的输出说明）；smoke test 单独输出到 `besthit/smoke_test/`。

## 3. 数据库

WGS 真核数据库：`/home/database/ref20250728/cph_euk/wgs_eukaryota.{1..129}.fas.gz`（129个shard）

亚洲水稻 panel：`/home/scratch/yinmt202607/db/asian_rice_panel_index/asian_rice_panel.fa`
（12个基因组：np7/mh63/X24_kas/azu/arc/liuxu 6个O. sativa代表品种 + 7个
G25_ruf_W*的O. rufipogon野生稻，选取理由/覆盖的遗传谱系/质量评估见
`docs/asian_rice_panel_reference_design_conversation.md`）

IRGSP：`/home/scratch/yinmt202607/db/asian_rice_panel_index/irgsp_bt2idx`

合并 accession→taxid：`/home/scratch/yinmt202607/db/asian_rice_panel_index/all_wgs_asian_irgsp.acc2taxid`
（标准 NCBI 三列格式：accession / accession.version / taxid。由三个源文件
`cat`拼接而成——`/home/database/ref20250728/taxonomy_CPH/wgs_eukaryota.
acc2taxid` + `asian_rice_panel.acc2taxid` + `irgsp.acc2taxid`，见0.5节。
**⚠️2026-08-08发现np7标签错误，修正+重建方案见0.5节；irgsp.acc2taxid是否
有同样问题还没验证，重建脚本会打印其taxid分布供人工确认**）

NCBI taxonomy dump：`/home/database/ref20250728/taxonomy_CPH/ncbi/20250530/{nodes.dmp,names.dmp}`

目标 Oryza taxid（species rank，**已用 --oryza-taxids 默认值写死在脚本里，
但 rank 是不是真的 species 没有校验过，见第7节待办1**）：
- 4529：*Oryza rufipogon*（野生稻/普通野生稻）
- 4530：*Oryza sativa*（亚洲栽培稻）
- 4536：*Oryza nivara*（尼瓦拉野生稻）

### ⚠️ 潜在的额外野生稻参考资源（还没用上，见第7/8节）

`file_path.md`（`main` 分支）提到服务器上 `/home/scratch/yinmt202607/db/3k/wild/`
下有 **140+ 个 `{SampleID}.transfer.merge.chr.fasta`（野生稻/近缘种染色体级
组装基因组）**，来自 3K Rice Genome Project 相关资源，**目前完全不在
`all_wgs_asian_irgsp.acc2taxid` / competitive mapping 的131个数据库体系里**。
这些样本具体是什么物种、有没有已知的 taxid，仓库文档里没写，需要上服务器查
（命令见第7节）。如果这批组装真的能用，是目前看来解决"Oryza 参考单一化"
问题最直接的路径。

## 4. Mapping 阶段（competitive mapping，besthit 的上游）

每样本比对131个数据库（129 WGS shard + asian_rice_panel + IRGSP），Bowtie2：
```bash
-k 100 -L 22 -i S,1,1.15 --mp 1,1 --rdg 0,1 --rfg 0,1 --score-min L,0,-0.1 --no-unal
```

截至 2026-08-05，`by_sample/` 下已完成（`.finished` 标记齐全）的样本：
- LV6000619499（唯一做过 besthit smoke test 的样本，见第6节）
- LV6000619917
- LV6000620016

其余样本仍在串行队列（`.../series/`）里陆续跑。**besthit 阶段不需要等全部
mapping 完成**——`submit_oryza_besthit.sh submit all` / `local all` 会自动只
处理当前已有 `.finished` mapping BAM 的样本，可以随 mapping 进度反复重跑
（已完成 besthit 的样本会被跳过，见第5.4节）。

## 5. Best-hit / aDNA 损伤校正过滤（本次工作的主体）

### 5.1 为什么没有直接复用 aeDNA_popgen 的通用分类器

`aeDNA_popgen/01_TaxAssign/workflow/scripts/besthit_competitive.py` 是一个
通用的 ngsLCA 风格层级分类器：沿谱系从叶子往根走，每层比较"子树内最优 NM"
vs"子树外最优 NM"，赢的最深节点即为归属。这套逻辑**不区分古 DNA 末端损伤和
真实突变**，也不是"Oryza vs 非 Oryza 二元竞争"这个具体需求。用户明确要求
（原话）：

> 参考代码不能直接照搬，因为我们需要增加古DNA末端损伤校正，并且只判断
> Oryza 与非 Oryza 的竞争关系。

所以这里是专门重写的一套逻辑，核心是：**扣除末端损伤后的 `adjusted_NM`，
Oryza 最优 hit 必须与非 Oryza 最优 hit 打平或更好，才 KEEP**。

脚本落地：
- `scripts/oryza_besthit/oryza_besthit_damage_filter.py` —— 主统计脚本（约550行）
- `scripts/oryza_besthit/submit_oryza_besthit.sh` —— SLURM 提交脚本（约430行）
- `scripts/oryza_besthit/rebuild_all_wgs_asian_irgsp_acc2taxid.sh` —— acc2taxid
  三源重建+taxid错标修正脚本（见0.5节，2026-08-08新增，取代已删除的
  `fix_asian_rice_panel_taxid.sh`）

**2026-08-08补充：与师兄（同实验室）的类似脚本对比**（本仓库外，
`/home/zxr619/...`路径下，未纳入本项目）——师兄有一套两阶段架构：
`besthit_competitive_top10_showOryza_optimized.py`（通用ngsLCA风格全谱系
分类器，带损伤校正但`--damage-end-bases`默认**1bp**，比我们的5bp窄很多，
另外还有我们没有的质量预筛`--target-min-sim 95.0`/`--target-max-nm 2`）+
下游`select_oryza_competitive_reads.py`（从上一步的top10表里做真正的KEEP
筛选，`.smk`里配置`MAX_DIFFERENCE=0`，只认sativa/rufipogon两个种，
**不含nivara**，比我们当时`--oryza-taxids`默认的3个种少一个）。他的
competitive mapping数据库组成从给到的两个文件里**无法确认**（目录命名
线索是"rice_panel"+"k100_N1"的bowtie2风格参数，看起来跟我们用WGS+亚洲
panel的思路接近，但没有直接证据）。这个对比当时只在对话里讨论过，没有
写进代码——下面5.1b节是根据这个对比做的实际重写。

### 5.1b ⚠️2026-08-08 v2重写：Oryza范围改成全genus，取代3-species白名单

**用户明确要求（原话）**："不要像师兄那个只看sativa和rufipogon，我们的是
里面oryza genus的全部species全部留下来比较看就可以"——这是v2唯一的硬性
需求，其余"综合师兄脚本能调整什么"由我判断取舍，取舍理由都写在下面。

**核心改动**：`--oryza-taxids`不再默认硬编码`"4529 4530 4536"`(rufipogon/
sativa/nivara 3个种)，改为默认**动态解析整个Oryza属**——启动时从
`nodes.dmp`里找到taxid 4527(Oryza属)下面**全部rank=species的后代**(新增
`Taxonomy.genus_species_taxids()`，用children_map做BFS，只遍历Oryza属自己
的子树，不用扫全部taxonomy)，这个集合就是新的白名单。`--oryza-taxids`
仍然保留，作为**手动override**——显式传参数时完全按v1那样只用给定的taxid，
不做全genus解析，用来复现v1的窄范围行为。

**为什么这个改动本身就是在修一个真实的误判风险，不只是"更全面"**：旧版
`--oryza-taxids`默认只有3个种，其余~15个Oryza属物种(比如数据库里contig数
不少的*O. longistaminata* 905条、*O. coarctata* 450条，见6.1节)**不在白
名单里**——这意味着一条read如果真的是*O. longistaminata*，但数据库里同时
也命中了某个真正的非Oryza物种，这条read会被当成"非Oryza竞争者"去跟
sativa/rufipogon/nivara的最优hit比，如果sativa/rufipogon/nivara都没命中
（这条read其实压根不是这3个种），就会直接落进`best_oryza is None`分支，
判定`REJECT/no_oryza_hit`——**一条真正的Oryza属read被当成"没有Oryza命中"
拒绝掉了**。这不是理论担忧，本地用合成数据复现过这个场景（见下方"本地
验证"），v1语义下确实会把longistaminata-only的read错误REJECT。genus-wide
之后，`best_oryza`池子覆盖整个属，这类read会被正确识别为Oryza并纳入
adjusted_NM比较。

**综合师兄脚本，还考虑过什么，取舍理由**：

| 想法来源 | 采纳/不采纳 | 理由 |
|---|---|---|
| genus-wide范围（用户明确要求） | ✅采纳，核心改动 | 见上 |
| MD tag字符串直接解析(不用`pysam.get_aligned_pairs(with_seq=True)`) | ❌**不采纳**（本次） | 师兄这套快速实现有一个连带的"只对每个species的raw-NM最小值那几条alignment做MD解析"的剪枝——这个剪枝假设raw NM排序等于damage校正后的adjusted NM排序，**这个假设不总成立**：一条raw NM更高的alignment，如果它的额外错配恰好落在末端损伤窗口内，扣除损伤后反而可能比raw NM更低的alignment adjusted_NM更小。我们的`--damage-window`默认5bp/端(最多10个信用)，比师兄的1bp/端(最多2个信用)宽得多，这个"raw NM差距被损伤信用反超"的风险对我们更大。为避免在没有服务器真实数据比对验证的情况下引入静默的准确性回归，`alignment_metrics()`保持v1原样(对每条alignment都完整计算，不做raw-NM预剪枝)。如果以后单样本运行时间真的成为瓶颈，值得重新评估，但要先跟v1的输出做过并排验证 |
| `--damage-end-bases`默认1bp | ❌不采纳 | 用户没要求改，我们的5bp已经过4样本留存率验证("5bp暂时可用"，见7.5节)，不是这次改动范围 |
| `--target-min-sim`/`--target-max-nm`质量预筛 | ⚠️采纳但**默认关闭** | 新增`--min-best-similarity`/`--max-best-raw-nm`两个可选参数(OR逻辑，任一达标就放行，跟师兄脚本一致)，在per-species分类之前先看这条read全局最优raw-NM是否够格。默认不传(两者都是`None`)时完全不生效，不影响已验证的行为；想尝试的话可以显式加上，但目前对我们数据没有验证过阈值该设多少 |
| `select_oryza_competitive_reads.py`的`MAX_DIFFERENCE=0`(只认打平) | ❌不采纳 | 我们的判定本来就是`adjusted_NM(oryza) <= adjusted_NM(nonoryza)`，等价于打平或更好，逻辑已经一致，不需要额外参数 |
| 完整ngsLCA式全谱系归类(assign到最深胜出节点) | ❌不采纳 | 我们是"Oryza vs 非Oryza"这个具体的二元竞争需求(见5.1节用户原话)，不是通用谱系分类问题，不适合套用 |
| multiprocessing worker pool加速 | ⏸️推迟，未实现 | 有价值但改动面大、这次没有服务器条件实测验证，先不引入，记在这里供以后需要提速时参考 |

**已从仓库归档**：`scripts/oryza_besthit/oryza_besthit_damage_filter_v1.py`
——`git mv`保留完整历史，不是删除。`submit_oryza_besthit.sh`同步更新：
`ORYZA_TAXIDS`环境变量默认改成空字符串(触发genus-wide)，新增
`ORYZA_GENUS_TAXID`(默认`4527`)、`MIN_BEST_SIMILARITY`/`MAX_BEST_RAW_NM`
(默认都为空=关闭)。

**输出格式变化**：`summary.tsv`新增一列`rejected_low_quality`(质量预筛
关闭时恒为0)。⚠️**这意味着v1和v2产出的`summary.tsv`列数不一样**——已经
用v1跑完的4个样本(LV6000619499/619917/620016/620032)的`summary.tsv`是
8列(无`rejected_low_quality`)，v2产出的是9列。`submit_oryza_besthit.sh
merge`直接`head+tail`拼接多个样本的`summary.tsv`，**不要把v1和v2的
summary.tsv混在一起跑merge**，会列错位。

**本地验证**（本机Mac，无服务器数据/无真实BAM时做的，隔离venv装pysam
0.24.0）：手工构造了一个最小合成taxonomy(genus Oryza=4527，4个种：
rufipogon/sativa/nivara/longistaminata，加一个非Oryza外群物种)+对应BAM，
3条测试read，重点验证：
- **genus-wide默认模式**：一条只命中`longistaminata`(不在v1白名单里)和
  外群物种的read，`longistaminata`adjusted_NM更优——v2正确判定
  `KEEP/oryza_at_least_as_good`，`best_oryza`正确解析成
  `Oryza longistaminata`(taxid 4528)。启动日志正确打印出解析到的4个属
  内物种。
- **`--oryza-taxids`手动override**：显式传`4529 4530 4536`时，同一条
  longistaminata-only read被判`REJECT/no_oryza_hit`，且
  `best_nonoryza`栏位显示的正是`Oryza longistaminata`——**这正是上面
  描述的误判风险的实际复现**，证明override确实完整复现了v1语义，也
  证明了genus-wide改动解决的是一个真实存在、可复现的问题，不只是
  理论上更完备。
- sativa-vs-外群的baseline read在两种模式下都正确`KEEP`，行为不受影响。
- 质量预筛(`--min-best-similarity 99`)开启时，3条测试read(相似度都在
  93-97%之间)全部被正确判`REJECT/low_quality_pregate`，一致性校验
  (`input_reads == kept + rejected_nonoryza_better + rejected_no_oryza +
  rejected_low_quality + unclassified_reads`)照常通过。
- **没有验证过的地方**：这一切都是合成数据，没有在真实BAM(131个数据库
  的competitive mapping产出)上跑过，也没有拿v1已经产出的4个样本真实结果
  跟v2重跑做过并排对比。genus-wide范围扩大后，实际retention rate(留存率)
  大概率会变化(至少一部分之前被错误REJECT的genus-内非focal-3-species
  read现在会被KEEP)，具体变化多少不知道。

**待办（下一个session）**：
1. 下载新版本三个脚本(见下方下载命令)，先跑`check`确认环境正常
2. 建议先对已经跑过v1的4个样本之一做`smoke`(1000-read子集)，对比v2和
   v1的KEEP/REJECT分布差异，确认genus-wide改动的实际影响幅度符合预期
   (不应该暴增，因为smoke test的分析显示"惜败"read的头号竞争对手基本都是
   跟Oryza毫不相关的物种，genus内其他种在top20里几乎没出现过，见6节)
3. 确认没问题后，可以考虑重新跑一遍这4个已完成的样本(删除`.finished`
   标记后重跑，或者直接跑到新的`OUT_DIR`保留两份对比)，让"5bp可用"这类
   已有结论建立在v2的口径上，避免v1/v2结果混用造成解读混乱
4. 剩余样本直接用v2跑(`submit all`/`local all`默认行为已经是genus-wide)

**下载命令**：
```bash
cd /home/scratch/yinmt202607/gene/scripts
curl -fsSL -O https://raw.githubusercontent.com/Inmpain/rice_adna_pipeline/codex/oryza-competitive-mapping/scripts/oryza_besthit/oryza_besthit_damage_filter.py
curl -fsSL -O https://raw.githubusercontent.com/Inmpain/rice_adna_pipeline/codex/oryza-competitive-mapping/scripts/oryza_besthit/oryza_besthit_damage_filter_v1.py
curl -fsSL -O https://raw.githubusercontent.com/Inmpain/rice_adna_pipeline/codex/oryza-competitive-mapping/scripts/oryza_besthit/submit_oryza_besthit.sh
chmod +x submit_oryza_besthit.sh
bash submit_oryza_besthit.sh check
```

### 5.2 核心算法

每条 alignment 计算：

| 字段 | 含义 |
|---|---|
| `NM` | bowtie2 的 NM tag（原始编辑距离，含替换+indel） |
| `substitution_count` | 从 MD/CIGAR 精确统计的真实替换数（不含 indel），靠 `pysam.get_aligned_pairs(with_seq=True)` 拿到逐位点 (ref_base, read_base)，lowercase ref_base = 该位点是 mismatch |
| `terminal_damage_count` | 只在读长两端 `--damage-window`（默认5）范围内、且**换算回 read 自身 5'→3' 方向后**符合古 DNA 损伤特征（5'端 C→T，3'端 G→A）的替换才计入 |
| `adjusted_NM` | `NM - terminal_damage_count`（indel 永远不扣） |

**reverse-strand 的方向换算是全脚本最容易写错的地方**：BAM 里 SEQ 字段已经是
原始 read 的反向互补（SAM 规范），要拿到"read 自身 5'→3' 方向"，需要把
reference-forward 顺序的 (ref_base, read_base) 序列**倒序 + 互补**两步一起做——
只倒序不互补，会把参考链上的 G→A 信号错当成"3'端 G→A"留在错误的位置，
从而漏判或错判。已经用手工构造的 reverse-strand 测试 BAM 验证过这一步
（见5.3）。

按 species taxid 分组（BAM reference name → acc2taxid → 沿 nodes.dmp 上溯到
最近的 `rank=="species"` 祖先），每个 species 只保留最优 alignment，排序
优先级：`adjusted_NM` 最小 → `NM` 最小 → `substitution_count` 最小 → `AS`
最大 → reference 名称（稳定/可重复）。非 Oryza species 取前 `--top-n`
（默认10）写入审计表；Oryza（`--oryza-taxids`，默认4529/4530/4536）**不占
非 Oryza 的10个名额**，命中就必录。

判定（`best_nonoryza` / `best_oryza` 各自是排名第一的非Oryza/Oryza species）：

1. 没有 Oryza hit → **REJECT** (`no_oryza_hit`)
2. 有 Oryza hit、没有非 Oryza hit → **KEEP** (`oryza_only_no_competitor`)
3. `best_oryza.adjusted_NM <= best_nonoryza.adjusted_NM` → **KEEP** (`oryza_at_least_as_good`)
4. 否则 → **REJECT** (`nonoryza_better`)

不用 MAPQ（bowtie2 `-k 100` competitive mapping 下不可靠，`asian_rice_panel/
IRGSP` 那种多参考竞争场景下 MAPQ 基本是固定值，没有信息量）。BAM 没有去重、
也**不需要**在这一步去重——同一 read 的多行来自131个数据库 × `-k 100`，不是
PCR duplicate，去重逻辑放在更下游。

### 5.3 输出

每样本，全部 `.tmp` 写完原子 `os.rename`；**只有一致性检查通过才写
`<sample>.finished`**：

| 文件 | 内容 |
|---|---|
| `<sample>.besthit.top10_species.tsv.gz` | 审计表，一条 read 对应多行（非Oryza top10 + 必录的Oryza行），列：`read_name read_length species_hit_count top10_rank is_always_included species_taxid species_name reference_name NM substitution_count terminal_damage_count adjusted_NM AS` |
| `<sample>.oryza_filter.decisions.tsv.gz` | 每 read 一行的最终判定，列：`read_name best_nonoryza_taxid best_nonoryza_name best_nonoryza_NM best_nonoryza_damage best_nonoryza_adjusted_NM best_oryza_taxid best_oryza_name best_oryza_NM best_oryza_damage best_oryza_adjusted_NM decision reason` |
| `<sample>.besthit_oryza.fastq.gz` | 从候选 FASTQ 按 KEEP read name 抽出，每条只输出一次，序列/质量值原样保留 |
| `<sample>.summary.tsv` | 单样本一行：`sample input_reads reads_with_alignment reads_with_oryza_hit kept_reads rejected_nonoryza_better rejected_no_oryza unclassified_reads`，脚本内部自查 `input_reads == kept + rejected_nonoryza_better + rejected_no_oryza + unclassified_reads`，不一致会非零退出（但这个文件本身仍会被写出，方便debug；`.finished` 不会被写） |

**`unclassified_reads` 的定义要注意**：包含两类，(a) 候选 FASTQ 里的 read 但
BAM 里完全没有任何 alignment（bowtie2 `--no-unal` 直接没输出，这类 read
**不会**出现在 `decisions.tsv` 里，只能从 `summary.tsv` 的
`input_reads - reads_with_alignment` 间接算出数量，查不到具体是哪些
read）；(b) BAM 里有 alignment，但每一条都解析不到 taxid/species（会在
`decisions.tsv` 里留一行 `decision=UNCLASSIFIED, reason=no_resolvable_species`，
可以按名字查到）。

多样本汇总 `besthit_summary.tsv` **不**由 Python 脚本直接写（避免并行 SLURM
作业竞争同一个文件），而是 `submit_oryza_besthit.sh merge` 事后拼接每个样本
自己的 `<sample>.summary.tsv`。

### 5.4 用法（`submit_oryza_besthit.sh` 的5种模式）

```bash
cd /home/scratch/yinmt202607/gene/scripts

# 1. 校验路径 / python3+pysam / SLURM partition，不提交任何作业
bash submit_oryza_besthit.sh check

# 2a. SLURM 队列：单个样本前1000条 reads 做 smoke test
#     （独立输出目录 besthit/smoke_test/，不写 .finished，不影响正式产出）
bash submit_oryza_besthit.sh smoke LV6000619499

# 2b. 前台本地跑（不进 SLURM 队列，适合调试/小样本快速出结果；量产还是用3a/3b）
bash submit_oryza_besthit.sh local LV6000619499
bash submit_oryza_besthit.sh local all      # 全部已完成mapping的样本，依次前台跑

# 3a. SLURM 队列：全量跑（自动跳过 besthit/<sample>.finished 已存在的样本）
bash submit_oryza_besthit.sh submit LV6000619499 LV6000619917 LV6000620016

# 3b. 等价写法：自动发现所有已完成 mapping 的样本（可随 mapping 进度反复重跑）
bash submit_oryza_besthit.sh submit all

# 4. 所有作业跑完后，拼接每个样本的 summary.tsv
bash submit_oryza_besthit.sh merge
# -> .../oryza_competitive_mapping/besthit/besthit_summary.tsv
```

可覆盖的环境变量（默认值见脚本头部）：`BAM_DIR` `FASTQ_DIR` `ACC2TAXID`
`NODES` `NAMES` `ORYZA_TAXIDS`（默认 `"4529 4530 4536"`）`DAMAGE_WINDOW`
（默认5）`TOP_N`（默认10）`OUT_DIR` `SLURM_PARTITION`（默认`comp`）
`JOB_CPUS`（默认4）`JOB_MEM_MB`（默认16000，**未实测校准，见第7节待办4**）
`JOB_TIME`（默认`04:00:00`）。

若 `check` 报 `import pysam` 失败，需要先在这个 `python3` 所在环境装 pysam。

### 5.5 已完成的本地验证（本机 Mac，无服务器数据/无sbatch时做的）

用 `samtools view -bS` + `sort -n` 手工构造了一个最小 BAM（5条参考、5条候选
read，覆盖 KEEP-打平/KEEP-oryza独占/REJECT-惜败/UNCLASSIFIED/BAM里完全没有
alignment 这5种情况），在隔离 venv 里装 pysam 实测跑通全流程，重点验证：

- **reverse-strand 末端损伤校正确实会反转最终判定**：构造了一条 read，原始
  NM 下"非 Oryza 更优"应判 REJECT，扣除末端损伤后 `adjusted_NM` 打平，正确
  判成 KEEP（`oryza_at_least_as_good`）——这是本脚本存在的核心原因，是最关键
  的一条验证。
- `input_reads == kept + rejected_nonoryza_better + rejected_no_oryza +
  unclassified_reads` 自检通过；FASTQ 抽取数等于 `kept_reads`。
- `--limit-reads`（smoke test）路径确认不写 `.finished`。
- 重复跑同一个样本（模拟 SLURM 任务失败重提交）不会因为 `.tmp` 残留或
  rename 冲突出错（幂等）。
- `submit_oryza_besthit.sh run`（sbatch job body）、`local`、`merge` 三种
  调用路径都跑通。

### 5.6 Bug 修复历史（服务器实跑中发现，均已修复并推送）

| # | 问题 | 现象 | 修复 | commit |
|---|---|---|---|---|
| 1 | `module load python/ 2>/dev/null` 只重定向了 stderr，但这台集群的 `module` 命令把错误信息写到了 stdout | `check`/`smoke` 输出里混入 `ERROR: Unable to locate a modulefile for 'python/'`（无害噪音，不影响功能） | 改成 `>/dev/null 2>&1` | `a1738fc` |
| 2 | **最关键的一个**：`sbatch` 会把提交的脚本复制到每个作业专属的 spool 目录（`/var/spool/slurm/d/job<id>/`）并执行那份拷贝，不是原始文件；脚本内部用 `readlink -f "${BASH_SOURCE[0]}"` 推导自己所在目录来找同目录下的 `oryza_besthit_damage_filter.py`，这个推导在 SLURM job 里跑的时候会指向错误的 spool 目录 | smoke test 实际报错：`python3: can't open file '/var/spool/slurm/d/job1807625/oryza_besthit_damage_filter.py': No such file or directory` | 提交时把已经在提交端正确解析好的 `PY_SCRIPT` 绝对路径通过 `--export="ALL,PY_SCRIPT=${PY_SCRIPT}"` 显式传进 job 环境，job 内不再重新推导 | `a1738fc` |
| 3 | 无 `submit all` / `local` 模式，用户批量跑/本地调试不方便 | — | 新增 `submit all`（自动发现已完成mapping的样本）和 `local`/`local all`（前台不走 sbatch，复用同一套 worker + `.finished` 跳过逻辑） | `83a0f6f` |
| 4 | `asian_rice_panel.acc2taxid` **整份文件** 4529/4530 系统性反标（最初误判为只有`np7`等12个基因组按名字错标，用户截图纠正为全文件二元对调），且`all_wgs_asian_irgsp.acc2taxid`的三源合并配方一度以为已丢失 | 用户先用染色体MD5比对发现np7错标；后贴出文件截图证明是整文件反标（`G25_ruf_W*`标成4530、`IRGSP-1.0_genome`/`X24_kas`标成4529）；用户在服务器`taxonomy_CPH`目录下找到了三源合并的真实配方 | 见0.5节，`rebuild_all_wgs_asian_irgsp_acc2taxid.sh`完整重建+全文件4529/4530互换（取代已删除的`fix_asian_rice_panel_taxid.sh`补丁式方案和基于12基因组误判的第一版重建脚本）——**服务器dry-run已验证方向正确**，见0.6节 | `8ba66b4`（全文件互换，服务器验证过）、`3edf2dbc`（第一版重建，已废弃）、`b586778`（删除最早补丁脚本） |
| 5 | `irgsp.acc2taxid`有一行字面表头(`accession/accession.version/taxid`)，如果原样`cat`进合并文件，第3列会是非数字字符串`"taxid"`，besthit脚本解析taxid时会出问题——服务器真实dry-run才暴露出来，本地合成数据测试没覆盖这个场景 | 服务器dry-run输出里`irgsp.acc2taxid`第3列分布显示`1 taxid`（应该全是数字taxid） | 见0.6节，新增Step 2b检查三个源文件里第3列是否为纯数字，`--apply`阶段自动过滤非数字taxid行——**只在本地合成数据测试过，还没在服务器重新跑dry-run确认**，是明天第一件事 | `2965348` |

**重要经验**：bug #2 这一类"脚本自我定位"的坑，本地测试完全测不出来（本机
没有 sbatch，只能测 `run`/`local` 这两个不依赖 `sbatch` 拷贝脚本的路径），
必须在真实 SLURM 环境里才会暴露。以后改这个脚本、涉及"脚本找自己同目录下的
文件"这类逻辑时要格外小心，优先考虑用 `--export` 显式传值，而不是依赖
`$0`/`BASH_SOURCE` 在 job 里重新推导。

## 6. 已经拿到的实际结果：smoke test 分析（LV6000619499，前1000条 reads）

```
sample        input_reads  reads_with_alignment  reads_with_oryza_hit  kept_reads  rejected_nonoryza_better  rejected_no_oryza  unclassified_reads
LV6000619499  1000         1000                  956                   103         853                       44                 0
```

decision/reason 分布：
```
853 REJECT  nonoryza_better
100 KEEP    oryza_at_least_as_good
 44 REJECT  no_oryza_hit
  3 KEEP    oryza_only_no_competitor
```

**REJECT(nonoryza_better) 的输赢差距分布**（`best_nonoryza_adjusted_NM -
best_oryza_adjusted_NM`，负值=Oryza更差）：
```
差1: 541 (63%)   差2: 288 (34%)   差3: 23   差4: 1
```

**KEEP 的差距分布**：
```
打平(差0): 98个   Oryza真正赢(差1): 2个
```

**REJECT(nonoryza_better) 里排名第一的非Oryza物种 top20**（taxid → 学名，
已用 names.dmp 逐条核实）：
```
1711249  37  Eucalyptus dawsonii        桉树
145626   31  Bradybaena similaris       陆生蜗牛
519541   29  Geodia barretti            深海海绵
13415    26  Chamaecyparis obtusa       扁柏(针叶树)
223100   23  Orobanche coerulescens     列当(寄生植物)
2790670  19  Closterium sp. NIES-67     新月藻(藻类)
126911   19  Ammopiptanthus mongolicus  沙冬青(灌木)
103762   18  Zizania palustris          北美菰(野生稻近缘属，唯一说得通的)
4682     13  Allium sativum             大蒜
322858   13  Diplosoma virens           海鞘(海洋无脊椎动物)
3077     12  Chlorella vulgaris         小球藻
3759     11  Prunus yedoensis           樱花
192012    9  Mikania micrantha          薇甘菊(藤本)
9872      8  Odocoileus hemionus        骡鹿(哺乳动物)
56866     8  Luffa acutangula           丝瓜
1709936   8  Eurotiomycetes sp.         真菌
988163    7  Tiliacea citrago           蛾类
641091    7  Trididemnum clinides       海鞘
980011    6  Bidens hawaiensis          鬼针草
39329     6  Lavandula angustifolia     薰衣草
```
这份名单在生物学上没有一致性——深海海绵、两种海鞘、陆生蜗牛、北美骡鹿，跟
东南亚考古水稻沉积样本对不上号，除了 *Zizania palustris*（野生稻近缘属）之外
基本看不出真实污染的规律。**这是判断"惜败"read性质的关键证据**：如果真是
参考基因组缺失导致的误判，输给的应该是生物学上说得通的近缘种；现在这个
名单更像是"短片段靠概率撞上129个shard里某个不相关物种的某段序列"，见6.2。

### 6.1 acc2taxid 覆盖度诊断结果（第7.1节问题，已有答案）

⚠️**2026-08-08更新：下表统计于`asian_rice_panel.acc2taxid`整文件taxid反标
(见0.5节)被修正之前，需要在服务器跑完`rebuild_all_wgs_asian_irgsp_acc2taxid.sh
--apply`后重新统计**。这次修正的范围比最初以为的大很多——不是只有np7的16行，
而是`asian_rice_panel.acc2taxid`**全部**行的4529/4530标签对调（panel里全部
sativa代表基因组此前很可能被计入了下面的4529/rufipogon桶，全部野生rufipogon
基因组则可能被计入了4530/sativa桶）。**下面这张表在修正前的具体数字（包括
"sativa 111条"和"rufipogon 717条"这两个关键数字本身）都不可信，"sativa覆盖度
偏低"这个结论需要重新统计后才能确认是否成立、成立到什么程度——不能再假设
方向不变**。

Oryza 属18个已知种在 `all_wgs_asian_irgsp.acc2taxid` 里的 contig 数（用
`nodes.dmp` 里 parent==4527 找到全部种级taxid，逐个数 acc2taxid 里精确匹配
的行数；*sativa* 另外加了亚种级taxid 39946/39947/1050722/1080340/2998809
的计数，只有39947=japonica有15条，其余0条）：

```
4528  Oryza longistaminata    905
4529  Oryza rufipogon         717   <- 野生稻，覆盖很好 (⚠️待修正后重新统计)
4530  Oryza sativa             96 (+15 japonica = 111)  (⚠️待修正后重新统计)
4532  Oryza australiensis      60
4533  Oryza brachyantha        60
4534  Oryza latifolia          78
4535  Oryza officinalis        91
4536  Oryza nivara             26   <- 野生稻，有但不算多
4537  Oryza punctata           61
4538  Oryza glaberrima         60  (非洲栽培稻)
40148 Oryza glumipatula        47
40149 Oryza meridionalis       12
52545 Oryza alta               24
63629 Oryza minuta            189
65489 Oryza barthii            51
77588 Oryza coarctata         450
110451 Oryza schlechteri      119
127571 Oryza malampuzhaensis  198
```

**结论：这个假说被证伪**（⚠️此结论基于修正前、已知有整体反标问题的数字，
**需要重新统计后重新确认，不能再假设方向不变**——见上方警告）。原始数字
显示*sativa* (111条) 不但不缺，反而是覆盖度偏低的那一档——*rufipogon*
(717)、*longistaminata*(905)、*coarctata*(450) 等好几个种覆盖度都远高于
sativa；"数据库里野生稻/近缘种缺失"不是"差1惜败"这个模式的解释，第8.1节
原本"优先扩库"的建议已经不成立，改成"优先查损伤窗口"（见6.2、7.5、8.1）。
**但既然asian_rice_panel的sativa/rufipogon计数本身可能是对调的，"证伪"这个
结论也需要拿修正后的真实数字重新验证一遍，不是自动成立**。

### 6.2 差1惜败的read为什么会输：不是短，是长（新发现，推翻了最早的猜测）

最早怀疑"惜败"read是不是偏短（短片段信息量少、更容易被噪声压过），实测
结果**方向相反**：

```
                均值(bp)  中位数(bp)  n
差1惜败(REJECT)   ~68.5      ~68      541
KEEP(打平守住)    ~62.3      ~61      103
```

（原始 `top10_species.tsv.gz` 是一条read多行的审计表，第一次按行统计长度会
被"这条read命中了多少个物种"放大，必须先按read_name去重再统计——这是个
真实踩过的坑，两组的原始"行数"分别是541条read共约6778行、103条read共约
1225行，去重后总数才分别对上541/103，接手时如果要重跑这类分析记得先去重）

**长而不是短的read更容易惜败，机制上的解释**：`--damage-window` 固定只看
读长两端各5bp，超出这个范围的任何替换（不管是没修完的非末端损伤、真实的
生物学SNP、还是测序错误）都没资格被 `adjusted_NM` 扣除，会原样拉高分数。
读长越长，两端5bp加起来占全长的比例越小，"没资格被校正"的中段就越大，
累积一个额外错配、从而在adjusted_NM上被非Oryza竞争者反超的概率也就越高。
这个方向和实测数据完全吻合。

**这个假说还差一步验证**：损伤窗口5bp这个默认值有没有校准过？如果这批样本
的真实损伤信号（末端C→T/G→A频率）衰减到背景水平的距离明显超过5bp（比如
mapDamage风格分析常见的10-15bp甚至更远，取决于样本降解程度），那5bp的窗口
本身就设窄了，会系统性地对长读长不利——这不是"要不要加margin"的问题，是
损伤校正窗口参数本身需要重新校准。有没有现成的损伤曲线分析（不管是这批
angkor样本还是这个项目别的样本）？如果没有，可以从BAM直接算一个粗略版本
（末端N个位置的C→T/G→A频率 vs 距末端距离），不需要跑完整mapDamage2。
见7.5的具体诊断命令。

## 7. 待确认的开放问题

### 7.1 ~~acc2taxid 里 Oryza 属各物种的实际覆盖度~~ 已回答，见6.1

结论：Oryza属18个种普遍有几十到几百条contig，*sativa*(111条)不缺，反而是
偏低的一档。⚠️2026-08-08：这个"已解决"的结论建立在有整体taxid反标问题的
旧数字上（见0.5/6.1节），**修正后需要重新核实，不能预设方向不变**——严格
说这个问题应该重新打开，等服务器重新统计6.1节的表之后再关闭。

### 7.2 `db/3k/wild/` 140+ 野生稻组装的物种身份

```bash
ls /home/scratch/yinmt202607/db/3k/wild/ | head -20
ls /home/scratch/yinmt202607/db/3k/wild/ | wc -l
ls /home/scratch/yinmt202607/db/3k/       # 找有没有说明文档/manifest
head -3 /home/scratch/yinmt202607/db/3k/tmp/3kall_variety_map.tsv
for f in $(ls /home/scratch/yinmt202607/db/3k/wild/ | head -5); do
  sid="${f%%.*}"
  echo "=== $sid ==="
  grep -w "$sid" /home/scratch/yinmt202607/db/3k/tmp/3kall_variety_map.tsv
done
```
**给了用户、还没有结果**。这份 `file_path.md` 里只有一行目录清单提及
（`db/3k/wild/{SampleID}.transfer.merge.chr.fasta  # 140+野生稻/近缘种染色体级
组装基因组`），仓库里没有任何地方写明这140+个样本具体是哪个物种/亚群——
需要用 `3kall_variety_map.tsv`（3K RGP官方样本元数据表，含
VARIETY_INDEX/NAME/IRIS_ID）交叉查询才能确认。**注意**：3K Rice Genome
Project官方主体是3024份**栽培稻**（indica/aus/aromatic/japonica等亚群），
不是野生稻——这批 `wild/` 目录既然单独命名和别的3K数据分开放，大概率是
额外补充的野生近缘种outgroup，但没有manifest佐证之前不要假设。**另见
`docs/RESEARCH_ROADMAP.md`第2节C——这140+个组装是否等同于Guo et al. 2025
pangenome论文提供的145个组装，是main分支也在追的同一个开放问题。**

**优先级下调**：这条本来是基于"数据库野生稻缺失"假说排的第一优先级，
现在7.1已经证明acc2taxid里Oryza属覆盖面本身就很宽，这条不再紧急，可以
按需做，不是阻塞项。

### 7.3 `--oryza-taxids` 是否真的是 species rank

```bash
grep -P '^4529\t|^4530\t|^4536\t' /home/database/ref20250728/taxonomy_CPH/ncbi/20250530/nodes.dmp
```
脚本启动时只校验这三个 taxid 存在于 nodes.dmp，**不校验 rank 是不是
species**——如果 asian_rice_panel/IRGSP 的 acc2taxid 把 contig 指到了某个非
species 的 Oryza 内部节点，`species_of()` 会往上找到别的 species 节点甚至
找不到。还没确认过。

### 7.4 MD tag 覆盖率 / 实际内存占用

Smoke test 日志里如果有 `[warn] N alignments had no usable SEQ/MD` 这一行，
需要看 N 的数值——用户贴过的 smoke test 输出片段里没有看到这行（只贴到
`[load] acc2taxid...` 就中断了，后续 `[warn]`/`[summary]` 完整输出没有再贴
过一次完整版），**没有正式确认过 MD tag 覆盖率是否有问题**，只是根据最终
1000条里956条有 alignment 走完全流程、没有报错退出，间接推测大概率没问题。
接手时建议重新要一次完整 log。

内存：`JOB_MEM_MB` 默认16000（16GB）是本次工作开始时凭经验猜的，从未用
`seff <job_id>` 或 `wc -l` nodes.dmp/names.dmp/acc2taxid 校准过实际占用——
理论上 taxonomy dump（几百万taxid，3个dict）+ acc2taxid dict 全量常驻内存
是这个脚本的内存大头，量级可能到GB级别，具体数字取决于这两个文件的真实大小，
需要 `seff` 或 `wc -l` 实测。

### 7.5 `--damage-window=5` 是否偏窄——需要实测损伤衰减曲线

见0/6.2的推理：差1惜败的read系统性比KEEP的read长，机制上最可能的解释是
固定5bp窗口对长读长不利。这个假说要成立，前提是这批样本的真实损伤信号确实
延伸超过5bp。**先问用户/查项目历史**有没有现成的损伤曲线分析（mapDamage2
或类似工具跑过的5'C→T / 3'G→A频率-vs-位置图，针对angkor这批样本或者项目
里任何同批次测序的样本）。如果没有，可以从besthit的smoke test BAM直接估算
一个粗糙版本（只需要 KEEP 结果里"确认是水稻"的read，避免非水稻read混进来
干扰损伤信号统计）：

```bash
# 思路：besthit_oryza.fastq.gz 里的read名字，去原始 by_sample BAM 里找到它们
# 的primary alignment，统计末端N个位置(比如N=15，覆盖比默认5bp更远的范围)
# 的错配是不是C->T(5'端)/G->A(3'端)富集，并且看这个富集程度随距离衰减到
# 背景水平大概在第几个bp——这个可以另写一个小脚本，不在
# oryza_besthit_damage_filter.py现有功能范围内，需要新写。
```
如果衰减距离明显超过5bp（比如10-15bp），建议把 `--damage-window` 调大，
重新跑 smoke test 看"差1惜败"的541条里有多少会因此翻盘成KEEP，再判断这
是不是6.2这个模式的完整解释，还是只能解释一部分。

**2026-08-07进展**：4个样本(LV6000619499/619917/620016/620032)已跑完全量
besthit，留存率(kept/input)分别为5.00%/5.00%/4.37%/4.92%，四个样本高度
一致，没有出现某个样本被5bp窗口明显"冤枉"的迹象——这是"5bp目前没出大问题"
的数据支撑，用户据此判断"5bp暂时可用"，倾向于不再深挖这个问题，转向
确认参考基因组体系（0.5/7.2节）和ecotype-pca-panel分支的标签来源问题。
详见`docs/RESEARCH_ROADMAP.md`。

**参考对照**：师兄的类似脚本`--damage-end-bases`默认只有1bp（见5.1末尾的
对比），比我们的5bp窄得多——如果他那边也在类似样本上跑出了合理结果，可能
说明5bp已经偏宽而不是偏窄；但这只是间接旁证，他的数据库/样本降解程度都
和我们不一定可比，不能直接当结论用。

## 8. 后续三条工作线（讨论过，还没定最终方案）

这三条是用户提出的下一步方向，讨论了可行性和优先级，但都**还没有开始实际
写代码/跑分析**，需要先决定做哪个/什么顺序。

### 8.1（优先级已下调，见7.1）Oryza 参考基因组数据库整理/扩充

**这条原本排第一优先级，前提是"数据库野生稻缺失"——这个前提已经被7.1的
实测数据证伪，所以不再是紧急项**，但参考基因组体系本身"还是很乱"（用户
原话）这个观察依然成立，值得找机会理一遍，只是不再挡着best-hit这条线：

- `db/16/`（资源组A，NCBI datasets 16基因组+Liftoff注释）和 `db/3k/`（资源组B，
  3K Rice Genome Project数据）是**两个不同来源**，关系没有理清。
- `asn720data/`（720份现代/近现代品种PLINK面板）跟16个angkor古代样本的关系
  "尚未最终确认"（是否包含、是否capture panel来源）——⚠️2026-08-07更新：
  `ecotype-pca-panel`分支发现`asn720data/asn720.pop.fam`的FID列是`OrA-OrF`
  群体标签的关键来源，已不再是"可以忽略"的旧数据，详见
  `docs/RESEARCH_ROADMAP.md`。
- `db/3k/wild/` 这140+个野生稻/近缘种组装身份不明（见7.2，优先级已下调）。
- 已经有一个独立的 `results/07.wild_rice_alignment/`（minimap2 asm10 预设，
  野生稻组装比对到IRGSP，用于paftools.js call变异，标注"进行中"）——这是
  另一条已经在跑的、跟competitive mapping平行的分析线，不要跟best-hit阶段
  混淆，但如果这些wild genome的minimap2比对已经产出了变异，可能可以复用来
  确认这些野生稻组装的物种/谱系身份。

**现在真正优先级最高的是0.5节的taxid修正**，其次是7.5（损伤窗口校准，
已有初步数据支撑"5bp可用"）。这条可以往后放。

### 8.2 时间序列选择信号扫描（selection scan）

背景（用户和同门的聊天记录）：想找"某个位点等位基因频率随时间**定向、快速**
变化"的信号，同门评价"很容易看""碰运气"。

**可行性评估**：
- 全项目16个古代样本，年代跨度约公元1160-1981年（近800年）——分到各个
  时间段可用个体数很少，是这类分析最大的硬约束。
- "碰运气/抽的多了总能碰上"这个说法点出了真实的方法论风险：不做多重检验
  校正/零分布模拟，genome-wide扫描出来的"候选位点"里假阳性会很多，纯粹的
  遗传漂变在小样本下也能产生看起来"剧烈"的频率波动。
- **正面信息**：项目之前已经搭建过 NGSNGS 模拟环境（git log:
  "NGSNGS simulation environment setup"，commit `bf107e84`），这个可以用来
  生成"纯漂变、无选择"下的零分布，给观测到的频率变化做显著性检验，而不是
  只看"变化快不快"这种没有统计标定的判断。
- **建议**：与其一上来做无先验的genome-wide扫描，优先聚焦在已经有生物学
  先验的候选基因（比如8.3提到的57个开花基因，本身就是驯化/生态型相关性状
  候选），信噪比会更好，也更容易向外解释结果。
- **这条线依赖第8.1**：genotype矩阵/等位基因频率的计算，前提是有干净、
  一致的比对/genotype pipeline，参考基因组整理不清楚会直接影响这里的
  可信度。

**还没做的事**：具体用哪些位点/genotype来源（是besthit过滤后的read直接
call变异，还是走已有的 `results/02.irgsp/` 主流程产出）、零分布模拟的具体
参数、时间分段方式，都没有讨论到细节，需要下一步单独设计。

### 8.3 57个开花基因：探针内部/两侧变异判定生态型（旱稻 vs 水稻）

背景：57个开花基因清单在 `db/gene/flower_gene.txt`（含 Level/Group/MSU_id/
RAPDB_id/Name），探针设计覆盖基因两侧（`db/gene/flower_gene.flank1kb.bed`，
±1kb）和基因内部。想法是找基因内部/附近的"典型插入"（结构变异）来判断
古代样本是旱稻还是水稻生态型。

**方法论上是站得住的**：这不是拍脑袋——项目自己已经有先例，git log里查过
DTH8/Ghd8经典大片段缺失（3024份现代品种里5.4%携带这个缺失，`commit
75f7e81`），这正是"用已知SV有无判定表型/生态型"这套方法在本项目里的先例。
现代水稻遗传学里也确实有一批性状对应已知的causal SV（比如落粒基因sh4、
半矮秆sd1、芒退化An-1等），不是凭空设计。

**最大的坑，项目自己已经踩过一次**：git log里同一批commit还有一条
"DTH8 breakpoint zero-coverage result"（commit `bf107e84`）——**哪怕是这么一个
研究得很透的明星基因，古DNA在这个位点的覆盖度都不够看**，这个已经提前预警了
"低覆盖度能不能可靠call出SV有无"这个问题在本项目里是真实存在的，不是理论
担忧。盲目对57个基因逐一去找"典型插入"，大概率在多数基因上撞到同样的墙。

**建议的下一步**：
1. **先做覆盖度QC**，不要直接冲进去分析。用现有 `results/02.irgsp/`（或
   `results/igv_package/bam_q30/`）BAM，对57个基因（含±1kb flank）逐样本
   算深度，看到底有几个基因、几个样本的覆盖度真的够格尝试call SV有无。这个
   半天工作量就能出结果，能立刻告诉你这条路现实不现实，比直接逐基因分析
   划算得多。
2. 确认这57个基因里**哪些在文献里真的有已知的旱稻/水稻区分性SV**（不是
   每个开花基因都有已知causal variant）——`flower_gene.txt` 的 Level/Group
   字段可能已经有部分线索，需要人工核实。
3. 覆盖度QC + 已知SV筛选两步都做完之后，能分析的基因数量可能远小于57，
   到时候再具体设计怎么从探针捕获的reads里判定SV有无（比如split-read/
   discordant-pair证据，还是简单的深度骤降，取决于具体是哪种SV）。

## 9. 给下一个接手者的具体待办（按优先级，2026-08-11更新后已重排）

0. 【2026-08-11当前正卡住，等用户反馈，比下面所有条目都优先】16个样本
   mapping已经全部跑完。开始实测5.1b节的besthit v2脚本时，服务器上跑
   `bash submit_oryza_besthit.sh check`打出来的还是v1的老格式输出
   （一行`ORYZA_TAXIDS=... DAMAGE_WINDOW=... TOP_N=...`），不是v2应该
   打的`Oryza scope: whole genus, ORYZA_GENUS_TAXID=4527`——怀疑服务器
   `scripts/`目录里的`submit_oryza_besthit.sh`没有被8号推送的v2版本真正
   覆盖（`curl`那步没生效或没执行）。已经请用户跑
   `grep -c "genus" submit_oryza_besthit.sh oryza_besthit_damage_filter.py`
   确认文件内容，并给了带`-o`强制指定文件名的重新下载命令。**接手时先
   要这个grep结果，如果还没跑，先让用户跑，不要跳过这一步直接假设脚本
   已经是v2去分析besthit结果**——如果grep结果是0，说明确实没覆盖成功，
   需要再排查是不是目录/权限问题；如果不是0，那说明只是之前那次check
   凑巧读到了缓存/旧终端里的输出，重新跑一次`check`应该就正常了。
1. 【0.6节，acc2taxid taxid修正】服务器dry-run已经跑过一次，
   证实互换方向完全正确（详见0.6节的真实输出数字），但发现`irgsp.acc2taxid`
   有一行表头会污染合并文件，脚本已修复（commit `2965348`）但**这个修复
   本身还没在服务器上重新验证过**。步骤：①重新下载最新脚本；②再跑一次
   dry-run，确认Step 2b显示`irgsp.acc2taxid 非数字taxid行数 = 1`且被正确
   识别；③确认无误后`--apply`，完整重建`all_wgs_asian_irgsp.acc2taxid`；
   ④**重新统计第6.1节的物种contig数表并重新评估"sativa覆盖度偏低"这个
   结论是否仍然成立**（这次不是小幅数字修正，"sativa 111条/rufipogon 717条"
   这两个数字本身可能是对调的，结论有可能反转，不能假设方向不变），更新
   这份文档。
2. 【次优先】7.5节：损伤窗口诊断——已有4个样本的初步数据支撑"5bp可用"，
   如果用户认为证据已经够用，可以不再深挖，直接按当前参数继续跑剩余样本；
   如果想更严谨验证，见7.5节的诊断脚本思路。
3. 【诊断，优先级降低，非阻塞】7.2（`db/3k/wild/`物种身份，与main分支
   `RESEARCH_ROADMAP.md`的P0第2条是同一个问题）、7.3（taxid rank校验）、
   7.4（MD tag覆盖率完整log、`seff`内存实测）——都可以做，但不再挡着
   best-hit往前走。
4. 只有在7.5的损伤窗口问题排查/调整完之后，才回头评估要不要给
   `oryza_besthit_damage_filter.py` 加一个margin参数（当前硬编码`<=`，没有
   margin概念）——在没排除损伤窗口这个更根本的原因之前，先调margin是在
   错的层面上打补丁。
5. 对剩余13个还在mapping队列里的样本，随mapping进度用 `submit_oryza_besthit.sh
   submit all` 陆续跑besthit（这条不依赖1-4，可以随时做）。
6. 8.1（参考基因组体系整理）——`asn720data`标签发现后优先级有所回升，
   但仍不是阻塞项。
7. 8.2（selection scan）和8.3（57基因SV判生态型）都还没真正开始。8.3的
   覆盖度QC是个独立、低成本、能快速出结论的子任务，可以和上面几条并行推进；
   8.2依赖更完整的genotype pipeline，且同样会受损伤/参考质量问题影响，
   优先级排在0.5/7.5之后。
8. 【已完成，见5.1b】5.1末尾记录的师兄脚本对比已经落地成v2重写（genus-wide
   Oryza范围+可选质量预筛），逐条取舍理由见5.1b节的表格。质量预筛的具体
   阈值(`--min-best-similarity`/`--max-best-raw-nm`)如果想启用，仍然需要
   在我们数据上摸索合适的值，目前默认关闭。
9. 【当前新增，优先级紧随0.5/0.6之后】5.1b节：v2版besthit主脚本
   （genus-wide Oryza范围）只在本地合成数据上验证过，**没有跑过真实BAM**。
   下载新脚本后先`smoke`一个已有样本，对比v2和v1的KEEP/REJECT分布差异
   是否符合预期（不应该暴增——smoke test数据显示"惜败"read的头号竞争者
   基本都是跟Oryza无关的物种，见6节），确认没问题后再考虑要不要重新跑
   已完成的4个样本（v1/v2的`summary.tsv`列数不同，不要混着`merge`）。
