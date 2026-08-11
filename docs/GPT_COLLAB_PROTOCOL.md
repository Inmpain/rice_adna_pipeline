# 网页版 GPT 协作协议

> 目的：让网页版 ChatGPT（通过其自带 GitHub 插件）能对本仓库做只读调研 +
> 提案/产出（论文方案、图表），同时不和 Claude Code / Codex CLI 正在做的
> 事情冲突，也不需要把本地文件系统暴露给云端。

最后更新: 2026-08-11

---

## 一、核心原则

1. **GitHub 是唯一总线**：GPT 走它自带的 GitHub 插件读写本仓库（已连接），
   不使用把本地目录暴露给云端的第三方 MCP 桥——那条路的风险（未审查
   MCP server + 和 Claude Code 同时改同一批文件没有锁）对这个项目不划算。
2. **`gpt/<topic>` 分支，永远 PR，不直接推 `main` 或任何 `codex/*` 分支**。
   `codex/*` 是 Claude Code / Codex CLI 的活跃工作区，GPT 的产出必须走
   独立分支 + PR，由你决定何时合并，避免和进行中的工作冲突（本仓库已经
   发生过一次两个并行 session 同晚差点冲突的情况，见
   `docs/ECOTYPE_PCA_PANEL.md` 末尾的并发编辑记录）。
3. **GPT 产出前必须先读现有文档**，不能凭空写方案——具体读哪些文档见下面
   每个场景的模板第1步，对应 `file_path.md` 第八节的快速定位表。
4. **不需要手动拼指令**：以后有具体任务时，让 Claude Code 现场把下面的
   模板填好具体分支名/文档路径/任务内容，你整段复制粘贴进网页 GPT 就行，
   不用自己回忆该读哪个文档、该建什么分支。

---

## 二、场景 A：新论文 → 方案提案

用途：拿到一篇新论文，想知道它对现有研究框架有什么影响、该怎么用。

**Claude Code 现场生成的成品 prompt 长这样**（示例已填好占位）：

```
你现在通过 GitHub 插件连接的仓库是 Inmpain/rice_adna_pipeline。在做任何
判断之前，请先做完这几步：

1. 读 main 分支的 file_path.md（整体路径地图+文档索引），再读
   <分支名，例如 codex/ecotype-pca-panel> 分支的
   <该分支HANDOFF文档，例如 docs/ECOTYPE_PCA_PANEL.md>
2. 我给你的论文是：<论文标题/DOI/上传的PDF>
3. 判断这篇论文对 docs/RESEARCH_ROADMAP.md 里的证据阶梯/五条工作线有什么
   增量，是否支持或推翻现有设计决策（尤其是 <具体决策，例如"PCA-A/PCA-B
   两个独立坐标系而非合并panel">），给出一个具体、可执行的下一步

完成后：
4. 把方案写成新文件 docs/proposals/<topic>_PROPOSAL.md
5. 新建分支 gpt/<topic>-proposal（不要碰 main 或任何 codex/* 分支），
   commit，PR 到 main
6. 在 docs/LITERATURE.md 追加一行，注明这篇论文已处理 + 方案文档路径
7. 把 PR 链接发给我
```

---

## 三、场景 B：Claude 管数据/代码，GPT 画图

用途：Claude Code 已经把数据产出提交进仓库，需要 GPT 基于这些数据出图
（示意图、统计图）。

**Claude Code 现场生成的成品 prompt**：

```
你现在通过 GitHub 插件连接的仓库是 Inmpain/rice_adna_pipeline。

1. 从 <分支名> 拉取 <数据文件路径>
2. 画图要求：<具体要求，例如"PCA-A散点图，按3K RGP官方亚群标签
   IND/AUS/ARO/TRJ/TEJ/ADM着色，古稻样本单独标出">
3. 输出图片（PNG 或 SVG）+ 一句话 caption，提交到
   docs/figures/<filename>，新建分支 gpt/figures-<topic>（不要碰 main
   或任何 codex/* 分支），commit，PR 到 main
4. 把 PR 链接发给我
```

Claude Code 在下一次交接时只引用图片路径，不改图片内容本身。

---

## 四、Claude Code 这边收到 GPT 的 PR 之后怎么处理

按 `github-repo-protocol` 的 Rule 1，下次扫全分支时会自然看到 `gpt/*`
分支的 PR。Claude 读完提案/确认图表后：

- **不自动合并** PR——合并永远是你的决定。
- 如果内容采纳，把关键结论折叠进对应层级的文档（`docs/RESEARCH_ROADMAP.md`
  / 分支级 HANDOFF 文档），并在 `file_path.md` 第九节补上新文档的索引行。
- 如果内容不采纳或需要修改，直接在 PR 里留评论说明原因，不代替你合并
  或关闭。
