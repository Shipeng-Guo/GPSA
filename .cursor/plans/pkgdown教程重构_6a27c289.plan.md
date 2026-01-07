---
name: Pkgdown教程重构
overview: 重写并拆分当前vignettes与README，使pkgdown站点首页=README（项目背景/特性/场景），Get Started只包含GitHub安装与2个Hallmark小例子，其余功能拆成多个独立英文教程文章；同时清除包内出现的中文内容。
todos:
  - id: rewrite-readme-home
    content: 重写 `README.md` 为全英文首页文案（背景/场景/特性/核心概念），并将安装方式收敛为 GitHub 安装
    status: completed
  - id: revise-get-started
    content: 更新 `vignettes/GPSA.Rmd`：仅 GitHub 安装 + 2个Hallmark示例（Estrogen suppressed、Interferon activated），其余高级示例移出
    status: completed
    dependencies:
      - rewrite-readme-home
  - id: split-case-study
    content: 拆分 `vignettes/case_study.Rmd`：迁移功能块到多篇新教程；保留一个全英文“Complete Case Study/Overview”页并链接各教程
    status: completed
    dependencies:
      - revise-get-started
  - id: add-tutorial-vignettes
    content: 新增并撰写 4 篇全英文教程 vignette（Up/Down含固定500与自适应；Viral mimicry事件；Batch+Meta；Vector search）
    status: completed
    dependencies:
      - split-case-study
  - id: update-pkgdown-nav
    content: 更新 `_pkgdown.yml`：重组 Tutorial 菜单与 articles 列表，使多教程以专业方式展示
    status: completed
    dependencies:
      - add-tutorial-vignettes
  - id: englishize-example-script
    content: 翻译并整理 `inst/examples/quick_start.R` 为全英文，与新文档口径一致
    status: completed
    dependencies:
      - update-pkgdown-nav
  - id: remove-all-chinese
    content: 最终核对并清除文档/教程相关文件中的中文字符，确保“全英文化”达标
    status: completed
    dependencies:
      - englishize-example-script
---

# GPSA 文档与教程重构计划

## 目标

- **首页专业化**：`README.md` 与 pkgdown 首页一致，聚焦项目背景、使用场景、核心概念与特性（不再堆叠教程细节）。
- **Get Started 聚焦入门**：仅保留 **GitHub 安装方式** + **2个Hallmark小例子**（雌激素通路“抑制”、干扰素信号“激活”）。
- **Tutorial 页面拆分**：把当前“巨型 case study”拆成多个独立教程文章，每个例子一篇，展示得更像专业文档站。
- **全英文化**：R 包与教程（README、vignettes、inst/examples）中不再出现中文。

## 将修改/新增的关键文件

- 更新首页内容：[`/home/shipeng/wsl_projects/GPSA_projects/GPSA/README.md`](/home/shipeng/wsl_projects/GPSA_projects/GPSA/README.md)
- 更新 Get Started：[`/home/shipeng/wsl_projects/GPSA_projects/GPSA/vignettes/GPSA.Rmd`](/home/shipeng/wsl_projects/GPSA_projects/GPSA/vignettes/GPSA.Rmd)
- 拆分/重写 case study：[`/home/shipeng/wsl_projects/GPSA_projects/GPSA/vignettes/case_study.Rmd`](/home/shipeng/wsl_projects/GPSA_projects/GPSA/vignettes/case_study.Rmd)
- 新增多篇教程 vignette（新文件）：
- [`/home/shipeng/wsl_projects/GPSA_projects/GPSA/vignettes/tutorial_updown.Rmd`](/home/shipeng/wsl_projects/GPSA_projects/GPSA/vignettes/tutorial_updown.Rmd)
- [`/home/shipeng/wsl_projects/GPSA_projects/GPSA/vignettes/tutorial_viral_mimicry.Rmd`](/home/shipeng/wsl_projects/GPSA_projects/GPSA/vignettes/tutorial_viral_mimicry.Rmd)
- [`/home/shipeng/wsl_projects/GPSA_projects/GPSA/vignettes/tutorial_batch.Rmd`](/home/shipeng/wsl_projects/GPSA_projects/GPSA/vignettes/tutorial_batch.Rmd)
- [`/home/shipeng/wsl_projects/GPSA_projects/GPSA/vignettes/tutorial_vector_search.Rmd`](/home/shipeng/wsl_projects/GPSA_projects/GPSA/vignettes/tutorial_vector_search.Rmd)
- 调整 pkgdown 导航与文章列表：[`/home/shipeng/wsl_projects/GPSA_projects/GPSA/_pkgdown.yml`](/home/shipeng/wsl_projects/GPSA_projects/GPSA/_pkgdown.yml)
- 翻译包内示例脚本：[`/home/shipeng/wsl_projects/GPSA_projects/GPSA/inst/examples/quick_start.R`](/home/shipeng/wsl_projects/GPSA_projects/GPSA/inst/examples/quick_start.R)

## 实施步骤

1. **重写 `README.md`（全英文）**

- 保留清晰结构：Motivation/Background、What GPSA does、When to use、Key features、How it works（rankZ/cMAP/GSEA-style 语义简述）、Links（Get Started/Tutorial/API Reference）。
- 安装方式仅保留 GitHub（例如 `remotes::install_github("Shipeng-Guo/GPSA")`），删除本地安装/多方法安装的中文段落。

2. **改造 `vignettes/GPSA.Rmd` 作为真正的 Get Started（全英文）**

- 安装仅 GitHub，并在文中明确项目地址 `https://github.com/Shipeng-Guo/GPSA/`。
- 入门示例改为 2 个 Hallmark gene set：
    - **Estrogen suppressed**：用 `msigdbr` 取 `HALLMARK_ESTROGEN_RESPONSE_EARLY`（默认）并以 `sign = -1` 查询。
    - **Interferon activated**：用 `HALLMARK_INTERFERON_ALPHA_RESPONSE`（默认）并以 `sign = +1` 查询。
- 其余高级用法移到 Tutorial 系列文章。

3. **把 `case_study.Rmd` 从“内容大杂烩”拆成“教程入口/完整串联”**

- 将当前 case study 中的各个功能块迁移到对应新 tutorial vignette。
- `case_study.Rmd` 本身改为全英文的“Complete Case Study”页：用一条连贯主线说明如何从 gene set → up/down → event → batch/meta → vector search，并链接到每篇独立教程。

4. **新增 4 篇独立教程 vignette（全英文）**

- **Up/Down tutorial**：固定 `n_each = 500` 的默认提取 + `gpsa_make_updown_query_adaptive()` 的能量阈值（eta/min_n/max_n）自适应方案。
- **Viral mimicry tutorial**：多通道事件定义（IFN up + cell cycle down）+ 门控筛选 + `gpsa_postprocess_modes()`（Phenocopy vs Discovery）。
- **Batch tutorial**：`gpsa_batch_search()` + query 相关性矩阵 + `gpsa_meta_stouffer()`。
- **Vector search tutorial**：`gpsa_search_vector()`（完整 geneList/logFC），包含稀疏化参数与 compile 信息解读。

5. **更新 `_pkgdown.yml` 以“专业方式展示”教程**

- Navbar 的 `Tutorial` 下拉菜单：
    - `Get Started`
    - 分隔线
    - Up/Down、Viral mimicry、Batch、Vector search
    - （可选）Complete case study
- `articles:` 内容按上述顺序组织，保证 `Tutorial` 页面自动生成分组清晰。

6. **翻译并整理 `inst/examples/quick_start.R`（全英文）**

- 去掉中文注释；保留脚本结构与功能演示，措辞与新教程保持一致。

7. **一致性检查**

- 确保 repo 内不再出现中文字符（至少 README、vignettes、inst/examples）。