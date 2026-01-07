# GPSA Case Study: 完整功能演示

## 本地安装GPSA包

``` r
## ============================================================
## 方法1: 从源码目录安装（推荐）
## ============================================================
# 设置镜像加速（中国用户）
options("repos" = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror = "https://mirrors.westlake.edu.cn/bioconductor")

# 安装依赖
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("rhdf5")
install.packages("Matrix")

# 从本地源码安装
install.packages("/path/to/GPSA", repos = NULL, type = "source")

## ============================================================
## 方法2: 使用devtools安装
## ============================================================
# install.packages("devtools")
devtools::install_local("/path/to/GPSA")

## ============================================================
## 方法3: 直接加载（开发模式，不安装）
## ============================================================
devtools::load_all("/path/to/GPSA")
```

## GPSA功能清单

GPSA支持以下核心功能：

1.  **快速建库** - 从diffTable RDS构建HDF5索引，支持float32/压缩/分块
2.  **快速加载** - 可选预载Z矩阵到内存（查询提速35-50倍）
3.  **单基因集查询** - 输入gene set，返回全库signature的Score/tau/p/FDR
4.  **Up/Down双向查询** - 传统cMAP风格，找最像/最相反的扰动
5.  **自适应基因选择** - 用能量阈值η代替固定500，更可解释
6.  **完整geneList向量检索** - rankZ语义，支持稀疏化
7.  **多通道事件定义** - 任意组合正/负通道，门控筛选
8.  **批量查询** - 多query一次性矩阵乘法
9.  **Meta分析** - Stouffer方法合并多query
10. **GSEA可视化** - 快速提取signature并绑制running curve

## 初始化设置

``` r
library(GPSA)

## 设置镜像（中国用户加速）
options("repos" = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror = "https://mirrors.westlake.edu.cn/bioconductor")

## 设置HDF5数据库路径（三种方式任选其一）
## 方式1: R选项（推荐，session内有效）
options(GPSA.h5 = "/path/to/gpsa_db.h5")

## 方式2: 环境变量
# Sys.setenv(GPSA_H5 = "/path/to/gpsa_db.h5")

## 方式3: 每次调用时指定h5file参数
```

## 加载数据库

``` r
## 内存模式（推荐）：预载Z矩阵，查询极快（~0.15s/次）
## 加载约9秒，占用约1GB内存
db <- gpsa_load_db(preload_Z = TRUE, verbose = TRUE)

## 磁盘模式：不预载，查询较慢（~7s/次），适合内存受限
# db <- gpsa_load_db(preload_Z = FALSE, verbose = TRUE)

print(db)
```

## Case Study 1: 雌激素信号通路查询

使用MSigDB Hallmark基因集查询雌激素响应相关的扰动。

``` r
## 加载Hallmark基因集
library(msigdbr)
hall <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_genes <- function(gs_name) unique(hall$gene_symbol[hall$gs_name == gs_name])

## 获取雌激素响应基因集
er_early <- hallmark_genes("HALLMARK_ESTROGEN_RESPONSE_EARLY")
er_late <- hallmark_genes("HALLMARK_ESTROGEN_RESPONSE_LATE")

message(sprintf("ER_EARLY基因数: %d", length(er_early)))
message(sprintf("ER_LATE基因数: %d", length(er_late)))
```

### 单基因集查询

``` r
## 查询ER_EARLY
q_er <- list(ER_EARLY = list(genes = er_early, sign = +1, weight = 1))

out_er <- gpsa_search(db, q_er, return_components = TRUE)

## 查看Top 10结果
message("=== ER_EARLY信号通路 Top 10 ===")
print(head(out_er$res[, c("signature_id", "pbgene", "CellLineName",
                          "Score", "tau", "p_rank", "FDR_rank")], 10))

## 查看Bottom 10（负向富集）
message("=== ER_EARLY信号通路 Bottom 10 ===")
print(tail(out_er$res[, c("signature_id", "pbgene", "CellLineName",
                          "Score", "tau", "p_rank", "FDR_rank")], 10))
```

### GSEA可视化验证

``` r
## 对Top hit进行GSEA验证
top_hit <- out_er$res$signature_id[1]
gl_top <- gpsa_get_signature(db, top_hit, dataset = "logFC")

## 绑制GSEA running curve
gpsa_plot_gsea(gl_top, er_early, 
               main = paste0(top_hit, " vs HALLMARK_ESTROGEN_RESPONSE_EARLY"))
```

## Case Study 2: Up/Down双向相似性查询

从某个扰动signature提取top/bottom基因，寻找产生相似/相反效应的扰动。

``` r
## 选择一个参考signature（如ESR1扰动）
sig_ref <- "D21455"  # 替换为实际signature ID

## 提取logFC
gl_ref <- gpsa_get_signature(db, sig_ref, dataset = "logFC")

## 创建Up/Down查询（top 500 + bottom 500）
query_ud <- gpsa_make_updown_query(gl_ref, n_each = 500)

message(sprintf("Up基因数: %d", length(query_ud$Up$genes)))
message(sprintf("Down基因数: %d", length(query_ud$Down$genes)))
```

### 执行查询

``` r
out_ud <- gpsa_search(db, query_ud, return_components = TRUE)

## 查看自身排名（应该很高）
self_rank <- gpsa_rank_in(out_ud$res, sig_ref)
message(sprintf("参考signature自身排名: %d", self_rank))

## Top 10最相似扰动
message("=== Top 10 最相似扰动 ===")
print(head(out_ud$res[, c("signature_id", "pbgene", "CellLineName",
                          "Score", "tau", "FDR_rank")], 10))

## Bottom 10最相反扰动
message("=== Bottom 10 最相反扰动 ===")
print(tail(out_ud$res[, c("signature_id", "pbgene", "CellLineName",
                          "Score", "tau", "FDR_rank")], 10))

## 排除自身后的最相似扰动
top_other <- out_ud$res[out_ud$res$signature_id != sig_ref, ][1, ]
message("=== 最相似的OTHER扰动 ===")
print(top_other[, c("signature_id", "pbgene", "CellLineName", "Score", "tau")])
```

### 查看Up/Down通道得分

``` r
## 附加通道分数
res_with_channels <- gpsa_attach_channel_scores(out_ud$res, 
                                                 out_ud$NES_raw, 
                                                 out_ud$NES_signed)

message("=== Top 5 with Up/Down channel scores ===")
print(head(res_with_channels[, c("signature_id", "pbgene", "Score",
                                  "Up_raw", "Down_raw", 
                                  "Up_signed", "Down_signed")], 5))
```

## Case Study 3: Viral Mimicry多通道事件定义

Viral mimicry是一种细胞状态，特征为： - 干扰素响应 **上调** (IFN-α,
IFN-γ) - 细胞周期 **下调** (E2F targets, G2M checkpoint)

``` r
## 获取相关基因集
ifna <- hallmark_genes("HALLMARK_INTERFERON_ALPHA_RESPONSE")
ifng <- hallmark_genes("HALLMARK_INTERFERON_GAMMA_RESPONSE")
e2f  <- hallmark_genes("HALLMARK_E2F_TARGETS")
g2m  <- hallmark_genes("HALLMARK_G2M_CHECKPOINT")

message(sprintf("IFNA: %d genes, IFNG: %d genes", length(ifna), length(ifng)))
message(sprintf("E2F: %d genes, G2M: %d genes", length(e2f), length(g2m)))

## 定义Viral mimicry查询（4通道）
q_vm <- list(
  IFNA = list(genes = ifna, sign = +1),  # IFN-α上调
  IFNG = list(genes = ifng, sign = +1),  # IFN-γ上调
  E2F  = list(genes = e2f,  sign = -1),  # E2F下调
  G2M  = list(genes = g2m,  sign = -1)   # G2M下调
)
```

### 执行多通道查询

``` r
out_vm <- gpsa_search(db, q_vm, return_components = TRUE)

message("=== Viral mimicry Top 10 (总分排序) ===")
print(head(out_vm$res[, c("signature_id", "pbgene", "CellLineName",
                          "Score", "tau", "FDR_rank")], 10))
```

### 门控筛选（所有通道都满足方向）

``` r
## 附加通道分数并进行门控
res_vm <- gpsa_attach_channel_scores(out_vm$res, out_vm$NES_raw, out_vm$NES_signed)

## 门控条件：所有NES_signed > 0
gated <- res_vm[res_vm$IFNA_signed > 0 & 
                res_vm$IFNG_signed > 0 & 
                res_vm$E2F_signed > 0 & 
                res_vm$G2M_signed > 0, ]
gated <- gated[order(gated$Score, decreasing = TRUE), ]

message(sprintf("=== 门控后命中数: %d ===", nrow(gated)))
message("=== Viral mimicry门控Top 20 ===")
print(head(gated[, c("signature_id", "pbgene", "CellLineName", "Score",
                     "IFNA_signed", "IFNG_signed", "E2F_signed", "G2M_signed")], 20))
```

## Case Study 4: 双模式评分（Phenocopy vs Discovery）

双模式评分区分两种场景： - **Phenocopy (Score_bi)**:
min(所有通道NES_signed)，要求所有通道都满足 - **Discovery (Score_uni)**:
max(所有通道NES_signed)，只要某通道特别强

``` r
## 继续使用Up/Down查询结果
out_modes <- gpsa_postprocess_modes(out_ud, 
                                     attach_channel_cols = TRUE,
                                     p_mode = "norm",
                                     filter_positive = TRUE)

## Phenocopy模式：严格表型复制
message("=== Phenocopy模式 Top 10 (Score_bi排序) ===")
print(head(out_modes$res_similarity[, c("signature_id", "pbgene", 
                                         "Score_bi", "Score", "Score_uni", 
                                         "Balance", "Driver")], 10))

## Discovery模式：发现有一个通道特别强的命中
message("=== Discovery模式 Top 10 (Score_uni排序) ===")
print(head(out_modes$res_discovery[, c("signature_id", "pbgene",
                                        "Score_uni", "Score", "Score_bi",
                                        "Balance", "Driver")], 10))
```

### 案例：为什么某个hit是Discovery而非Phenocopy

``` r
## 选择一个感兴趣的候选（如OVOL2）
candidate <- gpsa_pick_signature(db, pbgene = "OVOL2", cellline = "MCF7")

if (!is.null(candidate)) {
  ## 查看该候选在不同排序下的排名
  r_score <- gpsa_rank_in(out_modes$res, candidate)
  r_bi <- gpsa_rank_in(out_modes$res_similarity, candidate)
  r_uni <- gpsa_rank_in(out_modes$res_discovery, candidate)
  
  message(sprintf("候选 %s 排名:", candidate))
  message(sprintf("  总分(Score)排名: %s", ifelse(is.na(r_score), "NA", r_score)))
  message(sprintf("  Phenocopy(Score_bi)排名: %s", ifelse(is.na(r_bi), "NA", r_bi)))
  message(sprintf("  Discovery(Score_uni)排名: %s", ifelse(is.na(r_uni), "NA", r_uni)))
  
  ## 查看详细通道分数
  detail <- out_modes$res[out_modes$res$signature_id == candidate, ]
  message("=== 候选详细信息 ===")
  print(detail[, c("signature_id", "pbgene", "Score", 
                   "Up_signed", "Down_signed",
                   "Score_bi", "Score_uni", "Driver")])
}
```

## Case Study 5: 完整geneList向量检索

使用完整的gene expression profile进行检索，而非仅top/bottom基因。

``` r
## 使用完整logFC向量
gl_full <- gpsa_get_signature(db, sig_ref, dataset = "logFC")

## 向量检索（rankZ变换 + 自适应稀疏化）
out_vec <- gpsa_search_vector(db, gl_full,
                               eta = 0.90,       # 捕获90%能量
                               min_k = 300,      # 最少300基因
                               max_k = 5000,     # 最多5000基因
                               transform = "rankZ",
                               normalize = TRUE)

## 查看编译信息
message("=== 向量检索编译信息 ===")
print(out_vec$compile)

## 检查自身排名
self_rank_vec <- gpsa_rank_in(out_vec$res, sig_ref)
message(sprintf("向量检索中自身排名: %d", self_rank_vec))

## Top 10
message("=== 向量检索 Top 10 ===")
print(head(out_vec$res[, c("signature_id", "pbgene", "CellLineName",
                           "Score", "tau", "FDR_rank")], 10))
```

### 自适应Up/Down选择（用能量阈值代替固定500）

``` r
## 自适应选择（基于能量捕获）
q_adaptive <- gpsa_make_updown_query_adaptive(gl_full, 
                                               eta = 0.90,
                                               min_n = 200,
                                               max_n = 2000)

message("=== 自适应Up/Down信息 ===")
print(q_adaptive$.info)

out_adaptive <- gpsa_search(db, q_adaptive[c("Up", "Down")], 
                            return_components = FALSE)

self_rank_ad <- gpsa_rank_in(out_adaptive$res, sig_ref)
message(sprintf("自适应Up/Down中自身排名: %d", self_rank_ad))
```

## Case Study 6: 批量查询与Meta分析

同时运行多个query，并进行相关性分析和meta合并。

``` r
## 选择多个signatures
sig_batch <- c("D21455", "D20287", "D23923")  # 替换为实际IDs
sig_batch <- intersect(sig_batch, db$sigs)

## 为每个signature构建Up/Down查询
queries <- list()
for (sid in sig_batch) {
  gl <- gpsa_get_signature(db, sid, dataset = "logFC")
  queries[[sid]] <- gpsa_make_updown_query(gl, n_each = 300)
}

## 批量搜索
bat <- gpsa_batch_search(db, queries, compute_tau = TRUE)

message("=== 批量查询耗时 ===")
print(bat$timing)
```

### Query间相似性分析

``` r
## 计算query之间的相关性
cor_mat <- cor(bat$Score, use = "pairwise.complete.obs")
message("=== Query相似性矩阵 ===")
print(round(cor_mat, 3))
```

### Stouffer Meta分析

``` r
## 合并多个query的结果
meta_res <- gpsa_meta_stouffer(bat$z_ref)

message("=== Meta分析 Top 20 ===")
print(head(meta_res[, c("signature_id", "z_meta", "p_meta", "FDR_meta")], 20))
```

## 性能基准

在19,613 genes × 7,103 signatures规模下的性能：

| 模式     | 查询类型      | 平均耗时 | 标准差 |
|----------|---------------|----------|--------|
| 内存模式 | 单基因集      | 0.14s    | 0.004s |
| 内存模式 | Up/Down双通道 | 0.20s    | 0.04s  |
| 内存模式 | 4通道事件     | 0.25s    | 0.02s  |
| 内存模式 | 批量5查询     | 0.47s    | 0.07s  |
| 磁盘模式 | 单基因集      | 7.2s     | 0.1s   |
| 磁盘模式 | Up/Down双通道 | 8.5s     | 2.4s   |

**建议**：交互分析时使用 `preload_Z = TRUE`。

## 最佳实践

1.  **设置镜像加速**（中国用户）：

``` r
options("repos" = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror = "https://mirrors.westlake.edu.cn/bioconductor")
```

2.  **使用内存模式**进行交互分析：

``` r
db <- gpsa_load_db(preload_Z = TRUE)
```

3.  **批量查询**共享I/O开销：

``` r
bat <- gpsa_batch_search(db, list(Q1 = q1, Q2 = q2, Q3 = q3))
```

4.  **自适应基因选择**避免任意cutoff：

``` r
query <- gpsa_make_updown_query_adaptive(gl, eta = 0.90)
```

5.  **双模式评分**区分表型复制和功能发现：

``` r
out <- gpsa_postprocess_modes(out, p_mode = "norm")
head(out$res_similarity)  # Phenocopy hits
head(out$res_discovery)   # Discovery hits
```

## Session Info

``` r
sessionInfo()
```
