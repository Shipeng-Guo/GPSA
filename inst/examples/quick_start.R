## ============================================================
## GPSA Quick Start Example
## ============================================================

## 设置镜像（中国用户）
options("repos" = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror = "https://mirrors.westlake.edu.cn/bioconductor")

library(GPSA)

## ============================================================
## 1) 设置HDF5数据库路径
## ============================================================
## 首次使用需下载HDF5数据库文件，然后设置路径

## 方法1: 设置R选项（推荐）
options(GPSA.h5 = "/path/to/gpsa_db.h5")

## 方法2: 设置环境变量
# Sys.setenv(GPSA_H5 = "/path/to/gpsa_db.h5")

## 方法3: 每次调用时指定
# db <- gpsa_load_db(h5file = "/path/to/gpsa_db.h5")

## ============================================================
## 2) 加载数据库
## ============================================================

## 内存模式（推荐）：预载Z矩阵，查询快（~0.15s/次）
## 加载约9秒，占用约1GB内存
db <- gpsa_load_db(preload_Z = TRUE, verbose = TRUE)

## 磁盘模式：不预载，查询慢（~7s/次），适合内存受限或单次查询
# db <- gpsa_load_db(preload_Z = FALSE, verbose = TRUE)

print(db)

## ============================================================
## 3) 单基因集查询
## ============================================================

## 示例：查询与某基因集相关的扰动
my_genes <- sample(db$genes, 300)  # 替换为实际基因集
query <- list(MySet = list(genes = my_genes, sign = +1))

out <- gpsa_search(db, query, return_components = FALSE)

message("\n=== 单基因集查询结果 ===")
print(head(out$res[, c("signature_id", "pbgene", "CellLineName", 
                        "Score", "tau", "FDR_rank")], 10))

## ============================================================
## 4) Up/Down双向查询
## ============================================================

## 从某个signature提取logFC，构建Up/Down查询
sig_id <- db$sigs[1]  # 替换为实际signature ID
gl <- gpsa_get_signature(db, sig_id, dataset = "logFC")

## 创建Up/Down查询（top/bottom 300基因）
query_ud <- gpsa_make_updown_query(gl, n_each = 300)

out_ud <- gpsa_search(db, query_ud, return_components = TRUE)

message("\n=== Up/Down查询 - Top 10 (最相似) ===")
print(head(out_ud$res[, c("signature_id", "pbgene", "CellLineName",
                          "Score", "tau", "FDR_rank")], 10))

message("\n=== Up/Down查询 - Bottom 10 (最相反) ===")
print(tail(out_ud$res[, c("signature_id", "pbgene", "CellLineName",
                          "Score", "tau", "FDR_rank")], 10))

## ============================================================
## 5) 双模式评分（Phenocopy vs Discovery）
## ============================================================

out_ud <- gpsa_postprocess_modes(out_ud, 
                                  attach_channel_cols = TRUE,
                                  p_mode = "norm",
                                  filter_positive = TRUE)

message("\n=== Phenocopy模式 (所有通道都满足) ===")
print(head(out_ud$res_similarity[, c("signature_id", "pbgene",
                                      "Score_bi", "Up_signed", "Down_signed")], 10))

message("\n=== Discovery模式 (至少一个通道强) ===")
print(head(out_ud$res_discovery[, c("signature_id", "pbgene",
                                     "Score_uni", "Driver", "Up_signed", "Down_signed")], 10))

## ============================================================
## 6) 批量查询
## ============================================================

## 准备多个查询
sig_ids <- db$sigs[1:3]
queries <- list()
for (sid in sig_ids) {
  gl <- gpsa_get_signature(db, sid, dataset = "logFC")
  queries[[sid]] <- gpsa_make_updown_query(gl, n_each = 300)
}

## 批量搜索
bat <- gpsa_batch_search(db, queries, compute_tau = TRUE)

message("\n=== 批量查询完成 ===")
print(bat$timing)

## 查询相似度矩阵
message("\n查询间相关性:")
print(round(cor(bat$Score, use = "pairwise.complete.obs"), 3))

## Meta分析
meta <- gpsa_meta_stouffer(bat$z_ref)
message("\nMeta分析 Top 10:")
print(head(meta, 10))

## ============================================================
## 7) GSEA可视化
## ============================================================

## 验证某个候选hit
candidate <- out_ud$res$signature_id[2]  # 第二名hit
gl_cand <- gpsa_get_signature(db, candidate, dataset = "logFC")

## 绘制Up基因富集
gpsa_plot_gsea(gl_cand, query_ud$Up$genes, 
               main = paste0(candidate, " vs Query Up genes"))

## 绘制Down基因富集
gpsa_plot_gsea(gl_cand, query_ud$Down$genes,
               main = paste0(candidate, " vs Query Down genes"))

message("\n=== Quick Start 完成 ===")

