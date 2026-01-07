# GPSA: Gene-set Perturbation Signature Analysis

GPSA provides fast similarity search over large-scale perturbation signature databases using cMAP/GSEA-style rankZ semantics.

## Features

- **Fast Search**: ~0.15s per query in memory mode, ~7s in disk mode
- **Multiple Query Types**: Single gene set, Up/Down bidirectional, multi-channel events, batch queries
- **HDF5 Storage**: Memory-efficient storage for 20k genes × 7k+ signatures
- **Dual-mode Scoring**: Phenocopy (strict) and Discovery (one-sided) modes
- **Full Vector Search**: Query with complete gene expression profiles

## Installation

### 方法1: 从本地源码安装（推荐）

```r
# 设置镜像加速（中国用户）
options("repos" = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror = "https://mirrors.westlake.edu.cn/bioconductor")

# 安装依赖
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("rhdf5")
install.packages("Matrix")

# 从本地源码安装（将路径替换为实际GPSA目录）
install.packages("/path/to/GPSA", repos = NULL, type = "source")
```

### 方法2: 使用devtools安装

```r
# 安装devtools（如果没有）
install.packages("devtools")

# 从本地目录安装
devtools::install_local("/path/to/GPSA")

# 或直接加载（开发模式，无需安装）
devtools::load_all("/path/to/GPSA")
```

### 方法3: 终端命令行安装

```bash
cd /path/to/GPSA  # 进入包目录
R CMD build .     # 构建tar.gz包
R CMD INSTALL GPSA_*.tar.gz  # 安装
```

### 可选依赖

```r
# GSEA绘图（推荐）
BiocManager::install("fgsea")

# MSigDB基因集
install.packages("msigdbr")

# 高级绘图
install.packages("ggplot2")
```

## Quick Start

### 1. First-time Setup: Download HDF5 Database

Download the pre-built HDF5 database file and set the path:

```r
library(GPSA)

# Set path once (persists for session)
options(GPSA.h5 = "/path/to/gpsa_db.h5")

# Or set via environment variable
Sys.setenv(GPSA_H5 = "/path/to/gpsa_db.h5")
```

### 2. Load Database

```r
# Disk mode (slower queries, less memory)
db <- gpsa_load_db(preload_Z = FALSE)

# Memory mode (fast queries, ~1GB RAM)
db <- gpsa_load_db(preload_Z = TRUE)
```

### 3. Basic Query Examples

#### Single Gene Set Query

```r
# Query with a gene set
query <- list(MyGenes = list(genes = c("TP53", "BRCA1", "EGFR", ...), sign = +1))
out <- gpsa_search(db, query)
head(out$res)
```

#### Up/Down Bidirectional Query

```r
# Find perturbations similar to a signature
gl <- gpsa_get_signature(db, "D21455", dataset = "logFC")
query <- gpsa_make_updown_query(gl, n_each = 300)
out <- gpsa_search(db, query, return_components = TRUE)

# Most similar perturbations
head(out$res)

# Most opposite perturbations
tail(out$res)
```

#### Multi-channel Event Query

```r
# Define a complex event (e.g., viral mimicry: IFN up + cell cycle down)
query <- list(
  IFN_UP = list(genes = ifn_genes, sign = +1),
  CYCLE_DOWN = list(genes = cycle_genes, sign = -1)
)
out <- gpsa_search(db, query, return_components = TRUE)
out <- gpsa_postprocess_modes(out, p_mode = "norm")

# Phenocopy hits (all channels satisfied)
head(out$res_similarity)

# Discovery hits (at least one channel strong)
head(out$res_discovery)
```

#### Batch Query

```r
# Multiple queries at once
queries <- list(
  Q1 = gpsa_make_updown_query(gl1),
  Q2 = gpsa_make_updown_query(gl2),
  Q3 = gpsa_make_updown_query(gl3)
)
bat <- gpsa_batch_search(db, queries)
```

## Performance Benchmarks

Database: 19,613 genes × 7,103 signatures

| Mode | Query Type | Time (avg) |
|------|------------|------------|
| Memory | Single gene set | 0.14s |
| Memory | Up/Down 2-channel | 0.20s |
| Memory | Multi 4-channel | 0.25s |
| Memory | Batch 5 queries | 0.47s |
| Disk | Single gene set | 7.2s |
| Disk | Up/Down 2-channel | 8.5s |

**Recommendation**: Use `preload_Z = TRUE` for interactive analysis.

## Building Your Own Database

```r
# From a diffTable RDS (genes × signatures logFC matrix)
gpsa_build_db(
  rds_file = "diffTable.rds",
  h5file = "gpsa_db.h5",
  meta_file = "metadata.rds",  # optional
  compress_level = 1,
  store_logFC = TRUE,
  verbose = TRUE
)
```

## API Reference

### Core Functions

- `gpsa_load_db()` - Load database
- `gpsa_search()` - Gene-set query
- `gpsa_search_vector()` - Full geneList query
- `gpsa_batch_search()` - Multiple queries
- `gpsa_build_db()` - Build database

### Query Helpers

- `gpsa_make_updown_query()` - Create Up/Down query
- `gpsa_make_updown_query_adaptive()` - Adaptive gene selection
- `gpsa_get_signature()` - Extract signature vector

### Post-processing

- `gpsa_postprocess_modes()` - Dual-mode scoring
- `gpsa_attach_channel_scores()` - Add channel columns
- `gpsa_meta_stouffer()` - Meta-analysis

### Visualization

- `gpsa_plot_gsea()` - GSEA running plot
- `gpsa_gsea_running()` - Compute running scores

## License

MIT License

