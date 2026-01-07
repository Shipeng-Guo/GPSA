# Get Started with GPSA

## What is GPSA?

GPSA (Gene-set Perturbation Signature Analysis) provides **fast
similarity search** over large-scale perturbation signature databases
using cMAP/GSEA-style rankZ semantics.

**Key Features:**

- Single gene-set queries
- Up/Down bidirectional queries (traditional cMAP-style)
- Multi-channel event definitions (e.g., viral mimicry)
- Dual-mode scoring: Phenocopy (strict) vs Discovery (one-sided)
- Batch queries with meta-analysis
- Full geneList vector search with adaptive sparsification
- HDF5-based storage for memory efficiency

## Installation

``` r
# Set mirrors (China users)
options("repos" = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror = "https://mirrors.westlake.edu.cn/bioconductor")

# Install dependencies
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("rhdf5")
install.packages("Matrix")

# Install GPSA from local source
install.packages("/path/to/GPSA", repos = NULL, type = "source")

# Or use devtools
# devtools::install_local("/path/to/GPSA")
```

## First-time Setup: HDF5 Database

GPSA requires a pre-built HDF5 database containing perturbation
signatures. Download the database file and set its path:

``` r
library(GPSA)

# Method 1: R option (recommended, persists for session)
options(GPSA.h5 = "/path/to/gpsa_db.h5")

# Method 2: Environment variable
# Sys.setenv(GPSA_H5 = "/path/to/gpsa_db.h5")

# Method 3: Specify h5file parameter in each call
# db <- gpsa_load_db(h5file = "/path/to/gpsa_db.h5")
```

## Loading the Database

GPSA offers two loading modes:

### Memory Mode (Recommended for interactive analysis)

``` r
# Preload Z matrix into RAM (~1GB, ~9s load time)
# Subsequent queries: ~0.15s each
db <- gpsa_load_db(preload_Z = TRUE, verbose = TRUE)
```

### Disk Mode (For limited memory)

``` r
# No preload, read from disk on each query
# Load time: ~0.04s, Query time: ~7s each
db <- gpsa_load_db(preload_Z = FALSE, verbose = TRUE)
```

## Quick Example: Single Gene Set Query

``` r
# Define a query gene set
query <- list(
  MyGenes = list(
    genes = c("TP53", "BRCA1", "EGFR", "MYC", "KRAS"),  # Your genes
    sign = +1,    # +1 for up-regulation, -1 for down-regulation
    weight = 1    # Optional weight
  )
)

# Search
out <- gpsa_search(db, query, return_components = FALSE)

# View results
head(out$res[, c("signature_id", "pbgene", "CellLineName", "Score", "tau", "FDR_rank")])
```

## Quick Example: Up/Down Query

Find perturbations similar to a reference signature:

``` r
# Get a reference signature's logFC
gl <- gpsa_get_signature(db, "D21455", dataset = "logFC")

# Create Up/Down query (top 300 + bottom 300 genes)
query_ud <- gpsa_make_updown_query(gl, n_each = 300)

# Search
out <- gpsa_search(db, query_ud, return_components = TRUE)

# Most similar perturbations
head(out$res)

# Most opposite perturbations
tail(out$res)
```

## Performance Summary

| Mode   | Query Type        | Time (avg) |
|--------|-------------------|------------|
| Memory | Single gene set   | 0.14s      |
| Memory | Up/Down 2-channel | 0.20s      |
| Memory | Multi 4-channel   | 0.25s      |
| Memory | Batch 5 queries   | 0.47s      |
| Disk   | Single gene set   | 7.2s       |

**Recommendation**: Use `preload_Z = TRUE` for interactive analysis.

## Next Steps

For complete examples including:

- **Estrogen pathway query** - Single gene set with GSEA validation
- **Viral mimicry event definition** - Multi-channel gated search
- **Dual-mode scoring** - Phenocopy vs Discovery
- **Batch queries and meta-analysis**
- **Full vector search with adaptive sparsification**

See the **[Complete Case
Study](https://shipeng-guo.github.io/GPSA/articles/case_study.md)**
vignette.

## API Reference

### Core Functions

| Function                                                                                     | Description              |
|----------------------------------------------------------------------------------------------|--------------------------|
| [`gpsa_load_db()`](https://shipeng-guo.github.io/GPSA/reference/gpsa_load_db.md)             | Load HDF5 database       |
| [`gpsa_search()`](https://shipeng-guo.github.io/GPSA/reference/gpsa_search.md)               | Gene-set query           |
| [`gpsa_search_vector()`](https://shipeng-guo.github.io/GPSA/reference/gpsa_search_vector.md) | Full geneList query      |
| [`gpsa_batch_search()`](https://shipeng-guo.github.io/GPSA/reference/gpsa_batch_search.md)   | Multiple queries at once |
| [`gpsa_build_db()`](https://shipeng-guo.github.io/GPSA/reference/gpsa_build_db.md)           | Build new database       |

### Query Helpers

| Function                                                                                                               | Description                            |
|------------------------------------------------------------------------------------------------------------------------|----------------------------------------|
| [`gpsa_make_updown_query()`](https://shipeng-guo.github.io/GPSA/reference/gpsa_make_updown_query.md)                   | Create Up/Down query                   |
| [`gpsa_make_updown_query_adaptive()`](https://shipeng-guo.github.io/GPSA/reference/gpsa_make_updown_query_adaptive.md) | Adaptive gene selection (energy-based) |
| [`gpsa_get_signature()`](https://shipeng-guo.github.io/GPSA/reference/gpsa_get_signature.md)                           | Extract signature vector               |

### Post-processing

| Function                                                                                                     | Description                             |
|--------------------------------------------------------------------------------------------------------------|-----------------------------------------|
| [`gpsa_postprocess_modes()`](https://shipeng-guo.github.io/GPSA/reference/gpsa_postprocess_modes.md)         | Dual-mode scoring (Phenocopy/Discovery) |
| [`gpsa_attach_channel_scores()`](https://shipeng-guo.github.io/GPSA/reference/gpsa_attach_channel_scores.md) | Add per-channel columns                 |
| [`gpsa_meta_stouffer()`](https://shipeng-guo.github.io/GPSA/reference/gpsa_meta_stouffer.md)                 | Meta-analysis across queries            |

### Visualization

| Function                                                                                   | Description                  |
|--------------------------------------------------------------------------------------------|------------------------------|
| [`gpsa_plot_gsea()`](https://shipeng-guo.github.io/GPSA/reference/gpsa_plot_gsea.md)       | GSEA running enrichment plot |
| [`gpsa_gsea_running()`](https://shipeng-guo.github.io/GPSA/reference/gpsa_gsea_running.md) | Compute running scores       |

## Session Info

``` r
sessionInfo()
```
