# GPSA: Gene-set Perturbation Signature Analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**GPSA** provides fast similarity search over large-scale perturbation signature databases using cMAP/GSEA-style rankZ semantics.

## Why GPSA?

Large-scale perturbation experiments (e.g., CRISPR knockouts, drug treatments) generate thousands of transcriptomic signatures. Given a query gene set or expression profile, researchers often need to find which perturbations produce **similar** or **opposite** effects. Traditional connectivity map (cMAP) approaches require iterating over thousands of signatures, which can be slow.

**GPSA** solves this by:

1. **Pre-computing rankZ-transformed signatures** stored in HDF5 format
2. **Enabling sub-second queries** via optimized matrix operations
3. **Supporting flexible query types** from simple gene sets to complex multi-channel events

## When to Use GPSA

- **Drug repurposing**: Find compounds that reverse a disease signature
- **Mechanism discovery**: Identify genes whose perturbation mimics your phenotype
- **Pathway analysis**: Query with pathway gene sets to find relevant perturbations
- **Cross-study comparison**: Compare your signature against public perturbation databases

## Key Features

| Feature | Description |
|---------|-------------|
| **Fast Search** | ~0.15s per query in memory mode, ~7s in disk mode |
| **Multiple Query Types** | Single gene set, Up/Down bidirectional, multi-channel events, batch queries |
| **HDF5 Storage** | Memory-efficient storage for 20k genes x 7k+ signatures |
| **Dual-mode Scoring** | Phenocopy (strict) and Discovery (one-sided) modes |
| **Full Vector Search** | Query with complete gene expression profiles |
| **GSEA Visualization** | Built-in running enrichment plots |

## How It Works

GPSA uses **rankZ semantics** similar to the Connectivity Map methodology:

1. **Database signatures** are rank-transformed and Z-scored (rankZ)
2. **Query gene sets** are compiled into sparse indicator vectors with direction signs
3. **Similarity scores** are computed via efficient matrix multiplication
4. **Statistical significance** is assessed using empirical null distributions (tau scores)

This approach captures the **relative ordering** of genes rather than absolute expression values, making it robust across platforms and batches.

## Installation

Install GPSA directly from GitHub:

```r
# Install remotes if needed
if (!requireNamespace("remotes", quietly = TRUE)) {

install.packages("remotes")
}

# Install GPSA
remotes::install_github("Shipeng-Guo/GPSA")
```

### Dependencies

GPSA requires the following packages (installed automatically):

- **rhdf5** (Bioconductor) - HDF5 file handling
- **Matrix** - Sparse matrix operations

Optional packages for enhanced functionality:

```r
# For GSEA plots (recommended)
BiocManager::install("fgsea")

# For MSigDB gene sets
install.packages("msigdbr")

# For advanced plotting
install.packages("ggplot2")
```

## Quick Start

```r
library(GPSA)

# Set path to HDF5 database (download separately)
options(GPSA.h5 = "/path/to/gpsa_db.h5")

# Load database (memory mode for fast queries)
db <- gpsa_load_db(preload_Z = TRUE)

# Query with a gene set
query <- list(
  MyGenes = list(genes = c("TP53", "BRCA1", "EGFR"), sign = +1)
)
out <- gpsa_search(db, query)

# View top hits
head(out$res)
```

## Documentation

- **[Get Started](https://shipeng-guo.github.io/GPSA/articles/GPSA.html)** - Installation and basic examples
- **[Up/Down Query Tutorial](https://shipeng-guo.github.io/GPSA/articles/tutorial_updown.html)** - Bidirectional similarity search
- **[Multi-channel Events](https://shipeng-guo.github.io/GPSA/articles/tutorial_viral_mimicry.html)** - Complex event definitions
- **[Batch Queries](https://shipeng-guo.github.io/GPSA/articles/tutorial_batch.html)** - Multiple queries and meta-analysis
- **[Vector Search](https://shipeng-guo.github.io/GPSA/articles/tutorial_vector_search.html)** - Full geneList queries
- **[API Reference](https://shipeng-guo.github.io/GPSA/reference/)** - Complete function documentation

## Performance

Database: 19,613 genes x 7,103 signatures

| Mode | Query Type | Time (avg) |
|------|------------|------------|
| Memory | Single gene set | 0.14s |
| Memory | Up/Down 2-channel | 0.20s |
| Memory | Multi 4-channel | 0.25s |
| Memory | Batch 5 queries | 0.47s |
| Disk | Single gene set | 7.2s |

**Recommendation**: Use `preload_Z = TRUE` for interactive analysis.

## Citation

If you use GPSA in your research, please cite:

> Shao S. (2025). GPSA: Gene-set Perturbation Signature Analysis. R package version 0.1.0. https://github.com/Shipeng-Guo/GPSA

## Building the Documentation Site

The pkgdown documentation site includes tutorials with **real executed results** from the HDF5 database. Since the database file is too large for GitHub, the site must be built locally.

### Prerequisites

1. Have the HDF5 database file available locally (e.g., `output/gpsa_rankZ_logFC_optimized.h5`)
2. Set the database path before building:

```r
# Option 1: Set R option
options(GPSA.h5 = "/path/to/gpsa_db.h5")

# Option 2: Set environment variable
Sys.setenv(GPSA_H5 = "/path/to/gpsa_db.h5")
```

### Build and Deploy

```r
# Install pkgdown if needed
install.packages("pkgdown")

# Build site locally (vignettes will execute with real data)
pkgdown::build_site()

# Deploy to gh-pages branch
pkgdown::deploy_to_branch()
```

The `deploy_to_branch()` function will:
1. Build the site with executed vignettes (tables, figures)
2. Commit the generated `docs/` to the `gh-pages` branch
3. Push to GitHub (requires write access)

**Note**: The GitHub Actions workflow is configured for manual trigger only to prevent overwriting the site with non-executed vignettes.

## License

MIT License

## Links

- [GitHub Repository](https://github.com/Shipeng-Guo/GPSA)
- [Documentation Site](https://shipeng-guo.github.io/GPSA/)
- [Report Issues](https://github.com/Shipeng-Guo/GPSA/issues)
