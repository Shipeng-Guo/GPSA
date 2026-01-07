# Get Started with GPSA

## What is GPSA?

GPSA (Gene-set Perturbation Signature Analysis) provides **fast
similarity search** over large-scale perturbation signature databases
using cMAP/GSEA-style rankZ semantics.

Given a query gene set, GPSA finds perturbations (e.g., gene knockouts,
drug treatments) that produce **similar** or **opposite** transcriptomic
effects. This enables:

- Drug repurposing (find compounds that reverse a disease signature)
- Mechanism discovery (identify genes whose perturbation mimics your
  phenotype)
- Cross-study comparison (compare your results against public databases)

## Installation

Install GPSA from GitHub:

``` r
# Install remotes if needed
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Install GPSA from GitHub
remotes::install_github("Shipeng-Guo/GPSA")
```

### Optional Dependencies

``` r
# For GSEA visualization (recommended)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("fgsea")

# For MSigDB gene sets (used in examples)
install.packages("msigdbr")

# For advanced plotting
install.packages("ggplot2")
```

## First-time Setup: HDF5 Database

GPSA requires a pre-built HDF5 database containing perturbation
signatures. Download the database file and set its path:

``` r
library(GPSA)

# Set path once (persists for session)
options(GPSA.h5 = "/path/to/gpsa_db.h5")

# Alternative: environment variable
# Sys.setenv(GPSA_H5 = "/path/to/gpsa_db.h5")
```

## Loading the Database

GPSA offers two loading modes:

``` r
# Memory mode (recommended): ~0.15s per query, ~1GB RAM
db <- gpsa_load_db(preload_Z = TRUE, verbose = TRUE)

# Disk mode: ~7s per query, minimal RAM
# db <- gpsa_load_db(preload_Z = FALSE)
```

## Example 1: Estrogen Response (Suppressed)

Let’s find perturbations that **suppress** estrogen signaling. We’ll use
the Hallmark Estrogen Response Early gene set from MSigDB.

``` r
library(msigdbr)

# Get Hallmark gene sets for human
hallmark <- msigdbr(species = "Homo sapiens", category = "H")

# Helper function to extract genes
get_hallmark_genes <- function(gs_name) {
  unique(hallmark$gene_symbol[hallmark$gs_name == gs_name])
}

# Get estrogen response genes
er_genes <- get_hallmark_genes("HALLMARK_ESTROGEN_RESPONSE_EARLY")
message(sprintf("Estrogen Response Early: %d genes", length(er_genes)))
```

### Query for Estrogen Suppression

To find perturbations that **suppress** estrogen signaling, we use
`sign = -1`:

``` r
# Define query: negative sign means we want perturbations where these genes go DOWN
query_er_down <- list(
  ER_SUPPRESSED = list(
    genes = er_genes,
    sign = -1,      # Looking for suppression
    weight = 1
  )
)

# Run search
out_er <- gpsa_search(db, query_er_down, return_components = TRUE)

# View top hits (perturbations that suppress estrogen response)
message("=== Top 10 Estrogen Suppressors ===")
print(head(out_er$res[, c("signature_id", "pbgene", "CellLineName", 
                          "Score", "tau", "FDR_rank")], 10))
```

### Validate with GSEA Plot

``` r
# Get the logFC profile of the top hit
top_hit <- out_er$res$signature_id[1]
gl_top <- gpsa_get_signature(db, top_hit, dataset = "logFC")

# Plot GSEA running enrichment
gpsa_plot_gsea(gl_top, er_genes, 
               main = paste0(top_hit, " vs Estrogen Response Early"))
```

A negative enrichment score confirms that estrogen response genes are
indeed down-regulated in this perturbation.

## Example 2: Interferon Response (Activated)

Now let’s find perturbations that **activate** interferon signaling.

``` r
# Get interferon alpha response genes
ifn_genes <- get_hallmark_genes("HALLMARK_INTERFERON_ALPHA_RESPONSE")
message(sprintf("Interferon Alpha Response: %d genes", length(ifn_genes)))
```

### Query for Interferon Activation

To find perturbations that **activate** interferon signaling, we use
`sign = +1`:

``` r
# Define query: positive sign means we want perturbations where these genes go UP
query_ifn_up <- list(
  IFN_ACTIVATED = list(
    genes = ifn_genes,
    sign = +1,      # Looking for activation
    weight = 1
  )
)

# Run search
out_ifn <- gpsa_search(db, query_ifn_up, return_components = TRUE)

# View top hits (perturbations that activate interferon response)
message("=== Top 10 Interferon Activators ===")
print(head(out_ifn$res[, c("signature_id", "pbgene", "CellLineName",
                           "Score", "tau", "FDR_rank")], 10))
```

### Validate with GSEA Plot

``` r
# Get the logFC profile of the top hit
top_ifn <- out_ifn$res$signature_id[1]
gl_ifn <- gpsa_get_signature(db, top_ifn, dataset = "logFC")

# Plot GSEA running enrichment
gpsa_plot_gsea(gl_ifn, ifn_genes,
               main = paste0(top_ifn, " vs Interferon Alpha Response"))
```

A positive enrichment score confirms that interferon response genes are
up-regulated.

## Understanding the Results

The result table contains:

| Column         | Description                                           |
|----------------|-------------------------------------------------------|
| `signature_id` | Unique identifier for the perturbation                |
| `pbgene`       | Perturbed gene name                                   |
| `CellLineName` | Cell line where perturbation was performed            |
| `Score`        | Raw similarity score (higher = more similar to query) |
| `tau`          | Normalized score (percentile rank, -100 to +100)      |
| `p_rank`       | Rank-based p-value                                    |
| `FDR_rank`     | False discovery rate                                  |

**Interpreting tau scores:** - `tau > 90`: Strong positive connection
(similar to query) - `tau < -90`: Strong negative connection (opposite
of query) - `tau ~ 0`: No significant connection

## Next Steps

This tutorial covered single gene set queries. GPSA supports many more
advanced features:

- **[Up/Down Bidirectional
  Query](https://shipeng-guo.github.io/GPSA/articles/tutorial_updown.md)** -
  Find perturbations similar to a reference signature using both up- and
  down-regulated genes
- **[Multi-channel
  Events](https://shipeng-guo.github.io/GPSA/articles/tutorial_viral_mimicry.md)** -
  Define complex biological events (e.g., IFN up AND cell cycle down)
- **[Batch
  Queries](https://shipeng-guo.github.io/GPSA/articles/tutorial_batch.md)** -
  Run multiple queries at once with meta-analysis
- **[Full Vector
  Search](https://shipeng-guo.github.io/GPSA/articles/tutorial_vector_search.md)** -
  Query with complete expression profiles

## Session Info

``` r
sessionInfo()
```
