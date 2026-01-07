# Tutorial: Batch Queries and Meta-analysis

## Overview

When you have multiple related signatures (e.g., biological replicates,
time points, or related conditions), running them individually is
inefficient. GPSA’s batch query feature:

- Runs multiple queries in a single matrix operation
- Computes query-to-query correlations
- Enables meta-analysis to find consensus hits

## Setup

``` r
library(GPSA)

# Set database path
options(GPSA.h5 = "/path/to/gpsa_db.h5")

# Load database
db <- gpsa_load_db(preload_Z = TRUE, verbose = TRUE)
```

## Preparing Multiple Queries

### Select Related Signatures

``` r
# Choose multiple related signatures
# Example: different perturbations of related genes
sig_ids <- c("D21455", "D20287", "D23923")  # Replace with actual IDs

# Filter to those present in database
sig_ids <- intersect(sig_ids, db$sigs)
message(sprintf("Using %d signatures for batch query", length(sig_ids)))
```

### Build Query List

``` r
# Create Up/Down query for each signature
queries <- list()
for (sid in sig_ids) {
  gl <- gpsa_get_signature(db, sid, dataset = "logFC")
  queries[[sid]] <- gpsa_make_updown_query(gl, n_each = 300)
}

# Preview query structure
names(queries)
```

## Running Batch Search

``` r
# Run all queries at once
bat <- gpsa_batch_search(db, queries, compute_tau = TRUE)

# View timing information
message("=== Batch Query Timing ===")
print(bat$timing)
```

### Batch Results Structure

``` r
# The result contains:
# - Score: matrix of scores (signatures x queries)
# - tau: matrix of tau scores
# - z_ref: matrix of z-scores for meta-analysis
# - timing: execution time

dim(bat$Score)  # n_signatures x n_queries
colnames(bat$Score)  # Query names
```

### View Results for Each Query

``` r
# Get results for a specific query
query_name <- names(queries)[1]
scores_q1 <- bat$Score[, query_name]

# Sort to get top hits for this query
top_q1 <- sort(scores_q1, decreasing = TRUE)[1:10]
message(sprintf("=== Top 10 for Query: %s ===", query_name))
print(top_q1)
```

## Query-to-Query Correlation

Understanding how similar your queries are helps interpret the results.

``` r
# Compute correlation between query results
cor_mat <- cor(bat$Score, use = "pairwise.complete.obs")

message("=== Query Similarity Matrix ===")
print(round(cor_mat, 3))
```

### Interpreting Correlations

| Correlation | Interpretation                            |
|-------------|-------------------------------------------|
| \> 0.7      | Highly similar queries (may be redundant) |
| 0.3 - 0.7   | Related but distinct queries              |
| \< 0.3      | Independent queries                       |
| Negative    | Opposite effects                          |

## Meta-analysis: Stouffer’s Method

When queries represent related conditions, combine their results to find
consensus hits.

``` r
# Stouffer's method combines z-scores across queries
meta_res <- gpsa_meta_stouffer(bat$z_ref)

# View meta-analysis results
message("=== Meta-analysis Top 20 ===")
print(head(meta_res[, c("signature_id", "z_meta", "p_meta", "FDR_meta")], 20))
```

### Understanding Meta-analysis Output

| Column     | Description                          |
|------------|--------------------------------------|
| `z_meta`   | Combined z-score (Stouffer’s method) |
| `p_meta`   | P-value from combined z-score        |
| `FDR_meta` | False discovery rate (BH correction) |

**Interpretation:** - High `z_meta`: Consistently high-scoring across
all queries - Low `FDR_meta`: Statistically significant consensus hit

## Finding Consensus Hits

``` r
# Define consensus: significant in meta-analysis
consensus <- meta_res[meta_res$FDR_meta < 0.05, ]
consensus <- consensus[order(consensus$z_meta, decreasing = TRUE), ]

message(sprintf("=== Consensus Hits (FDR < 0.05): %d ===", nrow(consensus)))
print(head(consensus, 20))
```

### Compare Individual vs Meta Results

``` r
# Check if a hit is consistent across queries
candidate <- consensus$signature_id[1]

# Get individual scores
individual_scores <- bat$Score[candidate, ]
message(sprintf("=== Individual Scores for %s ===", candidate))
print(individual_scores)

# Get individual tau scores
individual_tau <- bat$tau[candidate, ]
message(sprintf("=== Individual Tau Scores ==="))
print(individual_tau)
```

## Advanced: Weighted Meta-analysis

If queries have different reliability, apply weights:

``` r
# Example: weight by sample size or confidence
weights <- c(1.0, 0.8, 1.2)  # Adjust based on your criteria
names(weights) <- names(queries)

# Manual weighted Stouffer
z_weighted <- bat$z_ref %*% weights / sqrt(sum(weights^2))
```

## Batch Query with Different Query Types

You can mix different query types in a batch:

``` r
# Prepare mixed queries
mixed_queries <- list(
  # Up/Down query
  UD_query = gpsa_make_updown_query(
    gpsa_get_signature(db, sig_ids[1], dataset = "logFC"), 
    n_each = 300
  ),
  # Single gene set query (example)
  ER_query = list(
    ER = list(
      genes = c("ESR1", "PGR", "GREB1", "TFF1"),  # Example genes
      sign = +1
    )
  )
)

# Run mixed batch
bat_mixed <- gpsa_batch_search(db, mixed_queries)
```

## Visualizing Query Relationships

``` r
# Heatmap of query correlations (if ggplot2 available)
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  
  # Melt correlation matrix for plotting
  cor_df <- as.data.frame(as.table(cor_mat))
  names(cor_df) <- c("Query1", "Query2", "Correlation")
  
  ggplot(cor_df, aes(Query1, Query2, fill = Correlation)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = 0, limits = c(-1, 1)) +
    theme_minimal() +
    labs(title = "Query-to-Query Correlation")
}
```

## Performance Considerations

### Batch Size

``` r
# Batch queries share I/O overhead
# Running 5 queries in batch is much faster than 5 individual searches

# Timing comparison (illustrative)
# Individual: 5 x 0.20s = 1.0s
# Batch: ~0.47s (shared matrix operations)
```

### Memory Usage

``` r
# Each query adds minimal overhead in batch mode
# The database Z matrix (~1GB) is loaded once and reused
```

## Best Practices

1.  **Group related queries**: Batch queries that make biological sense
    together
2.  **Check correlations first**: Highly correlated queries may not add
    information
3.  **Use meta-analysis for consensus**: Stouffer’s method finds robust
    hits
4.  **Validate consensus hits**: They’re more likely to be biologically
    meaningful
5.  **Consider weights**: If query reliability varies, use weighted
    meta-analysis

## Summary

| Task              | Function                         |
|-------------------|----------------------------------|
| Batch search      | `gpsa_batch_search(db, queries)` |
| Query correlation | `cor(bat$Score)`                 |
| Meta-analysis     | `gpsa_meta_stouffer(bat$z_ref)`  |

## Session Info

``` r
sessionInfo()
```
