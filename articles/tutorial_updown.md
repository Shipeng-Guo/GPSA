# Tutorial: Up/Down Bidirectional Query

## Overview

When you have a reference perturbation signature (with both up-regulated
and down-regulated genes), you can search for other perturbations that
produce **similar** or **opposite** transcriptomic effects. This is the
classic Connectivity Map (cMAP) approach.

GPSA supports two methods for selecting query genes: 1. **Fixed
selection**: Take top N and bottom N genes by logFC 2. **Adaptive
selection**: Use energy threshold to determine optimal gene counts

## Setup

``` r
library(GPSA)

# Set database path
options(GPSA.h5 = "/path/to/gpsa_db.h5")

# Load database in memory mode
db <- gpsa_load_db(preload_Z = TRUE, verbose = TRUE)
```

## Method 1: Fixed Gene Selection

The classic approach uses a fixed number of genes from each tail of the
expression distribution.

### Extract a Reference Signature

``` r
# Choose a reference signature (e.g., an ESR1 perturbation)
sig_ref <- "D21455"  # Replace with actual signature ID

# Get the logFC vector
gl_ref <- gpsa_get_signature(db, sig_ref, dataset = "logFC")

# Preview the distribution
summary(gl_ref)
```

### Create Up/Down Query (Default: 500 each)

``` r
# Create query with top 500 up-regulated and bottom 500 down-regulated genes
query_ud <- gpsa_make_updown_query(gl_ref, n_each = 500)

message(sprintf("Up genes: %d", length(query_ud$Up$genes)))
message(sprintf("Down genes: %d", length(query_ud$Down$genes)))

# Preview the query structure
str(query_ud)
```

### Run the Search

``` r
out_ud <- gpsa_search(db, query_ud, return_components = TRUE)

# Check if reference ranks highly (should be near top)
self_rank <- gpsa_rank_in(out_ud$res, sig_ref)
message(sprintf("Reference signature self-rank: %d", self_rank))

# Top 10 most similar perturbations
message("=== Top 10 Most Similar ===")
print(head(out_ud$res[, c("signature_id", "pbgene", "CellLineName",
                          "Score", "tau", "FDR_rank")], 10))

# Bottom 10 most opposite perturbations
message("=== Top 10 Most Opposite ===")
print(tail(out_ud$res[, c("signature_id", "pbgene", "CellLineName",
                          "Score", "tau", "FDR_rank")], 10))
```

### Exclude Self to Find Similar Others

``` r
# Find the most similar perturbation that isn't the reference itself
top_other <- out_ud$res[out_ud$res$signature_id != sig_ref, ][1, ]
message("=== Most Similar Other Perturbation ===")
print(top_other[, c("signature_id", "pbgene", "CellLineName", "Score", "tau")])
```

## Method 2: Adaptive Gene Selection

Instead of arbitrary cutoffs (500 genes), use energy-based thresholding
to automatically determine how many genes capture a given fraction of
the signal.

### Adaptive Selection with Energy Threshold

``` r
# Select genes that capture 90% of the signal energy
q_adaptive <- gpsa_make_updown_query_adaptive(
  gl_ref,
  eta = 0.90,      # Capture 90% of energy
  min_n = 200,     # At least 200 genes per direction
  max_n = 2000     # At most 2000 genes per direction
)

# View selection information
message("=== Adaptive Selection Info ===")
print(q_adaptive$.info)

# The selected gene counts may differ from 500
message(sprintf("Adaptive Up genes: %d", length(q_adaptive$Up$genes)))
message(sprintf("Adaptive Down genes: %d", length(q_adaptive$Down$genes)))
```

### Run Search with Adaptive Query

``` r
# Note: Use only Up and Down, not the .info metadata
out_adaptive <- gpsa_search(db, q_adaptive[c("Up", "Down")], 
                            return_components = TRUE)

# Compare self-rank with fixed approach
self_rank_ad <- gpsa_rank_in(out_adaptive$res, sig_ref)
message(sprintf("Adaptive self-rank: %d (vs fixed: %d)", self_rank_ad, self_rank))
```

## Examining Channel-level Scores

GPSA computes separate enrichment scores for the Up and Down channels.
These can reveal whether a hit matches both directions or just one.

``` r
# Attach per-channel scores to results
res_channels <- gpsa_attach_channel_scores(
  out_ud$res, 
  out_ud$NES_raw, 
  out_ud$NES_signed
)

# View channel breakdown
message("=== Top 5 with Channel Scores ===")
print(head(res_channels[, c("signature_id", "pbgene", "Score",
                            "Up_raw", "Down_raw", 
                            "Up_signed", "Down_signed")], 5))
```

### Interpreting Channel Scores

| Score         | Meaning                                                     |
|---------------|-------------------------------------------------------------|
| `Up_raw`      | Raw NES for Up channel (unsigned)                           |
| `Down_raw`    | Raw NES for Down channel (unsigned)                         |
| `Up_signed`   | Signed NES for Up channel (+/- indicates direction match)   |
| `Down_signed` | Signed NES for Down channel (+/- indicates direction match) |

A hit with **both** `Up_signed > 0` AND `Down_signed > 0` shows
concordant regulation in both directions (true phenocopy).

## Dual-mode Scoring: Phenocopy vs Discovery

GPSA distinguishes two types of hits:

- **Phenocopy** (Score_bi): Both channels must show the expected
  direction (min of signed scores)
- **Discovery** (Score_uni): At least one channel shows strong signal
  (max of signed scores)

``` r
out_modes <- gpsa_postprocess_modes(
  out_ud,
  attach_channel_cols = TRUE,
  p_mode = "norm",
  filter_positive = TRUE
)

# Phenocopy hits: strict matches where both channels agree
message("=== Phenocopy Mode Top 10 (Score_bi) ===")
print(head(out_modes$res_similarity[, c("signature_id", "pbgene", 
                                         "Score_bi", "Score", "Score_uni",
                                         "Balance", "Driver")], 10))

# Discovery hits: strong in at least one channel
message("=== Discovery Mode Top 10 (Score_uni) ===")
print(head(out_modes$res_discovery[, c("signature_id", "pbgene",
                                        "Score_uni", "Score", "Score_bi",
                                        "Balance", "Driver")], 10))
```

### Interpreting Dual-mode Results

| Column      | Description                                           |
|-------------|-------------------------------------------------------|
| `Score_bi`  | Phenocopy score: min(Up_signed, Down_signed)          |
| `Score_uni` | Discovery score: max(Up_signed, Down_signed)          |
| `Balance`   | Ratio Score_bi / Score_uni (1 = perfectly balanced)   |
| `Driver`    | Which channel dominates (“Up”, “Down”, or “Balanced”) |

**Use cases:** - High `Score_bi` + High `Balance`: True phenocopy, both
directions match - High `Score_uni` + Low `Balance`: One-sided hit, only
one direction matches - `Driver` column tells you which channel is
driving the score

## GSEA Validation

Validate hits by examining the running enrichment:

``` r
# Get logFC for a top hit
candidate <- out_ud$res$signature_id[2]  # Second-ranked hit
gl_cand <- gpsa_get_signature(db, candidate, dataset = "logFC")

# Plot enrichment for Up genes
gpsa_plot_gsea(gl_cand, query_ud$Up$genes,
               main = paste0(candidate, " vs Query Up genes"))

# Plot enrichment for Down genes
gpsa_plot_gsea(gl_cand, query_ud$Down$genes,
               main = paste0(candidate, " vs Query Down genes"))
```

## Best Practices

1.  **Start with fixed selection** (`n_each = 500`) for initial
    exploration
2.  **Use adaptive selection** when you want interpretable, energy-based
    cutoffs
3.  **Always examine channel scores** to understand what’s driving the
    hit
4.  **Use dual-mode scoring** when you need to distinguish phenocopy
    from discovery
5.  **Validate top hits** with GSEA plots to confirm biological
    relevance

## Summary

| Approach     | Function                                          | When to Use                      |
|--------------|---------------------------------------------------|----------------------------------|
| Fixed 500    | `gpsa_make_updown_query(gl, n_each = 500)`        | Quick exploration, standard cMAP |
| Fixed custom | `gpsa_make_updown_query(gl, n_each = 300)`        | Smaller/larger gene sets         |
| Adaptive     | `gpsa_make_updown_query_adaptive(gl, eta = 0.90)` | Interpretable cutoffs            |

## Session Info

``` r
sessionInfo()
```
