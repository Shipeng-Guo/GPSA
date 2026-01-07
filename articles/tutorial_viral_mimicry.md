# Tutorial: Multi-channel Event Definition (Viral Mimicry)

## Overview

Some biological phenotypes cannot be captured by a single gene set
direction. For example, **viral mimicry** is a cellular state
characterized by:

- **Interferon response genes UP** (immune activation)
- **Cell cycle genes DOWN** (proliferation arrest)

GPSA allows you to define such **multi-channel events** by combining
multiple gene sets with different expected directions.

## Setup

``` r
library(GPSA)
library(msigdbr)

# Set database path
options(GPSA.h5 = "/path/to/gpsa_db.h5")

# Load database
db <- gpsa_load_db(preload_Z = TRUE, verbose = TRUE)

# Get Hallmark gene sets
hallmark <- msigdbr(species = "Homo sapiens", category = "H")
get_hallmark_genes <- function(gs_name) {
  unique(hallmark$gene_symbol[hallmark$gs_name == gs_name])
}
```

## Defining the Viral Mimicry Event

### Get Component Gene Sets

``` r
# Interferon response (should go UP)
ifna <- get_hallmark_genes("HALLMARK_INTERFERON_ALPHA_RESPONSE")
ifng <- get_hallmark_genes("HALLMARK_INTERFERON_GAMMA_RESPONSE")

# Cell cycle (should go DOWN)
e2f <- get_hallmark_genes("HALLMARK_E2F_TARGETS")
g2m <- get_hallmark_genes("HALLMARK_G2M_CHECKPOINT")

message(sprintf("IFN-alpha: %d genes", length(ifna)))
message(sprintf("IFN-gamma: %d genes", length(ifng)))
message(sprintf("E2F targets: %d genes", length(e2f)))
message(sprintf("G2M checkpoint: %d genes", length(g2m)))
```

### Construct Multi-channel Query

``` r
# Define viral mimicry as 4 channels
q_vm <- list(
  IFNA = list(genes = ifna, sign = +1),  # IFN-alpha should be UP
  IFNG = list(genes = ifng, sign = +1),  # IFN-gamma should be UP
  E2F  = list(genes = e2f,  sign = -1),  # E2F targets should be DOWN
  G2M  = list(genes = g2m,  sign = -1)   # G2M checkpoint should be DOWN
)

# The sign parameter indicates expected direction:
# +1 = genes should be up-regulated in hits
# -1 = genes should be down-regulated in hits
```

## Running the Multi-channel Search

``` r
out_vm <- gpsa_search(db, q_vm, return_components = TRUE)

# View top hits by total score
message("=== Viral Mimicry Top 10 (Total Score) ===")
print(head(out_vm$res[, c("signature_id", "pbgene", "CellLineName",
                          "Score", "tau", "FDR_rank")], 10))
```

## Gated Filtering: Require All Channels

The total score sums across channels, but a high score could come from
just one or two strong channels. To find **true viral mimicry** hits, we
need all 4 channels to show the expected direction.

### Attach Channel Scores

``` r
res_vm <- gpsa_attach_channel_scores(
  out_vm$res, 
  out_vm$NES_raw, 
  out_vm$NES_signed
)

# Preview channel scores
message("=== Sample with Channel Scores ===")
print(head(res_vm[, c("signature_id", "pbgene", "Score",
                      "IFNA_signed", "IFNG_signed", 
                      "E2F_signed", "G2M_signed")], 5))
```

### Apply Gate: All Channels Positive

``` r
# Gate: all signed scores must be positive (direction matches expectation)
gated <- res_vm[
  res_vm$IFNA_signed > 0 & 
  res_vm$IFNG_signed > 0 & 
  res_vm$E2F_signed > 0 & 
  res_vm$G2M_signed > 0, 
]
gated <- gated[order(gated$Score, decreasing = TRUE), ]

message(sprintf("=== Gated Hits: %d signatures ===", nrow(gated)))

message("=== Viral Mimicry Gated Top 20 ===")
print(head(gated[, c("signature_id", "pbgene", "CellLineName", "Score",
                     "IFNA_signed", "IFNG_signed", 
                     "E2F_signed", "G2M_signed")], 20))
```

## Dual-mode Scoring for Multi-channel Events

For events with many channels, GPSA provides automatic dual-mode
scoring:

- **Phenocopy (Score_bi)**:
  [`min()`](https://rdrr.io/r/base/Extremes.html) across all channel
  signed scores
- **Discovery (Score_uni)**:
  [`max()`](https://rdrr.io/r/base/Extremes.html) across all channel
  signed scores

``` r
out_modes <- gpsa_postprocess_modes(
  out_vm,
  attach_channel_cols = TRUE,
  p_mode = "norm",
  filter_positive = FALSE  # Don't pre-filter, we'll examine both
)

# Phenocopy mode: all channels must agree
message("=== Phenocopy Mode Top 10 ===")
print(head(out_modes$res_similarity[, c("signature_id", "pbgene",
                                         "Score_bi", "Balance", "Driver")], 10))

# Discovery mode: at least one channel strong
message("=== Discovery Mode Top 10 ===")
print(head(out_modes$res_discovery[, c("signature_id", "pbgene",
                                        "Score_uni", "Balance", "Driver")], 10))
```

### Interpreting Dual-mode Results

| Scenario                             | Phenocopy Score | Discovery Score | Interpretation |
|--------------------------------------|-----------------|-----------------|----------------|
| All channels strong                  | High            | High            | True phenocopy |
| All channels weak                    | Low             | Low             | Not a hit      |
| One channel very strong, others weak | Low             | High            | Partial match  |
| Balanced moderate signal             | Moderate        | Moderate        | Possible hit   |

## GSEA Validation of Multi-channel Hits

Validate that a hit truly shows the expected pattern:

``` r
# Pick a gated hit to validate
if (nrow(gated) > 0) {
  hit <- gated$signature_id[1]
  gl_hit <- gpsa_get_signature(db, hit, dataset = "logFC")
  
  # Check IFN-alpha (should be enriched at TOP)
  gpsa_plot_gsea(gl_hit, ifna,
                 main = paste0(hit, " vs IFN-alpha (expect UP)"))
  
  # Check E2F targets (should be enriched at BOTTOM)
  gpsa_plot_gsea(gl_hit, e2f,
                 main = paste0(hit, " vs E2F targets (expect DOWN)"))
}
```

## Customizing the Event Definition

### Two-channel Events

For simpler events, use just 2 channels:

``` r
# Simple viral mimicry: IFN up + cell cycle down (combined gene sets)
ifn_combined <- unique(c(ifna, ifng))
cycle_combined <- unique(c(e2f, g2m))

q_simple <- list(
  IFN_UP = list(genes = ifn_combined, sign = +1),
  CYCLE_DOWN = list(genes = cycle_combined, sign = -1)
)

out_simple <- gpsa_search(db, q_simple, return_components = TRUE)
```

### Weighted Channels

Assign different weights to channels:

``` r
q_weighted <- list(
  IFNA = list(genes = ifna, sign = +1, weight = 2),   # Double weight
  IFNG = list(genes = ifng, sign = +1, weight = 1),
  E2F  = list(genes = e2f,  sign = -1, weight = 1),
  G2M  = list(genes = g2m,  sign = -1, weight = 1)
)
```

### Different Gate Thresholds

Use thresholds other than 0:

``` r
# Require stronger signal (signed > 0.5)
strong_gated <- res_vm[
  res_vm$IFNA_signed > 0.5 & 
  res_vm$IFNG_signed > 0.5 & 
  res_vm$E2F_signed > 0.5 & 
  res_vm$G2M_signed > 0.5, 
]
```

## Other Biological Events to Explore

The multi-channel approach works for any complex phenotype:

| Event           | Channels                                                       |
|-----------------|----------------------------------------------------------------|
| Senescence      | p53 pathway UP + Cell cycle DOWN + SASP UP                     |
| EMT             | Epithelial markers DOWN + Mesenchymal markers UP               |
| Stem-like state | Stemness genes UP + Differentiation markers DOWN               |
| Ferroptosis     | Iron metabolism UP + Lipid peroxidation UP + GPX4 targets DOWN |

## Best Practices

1.  **Start with known biology**: Define channels based on established
    pathway knowledge
2.  **Use gated filtering**: Total score can be misleading for
    multi-channel events
3.  **Examine channel balance**: `Driver` column shows which channel
    dominates
4.  **Validate with GSEA**: Confirm direction matches for each channel
5.  **Adjust thresholds**: Use stricter gates (\> 0.5) for higher
    confidence

## Summary

Multi-channel queries allow you to: - Define complex biological events
as combinations of gene sets - Specify expected direction (up or down)
for each component - Filter for hits where all channels satisfy the
expected direction - Use dual-mode scoring to distinguish phenocopy from
partial matches

## Session Info

``` r
sessionInfo()
```
