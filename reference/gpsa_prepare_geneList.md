# Prepare a geneList for GSEA/search

Prepare a geneList for GSEA/search

## Usage

``` r
gpsa_prepare_geneList(
  x,
  decreasing = TRUE,
  collapse_duplicates = TRUE,
  na.rm = TRUE
)
```

## Arguments

- x:

  Named numeric vector (gene symbols as names)

- decreasing:

  Sort decreasing?

- collapse_duplicates:

  Collapse duplicates by max\|value\|?

- na.rm:

  Remove NA values?

## Value

Sorted named numeric vector
