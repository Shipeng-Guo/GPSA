# GSEA-style validation with fgsea

Runs fgsea on up and down gene sets for validation.

## Usage

``` r
gpsa_fgsea_updown(
  geneList,
  up_genes,
  down_genes,
  nperm = 2000,
  minSize = 10,
  maxSize = 5000
)
```

## Arguments

- geneList:

  Named numeric vector

- up_genes:

  Up-regulated gene set

- down_genes:

  Down-regulated gene set

- nperm:

  Number of permutations

- minSize:

  Minimum gene set size

- maxSize:

  Maximum gene set size

## Value

fgsea result data.frame
