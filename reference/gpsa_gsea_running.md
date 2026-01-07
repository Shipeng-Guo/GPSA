# Compute GSEA running enrichment score

Compute GSEA running enrichment score

## Usage

``` r
gpsa_gsea_running(geneList, geneSet, gsea_param = 1)
```

## Arguments

- geneList:

  Named numeric vector (sorted descending)

- geneSet:

  Character vector of gene symbols

- gsea_param:

  GSEA weight parameter (default: 1)

## Value

List with ES, peak, running scores, hit positions
