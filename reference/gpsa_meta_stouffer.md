# Meta-analysis across queries (Stouffer method)

Combines z-scores from multiple queries using Stouffer's method.

## Usage

``` r
gpsa_meta_stouffer(z_mat, weights = NULL)
```

## Arguments

- z_mat:

  Matrix of z-scores (signatures x queries)

- weights:

  Optional weights for each query

## Value

data.frame with meta z-scores, p-values, and FDR
