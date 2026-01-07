# Compile full geneList to vector query weights

Converts a full geneList (named numeric vector) into query weights using
rankZ transformation.

## Usage

``` r
gpsa_compile_signature_vector(
  db,
  geneList,
  transform = c("rankZ", "zscore"),
  ties_method = "average",
  normalize = TRUE,
  min_size = 200
)
```

## Arguments

- db:

  A gpsa_db object

- geneList:

  Named numeric vector (gene symbols as names)

- transform:

  Transformation: "rankZ" or "zscore"

- ties_method:

  Ties method for ranking

- normalize:

  L2-normalize weights?

- min_size:

  Minimum overlap with DB

## Value

Compiled vector query object
