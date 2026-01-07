# Search with vector query (full geneList)

Performs similarity search using a full gene expression profile.

## Usage

``` r
gpsa_search_vector(
  db,
  geneList,
  eta = 0.9,
  min_k = 300,
  max_k = 5000,
  transform = c("rankZ", "zscore"),
  ties_method = "average",
  normalize = TRUE,
  block_s = NULL,
  robust = TRUE
)
```

## Arguments

- db:

  A gpsa_db object

- geneList:

  Named numeric vector (gene symbols as names)

- eta:

  Energy capture threshold for sparsification (NULL = no sparsification)

- min_k:

  Minimum genes to use

- max_k:

  Maximum genes to use

- transform:

  Transformation: "rankZ" or "zscore"

- ties_method:

  Ties method for ranking

- normalize:

  L2-normalize weights?

- block_s:

  Block size for reading

- robust:

  Use robust statistics?

## Value

List with res (data.frame) and compile info
