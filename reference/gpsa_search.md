# Search with gene-set query

Performs similarity search using one or more gene sets.

## Usage

``` r
gpsa_search(
  db,
  query_sets,
  return_components = TRUE,
  block_s = NULL,
  robust = TRUE
)
```

## Arguments

- db:

  A gpsa_db object

- query_sets:

  Named list of gene sets (see gpsa_compile_query)

- return_components:

  Return per-channel NES scores?

- block_s:

  Block size for reading (NULL = auto)

- robust:

  Use robust statistics?

## Value

List with res (data.frame), compile info, and optionally
NES_raw/NES_signed

## Examples

``` r
if (FALSE) { # \dontrun{
db <- gpsa_load_db(preload_Z = TRUE)

# Single gene set query
query <- list(MySet = list(genes = my_genes, sign = +1))
out <- gpsa_search(db, query)
head(out$res)

# Up/Down query
query <- list(
  Up = list(genes = up_genes, sign = +1),
  Down = list(genes = down_genes, sign = -1)
)
out <- gpsa_search(db, query, return_components = TRUE)
} # }
```
