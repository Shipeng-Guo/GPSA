# Batch search (multiple queries at once)

Performs multiple queries efficiently by sharing I/O.

## Usage

``` r
gpsa_batch_search(
  db,
  query_list,
  block_s = NULL,
  robust = FALSE,
  compute_tau = TRUE
)
```

## Arguments

- db:

  A gpsa_db object

- query_list:

  Named list of queries (each is a query_sets list)

- block_s:

  Block size for reading

- robust:

  Use robust statistics?

- compute_tau:

  Compute tau scores?

## Value

List with Score matrix and statistics matrices

## Examples

``` r
if (FALSE) { # \dontrun{
db <- gpsa_load_db(preload_Z = TRUE)

queries <- list(
  Q1 = list(Up = list(genes = genes1, sign = +1)),
  Q2 = list(Up = list(genes = genes2, sign = +1))
)
bat <- gpsa_batch_search(db, queries)
} # }
```
