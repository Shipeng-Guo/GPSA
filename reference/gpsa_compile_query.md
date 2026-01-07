# Compile gene-set query

Compiles one or more gene sets into query weights for search.

## Usage

``` r
gpsa_compile_query(db, query_sets, min_size = 20)
```

## Arguments

- db:

  A gpsa_db object

- query_sets:

  Named list of gene sets, each with \$genes (character vector), \$sign
  (+1 or -1), and optional \$weight (default 1)

- min_size:

  Minimum genes after intersection with DB (default: 20)

## Value

Compiled query object

## Examples

``` r
if (FALSE) { # \dontrun{
db <- gpsa_load_db()
query <- list(
  Up = list(genes = up_genes, sign = +1),
  Down = list(genes = down_genes, sign = -1)
)
cq <- gpsa_compile_query(db, query)
} # }
```
