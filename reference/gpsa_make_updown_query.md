# Build Up/Down query from a signature geneList

Extracts top N and bottom N genes from a geneList to create an Up/Down
bidirectional query.

## Usage

``` r
gpsa_make_updown_query(geneList, n_each = 500)
```

## Arguments

- geneList:

  Named numeric vector (e.g., logFC values)

- n_each:

  Number of genes for each direction (default: 500)

## Value

List with Up and Down gene sets ready for gpsa_search

## Examples

``` r
if (FALSE) { # \dontrun{
db <- gpsa_load_db()
gl <- gpsa_get_signature(db, "D21455", dataset = "logFC")
query <- gpsa_make_updown_query(gl, n_each = 300)
out <- gpsa_search(db, query)
} # }
```
