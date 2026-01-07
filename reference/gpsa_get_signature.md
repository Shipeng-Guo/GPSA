# Extract a signature vector from the database

Fast extraction of a single signature's logFC or Z values.

## Usage

``` r
gpsa_get_signature(
  db,
  signature_id,
  dataset = c("logFC", "Z"),
  sort_decreasing = FALSE
)
```

## Arguments

- db:

  A gpsa_db object

- signature_id:

  Signature ID to extract

- dataset:

  Which dataset: "logFC" or "Z"

- sort_decreasing:

  Sort by value decreasing?

## Value

Named numeric vector (gene symbols as names)

## Examples

``` r
if (FALSE) { # \dontrun{
db <- gpsa_load_db()
gl <- gpsa_get_signature(db, "D21455", dataset = "logFC")
head(sort(gl, decreasing = TRUE))
} # }
```
