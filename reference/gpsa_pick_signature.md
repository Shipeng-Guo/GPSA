# Pick a signature by gene and cell line

Pick a signature by gene and cell line

## Usage

``` r
gpsa_pick_signature(db, pbgene, cellline = NULL, fallback = NULL)
```

## Arguments

- db:

  A gpsa_db object

- pbgene:

  Perturbation gene name

- cellline:

  Cell line name (optional)

- fallback:

  Fallback signature ID if not found

## Value

Signature ID or NULL
