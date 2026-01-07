# Adaptive Up/Down gene selection by energy

Selects up and down genes adaptively based on energy capture threshold,
avoiding arbitrary fixed cutoffs like 500.

## Usage

``` r
gpsa_make_updown_query_adaptive(geneList, eta = 0.9, min_n = 200, max_n = 2000)
```

## Arguments

- geneList:

  Named numeric vector

- eta:

  Energy capture threshold (default: 0.90)

- min_n:

  Minimum genes per direction

- max_n:

  Maximum genes per direction

## Value

List with Up, Down gene sets and .info
