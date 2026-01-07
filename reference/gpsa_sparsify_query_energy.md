# Adaptive sparsification by energy capture

Reduces query size while preserving a specified fraction of signal
energy.

## Usage

``` r
gpsa_sparsify_query_energy(
  union_idx,
  w_union,
  eta = 0.9,
  min_k = 200,
  max_k = NULL
)
```

## Arguments

- union_idx:

  Gene indices

- w_union:

  Weights

- eta:

  Energy capture threshold (default: 0.90)

- min_k:

  Minimum genes to keep

- max_k:

  Maximum genes to keep

## Value

Sparsified query object
