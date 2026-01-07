# Compute dual-mode scores (Phenocopy vs Discovery)

Computes Score_bi (strict phenocopy) and Score_uni (discovery) from
signed NES matrix.

## Usage

``` r
gpsa_compute_mode_scores(
  NES_signed,
  eps = 1e-09,
  p_mode = c("none", "norm"),
  adjust_method = "BH"
)
```

## Arguments

- NES_signed:

  Signed NES matrix (positive = desired direction)

- eps:

  Small value to avoid division by zero

- p_mode:

  P-value mode: "none" or "norm"

- adjust_method:

  P-value adjustment method

## Value

List with Score_bi, Score_uni, Balance, Driver, and optional p/FDR
