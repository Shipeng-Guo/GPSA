# Install GPSA dependencies

Helper to install required and optional packages.

## Usage

``` r
gpsa_install_deps(
  install_minimal = TRUE,
  install_optional = TRUE,
  ask = FALSE,
  update = FALSE
)
```

## Arguments

- install_minimal:

  Install minimal dependencies (rhdf5, Matrix)?

- install_optional:

  Install optional packages (matrixStats, msigdbr, fgsea, ggplot2)?

- ask:

  Ask before installing?

- update:

  Update existing packages?

## Value

TRUE invisibly
