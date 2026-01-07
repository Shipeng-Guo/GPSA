# Load GPSA Database

Loads an HDF5 database for querying. Optionally preloads the Z matrix
into memory for maximum query speed (~0.15s vs ~7s per query).

## Usage

``` r
gpsa_load_db(h5file = NULL, verbose = TRUE, preload_Z = FALSE)
```

## Arguments

- h5file:

  Path to HDF5 file (optional if set via options/env)

- verbose:

  Print messages?

- preload_Z:

  Preload Z matrix into memory? (Recommended for frequent queries)

## Value

A db object (list) for use with search functions

## Examples

``` r
if (FALSE) { # \dontrun{
# Set path once
options(GPSA.h5 = "/path/to/gpsa_db.h5")

# Load with memory preloading (fast queries)
db <- gpsa_load_db(preload_Z = TRUE)

# Or load without preloading (slower queries, less memory)
db <- gpsa_load_db(preload_Z = FALSE)
} # }
```
