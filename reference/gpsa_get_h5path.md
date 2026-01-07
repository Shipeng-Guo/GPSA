# Get GPSA H5 file path

Returns the HDF5 file path from options, environment variable, or
parameter. Priority: parameter \> options("GPSA.h5") \>
Sys.getenv("GPSA_H5")

## Usage

``` r
gpsa_get_h5path(h5file = NULL)
```

## Arguments

- h5file:

  Optional explicit path to HDF5 file

## Value

Character path to HDF5 file
