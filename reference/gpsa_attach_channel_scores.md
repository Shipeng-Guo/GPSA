# Attach per-channel scores to result data.frame

Adds raw and signed NES scores for each channel to the result table.

## Usage

``` r
gpsa_attach_channel_scores(
  res,
  NES_raw,
  NES_signed,
  suffix_raw = "_raw",
  suffix_signed = "_signed"
)
```

## Arguments

- res:

  Result data.frame from gpsa_search

- NES_raw:

  Raw NES matrix (channels x signatures)

- NES_signed:

  Signed NES matrix

- suffix_raw:

  Suffix for raw score columns

- suffix_signed:

  Suffix for signed score columns

## Value

Updated data.frame with channel columns
