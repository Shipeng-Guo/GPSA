# Postprocess search results with dual-mode scoring

Adds phenocopy (Score_bi) and discovery (Score_uni) scores to search
results. Creates two ranked tables: res_similarity and res_discovery.

## Usage

``` r
gpsa_postprocess_modes(
  out,
  attach_channel_cols = TRUE,
  p_mode = c("none", "norm"),
  eps = 1e-09,
  filter_positive = TRUE
)
```

## Arguments

- out:

  Output from gpsa_search with return_components=TRUE

- attach_channel_cols:

  Attach per-channel score columns?

- p_mode:

  P-value computation mode

- eps:

  Epsilon for Balance calculation

- filter_positive:

  Filter to positive scores only?

## Value

Updated output with res, res_similarity, res_discovery

## Examples

``` r
if (FALSE) { # \dontrun{
db <- gpsa_load_db(preload_Z = TRUE)
gl <- gpsa_get_signature(db, "D21455", dataset = "logFC")
query <- gpsa_make_updown_query(gl)
out <- gpsa_search(db, query, return_components = TRUE)
out <- gpsa_postprocess_modes(out, p_mode = "norm")

# Phenocopy hits (all channels positive)
head(out$res_similarity)

# Discovery hits (at least one channel strong)
head(out$res_discovery)
} # }
```
