# Plot GSEA running enrichment score

Plot GSEA running enrichment score

## Usage

``` r
gpsa_plot_gsea(
  geneList,
  geneSet,
  main = NULL,
  method = c("base", "fgsea"),
  gsea_param = 1
)
```

## Arguments

- geneList:

  Named numeric vector

- geneSet:

  Character vector of gene symbols

- main:

  Plot title

- method:

  "base" (default) or "fgsea" if available

- gsea_param:

  GSEA weight parameter

## Value

Invisibly returns the running score object (base) or ggplot (fgsea)

## Examples

``` r
if (FALSE) { # \dontrun{
db <- gpsa_load_db()
gl <- gpsa_get_signature(db, "D21455", dataset = "logFC")
gpsa_plot_gsea(gl, my_gene_set, main = "My Gene Set")
} # }
```
