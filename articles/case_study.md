# Complete Case Study: From Gene Sets to Discovery

## Overview

This case study demonstrates a complete GPSA workflow, from simple gene
set queries to complex multi-channel event discovery. We’ll use a
realistic biological scenario to showcase the full power of GPSA.

**Scenario**: You have identified a set of genes differentially
expressed in your experiment and want to find perturbations in public
databases that produce similar or related effects.

## Tutorial Roadmap

GPSA supports increasingly sophisticated query types. Here’s how they
build on each other:

    Single Gene Set Query (Get Started)
             |
             v
    Up/Down Bidirectional Query -----> Find similar/opposite perturbations
             |
             v
    Multi-channel Event Definition --> Complex phenotype queries
             |
             v
    Batch Queries + Meta-analysis ---> Compare multiple signatures
             |
             v
    Full Vector Search --------------> Complete expression profile queries

## Detailed Tutorials

Each tutorial covers a specific capability in depth:

### 1. Up/Down Bidirectional Query

**[Tutorial: Up/Down
Query](https://shipeng-guo.github.io/GPSA/articles/tutorial_updown.md)**

When you have a reference perturbation signature (with both up- and
down-regulated genes), you can find other perturbations that produce
**similar** or **opposite** effects.

Key concepts covered: - Fixed gene selection (`n_each = 500`) - Adaptive
gene selection based on energy threshold - Phenocopy vs Discovery
scoring modes - Channel-level score interpretation

### 2. Multi-channel Event Definition

**[Tutorial: Viral
Mimicry](https://shipeng-guo.github.io/GPSA/articles/tutorial_viral_mimicry.md)**

Define complex biological events by combining multiple gene sets with
different expected directions. Example: “viral mimicry” state requires
both interferon activation AND cell cycle suppression.

Key concepts covered: - Multi-channel query construction - Gated
filtering (all channels must satisfy direction) - Dual-mode scoring
(Phenocopy vs Discovery) - Biological interpretation of multi-channel
hits

### 3. Batch Queries and Meta-analysis

**[Tutorial: Batch
Queries](https://shipeng-guo.github.io/GPSA/articles/tutorial_batch.md)**

When you have multiple related signatures (e.g., replicates or related
conditions), run them together for efficiency and perform meta-analysis.

Key concepts covered: -
[`gpsa_batch_search()`](https://shipeng-guo.github.io/GPSA/reference/gpsa_batch_search.md)
for multiple queries - Query-to-query correlation analysis - Stouffer’s
method for meta-analysis - Identifying consensus hits across queries

### 4. Full Vector Search

**[Tutorial: Vector
Search](https://shipeng-guo.github.io/GPSA/articles/tutorial_vector_search.md)**

Use complete gene expression profiles (e.g., full logFC vectors) as
queries, leveraging the entire transcriptomic information rather than
just top/bottom genes.

Key concepts covered: -
[`gpsa_search_vector()`](https://shipeng-guo.github.io/GPSA/reference/gpsa_search_vector.md)
with rankZ transformation - Adaptive sparsification (energy-based gene
selection) - Compile information interpretation - Comparison with
Up/Down approach

## Quick Reference: Choosing the Right Query Type

| Scenario                                  | Recommended Approach  | Tutorial                                                                               |
|-------------------------------------------|-----------------------|----------------------------------------------------------------------------------------|
| Query with a pathway gene set             | Single gene set query | [Get Started](https://shipeng-guo.github.io/GPSA/articles/GPSA.md)                     |
| Find perturbations similar to a reference | Up/Down query         | [Up/Down](https://shipeng-guo.github.io/GPSA/articles/tutorial_updown.md)              |
| Complex phenotype (multiple pathways)     | Multi-channel query   | [Viral Mimicry](https://shipeng-guo.github.io/GPSA/articles/tutorial_viral_mimicry.md) |
| Multiple related signatures               | Batch + meta-analysis | [Batch](https://shipeng-guo.github.io/GPSA/articles/tutorial_batch.md)                 |
| Complete expression profile               | Vector search         | [Vector Search](https://shipeng-guo.github.io/GPSA/articles/tutorial_vector_search.md) |

## Performance Tips

1.  **Use memory mode** for interactive analysis:

    ``` r
    db <- gpsa_load_db(preload_Z = TRUE)
    ```

2.  **Batch related queries** to share I/O overhead:

    ``` r
    bat <- gpsa_batch_search(db, list(Q1 = q1, Q2 = q2, Q3 = q3))
    ```

3.  **Use adaptive gene selection** to avoid arbitrary cutoffs:

    ``` r
    query <- gpsa_make_updown_query_adaptive(gl, eta = 0.90)
    ```

4.  **Apply dual-mode scoring** to distinguish phenocopy from discovery:

    ``` r
    out <- gpsa_postprocess_modes(out, p_mode = "norm")
    head(out$res_similarity)  # Phenocopy hits
    head(out$res_discovery)   # Discovery hits
    ```

## API Quick Reference

### Core Functions

| Function                                                                                     | Description              |
|----------------------------------------------------------------------------------------------|--------------------------|
| [`gpsa_load_db()`](https://shipeng-guo.github.io/GPSA/reference/gpsa_load_db.md)             | Load HDF5 database       |
| [`gpsa_search()`](https://shipeng-guo.github.io/GPSA/reference/gpsa_search.md)               | Gene-set query           |
| [`gpsa_search_vector()`](https://shipeng-guo.github.io/GPSA/reference/gpsa_search_vector.md) | Full geneList query      |
| [`gpsa_batch_search()`](https://shipeng-guo.github.io/GPSA/reference/gpsa_batch_search.md)   | Multiple queries at once |

### Query Helpers

| Function                                                                                                               | Description                    |
|------------------------------------------------------------------------------------------------------------------------|--------------------------------|
| [`gpsa_make_updown_query()`](https://shipeng-guo.github.io/GPSA/reference/gpsa_make_updown_query.md)                   | Create Up/Down query (fixed n) |
| [`gpsa_make_updown_query_adaptive()`](https://shipeng-guo.github.io/GPSA/reference/gpsa_make_updown_query_adaptive.md) | Adaptive gene selection        |
| [`gpsa_get_signature()`](https://shipeng-guo.github.io/GPSA/reference/gpsa_get_signature.md)                           | Extract signature vector       |

### Post-processing

| Function                                                                                                     | Description                  |
|--------------------------------------------------------------------------------------------------------------|------------------------------|
| [`gpsa_postprocess_modes()`](https://shipeng-guo.github.io/GPSA/reference/gpsa_postprocess_modes.md)         | Dual-mode scoring            |
| [`gpsa_attach_channel_scores()`](https://shipeng-guo.github.io/GPSA/reference/gpsa_attach_channel_scores.md) | Add per-channel columns      |
| [`gpsa_meta_stouffer()`](https://shipeng-guo.github.io/GPSA/reference/gpsa_meta_stouffer.md)                 | Meta-analysis across queries |

### Visualization

| Function                                                                                   | Description                  |
|--------------------------------------------------------------------------------------------|------------------------------|
| [`gpsa_plot_gsea()`](https://shipeng-guo.github.io/GPSA/reference/gpsa_plot_gsea.md)       | GSEA running enrichment plot |
| [`gpsa_gsea_running()`](https://shipeng-guo.github.io/GPSA/reference/gpsa_gsea_running.md) | Compute running scores       |

## Session Info

``` r
sessionInfo()
```
