# Package index

## Core Functions

Main search and database functions

- [`gpsa_load_db()`](https://shipeng-guo.github.io/GPSA/reference/gpsa_load_db.md)
  : Load GPSA Database
- [`print(`*`<gpsa_db>`*`)`](https://shipeng-guo.github.io/GPSA/reference/print.gpsa_db.md)
  : Print method for gpsa_db
- [`gpsa_search()`](https://shipeng-guo.github.io/GPSA/reference/gpsa_search.md)
  : Search with gene-set query
- [`gpsa_search_vector()`](https://shipeng-guo.github.io/GPSA/reference/gpsa_search_vector.md)
  : Search with vector query (full geneList)
- [`gpsa_batch_search()`](https://shipeng-guo.github.io/GPSA/reference/gpsa_batch_search.md)
  : Batch search (multiple queries at once)
- [`gpsa_build_db()`](https://shipeng-guo.github.io/GPSA/reference/gpsa_build_db.md)
  : Build GPSA HDF5 Database

## Query Helpers

Functions to construct queries

- [`gpsa_make_updown_query()`](https://shipeng-guo.github.io/GPSA/reference/gpsa_make_updown_query.md)
  : Build Up/Down query from a signature geneList
- [`gpsa_make_updown_query_adaptive()`](https://shipeng-guo.github.io/GPSA/reference/gpsa_make_updown_query_adaptive.md)
  : Adaptive Up/Down gene selection by energy
- [`gpsa_compile_query()`](https://shipeng-guo.github.io/GPSA/reference/gpsa_compile_query.md)
  : Compile gene-set query
- [`gpsa_compile_signature_vector()`](https://shipeng-guo.github.io/GPSA/reference/gpsa_compile_signature_vector.md)
  : Compile full geneList to vector query weights
- [`gpsa_sparsify_query_energy()`](https://shipeng-guo.github.io/GPSA/reference/gpsa_sparsify_query_energy.md)
  : Adaptive sparsification by energy capture

## Signature Extraction

Extract and prepare signature data

- [`gpsa_get_signature()`](https://shipeng-guo.github.io/GPSA/reference/gpsa_get_signature.md)
  : Extract a signature vector from the database
- [`gpsa_pick_signature()`](https://shipeng-guo.github.io/GPSA/reference/gpsa_pick_signature.md)
  : Pick a signature by gene and cell line
- [`gpsa_prepare_geneList()`](https://shipeng-guo.github.io/GPSA/reference/gpsa_prepare_geneList.md)
  : Prepare a geneList for GSEA/search
- [`gpsa_collapse_duplicates_maxabs()`](https://shipeng-guo.github.io/GPSA/reference/gpsa_collapse_duplicates_maxabs.md)
  : Collapse duplicate gene names by max absolute value

## Post-processing

Score computation and result processing

- [`gpsa_postprocess_modes()`](https://shipeng-guo.github.io/GPSA/reference/gpsa_postprocess_modes.md)
  : Postprocess search results with dual-mode scoring
- [`gpsa_compute_mode_scores()`](https://shipeng-guo.github.io/GPSA/reference/gpsa_compute_mode_scores.md)
  : Compute dual-mode scores (Phenocopy vs Discovery)
- [`gpsa_attach_channel_scores()`](https://shipeng-guo.github.io/GPSA/reference/gpsa_attach_channel_scores.md)
  : Attach per-channel scores to result data.frame
- [`gpsa_score_stats()`](https://shipeng-guo.github.io/GPSA/reference/gpsa_score_stats.md)
  : Compute score statistics
- [`gpsa_meta_stouffer()`](https://shipeng-guo.github.io/GPSA/reference/gpsa_meta_stouffer.md)
  : Meta-analysis across queries (Stouffer method)
- [`gpsa_rank_in()`](https://shipeng-guo.github.io/GPSA/reference/gpsa_rank_in.md)
  : Get rank of a signature in a result data.frame

## Visualization

GSEA plots and enrichment visualization

- [`gpsa_plot_gsea()`](https://shipeng-guo.github.io/GPSA/reference/gpsa_plot_gsea.md)
  : Plot GSEA running enrichment score
- [`gpsa_gsea_running()`](https://shipeng-guo.github.io/GPSA/reference/gpsa_gsea_running.md)
  : Compute GSEA running enrichment score
- [`gpsa_fgsea_updown()`](https://shipeng-guo.github.io/GPSA/reference/gpsa_fgsea_updown.md)
  : GSEA-style validation with fgsea

## Utilities

Configuration and helper functions

- [`gpsa_get_h5path()`](https://shipeng-guo.github.io/GPSA/reference/gpsa_get_h5path.md)
  : Get GPSA H5 file path
- [`gpsa_set_h5path()`](https://shipeng-guo.github.io/GPSA/reference/gpsa_set_h5path.md)
  : Set GPSA H5 file path
- [`gpsa_install_deps()`](https://shipeng-guo.github.io/GPSA/reference/gpsa_install_deps.md)
  : Install GPSA dependencies
