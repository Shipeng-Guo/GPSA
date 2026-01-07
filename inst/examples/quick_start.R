## ============================================================
## GPSA Quick Start Example
## ============================================================

library(GPSA)

## ============================================================
## 1) Set HDF5 Database Path
## ============================================================
## First-time users need to download the HDF5 database file,
## then set the path using one of these methods:

## Method 1: Set R option (recommended)
options(GPSA.h5 = "/path/to/gpsa_db.h5")

## Method 2: Set environment variable
# Sys.setenv(GPSA_H5 = "/path/to/gpsa_db.h5")

## Method 3: Specify path in each function call
# db <- gpsa_load_db(h5file = "/path/to/gpsa_db.h5")

## ============================================================
## 2) Load Database
## ============================================================

## Memory mode (recommended): Preload Z matrix, fast queries (~0.15s each)
## Load time: ~9 seconds, Memory usage: ~1GB
db <- gpsa_load_db(preload_Z = TRUE, verbose = TRUE)

## Disk mode: No preload, slower queries (~7s each)
## Suitable for limited memory or one-time queries
# db <- gpsa_load_db(preload_Z = FALSE, verbose = TRUE)

print(db)

## ============================================================
## 3) Single Gene Set Query
## ============================================================

## Example: Query perturbations related to a gene set
my_genes <- sample(db$genes, 300)  # Replace with your actual gene set
query <- list(MySet = list(genes = my_genes, sign = +1))

out <- gpsa_search(db, query, return_components = FALSE)

message("\n=== Single Gene Set Query Results ===")
print(head(out$res[, c("signature_id", "pbgene", "CellLineName", 
                        "Score", "tau", "FDR_rank")], 10))

## ============================================================
## 4) Up/Down Bidirectional Query
## ============================================================

## Extract logFC from a signature and build Up/Down query
sig_id <- db$sigs[1]  # Replace with actual signature ID
gl <- gpsa_get_signature(db, sig_id, dataset = "logFC")

## Create Up/Down query (top/bottom 300 genes)
query_ud <- gpsa_make_updown_query(gl, n_each = 300)

out_ud <- gpsa_search(db, query_ud, return_components = TRUE)

message("\n=== Up/Down Query - Top 10 (Most Similar) ===")
print(head(out_ud$res[, c("signature_id", "pbgene", "CellLineName",
                          "Score", "tau", "FDR_rank")], 10))

message("\n=== Up/Down Query - Bottom 10 (Most Opposite) ===")
print(tail(out_ud$res[, c("signature_id", "pbgene", "CellLineName",
                          "Score", "tau", "FDR_rank")], 10))

## ============================================================
## 5) Dual-mode Scoring (Phenocopy vs Discovery)
## ============================================================

out_ud <- gpsa_postprocess_modes(out_ud, 
                                  attach_channel_cols = TRUE,
                                  p_mode = "norm",
                                  filter_positive = TRUE)

message("\n=== Phenocopy Mode (All Channels Agree) ===")
print(head(out_ud$res_similarity[, c("signature_id", "pbgene",
                                      "Score_bi", "Up_signed", "Down_signed")], 10))

message("\n=== Discovery Mode (At Least One Channel Strong) ===")
print(head(out_ud$res_discovery[, c("signature_id", "pbgene",
                                     "Score_uni", "Driver", "Up_signed", "Down_signed")], 10))

## ============================================================
## 6) Batch Query
## ============================================================

## Prepare multiple queries
sig_ids <- db$sigs[1:3]
queries <- list()
for (sid in sig_ids) {
  gl <- gpsa_get_signature(db, sid, dataset = "logFC")
  queries[[sid]] <- gpsa_make_updown_query(gl, n_each = 300)
}

## Run batch search
bat <- gpsa_batch_search(db, queries, compute_tau = TRUE)

message("\n=== Batch Query Complete ===")
print(bat$timing)

## Query correlation matrix
message("\nQuery-to-Query Correlation:")
print(round(cor(bat$Score, use = "pairwise.complete.obs"), 3))

## Meta-analysis
meta <- gpsa_meta_stouffer(bat$z_ref)
message("\nMeta-analysis Top 10:")
print(head(meta, 10))

## ============================================================
## 7) GSEA Visualization
## ============================================================

## Validate a candidate hit
candidate <- out_ud$res$signature_id[2]  # Second-ranked hit
gl_cand <- gpsa_get_signature(db, candidate, dataset = "logFC")

## Plot enrichment for Up genes
gpsa_plot_gsea(gl_cand, query_ud$Up$genes, 
               main = paste0(candidate, " vs Query Up genes"))

## Plot enrichment for Down genes
gpsa_plot_gsea(gl_cand, query_ud$Down$genes,
               main = paste0(candidate, " vs Query Down genes"))

message("\n=== Quick Start Complete ===")
