#' Load GPSA Database
#'
#' Loads an HDF5 database for querying. Optionally preloads the Z matrix
#' into memory for maximum query speed (~0.15s vs ~7s per query).
#'
#' @param h5file Path to HDF5 file (optional if set via options/env)
#' @param verbose Print messages?
#' @param preload_Z Preload Z matrix into memory? (Recommended for frequent queries)
#'
#' @return A db object (list) for use with search functions
#' @export
#'
#' @examples
#' \dontrun{
#' # Set path once
#' options(GPSA.h5 = "/path/to/gpsa_db.h5")
#'
#' # Load with memory preloading (fast queries)
#' db <- gpsa_load_db(preload_Z = TRUE)
#'
#' # Or load without preloading (slower queries, less memory)
#' db <- gpsa_load_db(preload_Z = FALSE)
#' }
gpsa_load_db <- function(h5file = NULL, verbose = TRUE, preload_Z = FALSE) {
  h5file <- gpsa_get_h5path(h5file)
  if (!file.exists(h5file)) stop("HDF5 not found: ", h5file)
  
  genes <- rhdf5::h5read(h5file, "genes")
  sigs <- rhdf5::h5read(h5file, "sigs")
  
  pbgene <- tryCatch(rhdf5::h5read(h5file, "meta/pbgene"),
                     error = function(e) rep(NA_character_, length(sigs)))
  cell <- tryCatch(rhdf5::h5read(h5file, "meta/CellLineName"),
                   error = function(e) rep(NA_character_, length(sigs)))
  
  Z_mem <- NULL
  if (isTRUE(preload_Z)) {
    if (verbose) message("Preloading /Z into memory ...")
    Z_mem <- tryCatch({
      Z <- rhdf5::h5read(h5file, "Z", drop = FALSE, read.attributes = FALSE)
      if (!is.matrix(Z)) Z <- .gpsa_as_matrix2d(Z, nrow = length(genes), ncol = length(sigs))
      storage.mode(Z) <- "double"
      Z
    }, error = function(e) {
      if (verbose) message("Preload failed: ", conditionMessage(e))
      NULL
    })
    if (!is.null(Z_mem) && verbose) {
      megs <- round(object.size(Z_mem) / (1024^2), 1)
      message(sprintf("Preloaded /Z (%.1f MB) into RAM; subsequent queries avoid disk I/O.", megs))
    }
  }
  
  meta <- data.frame(
    signature_id = as.character(sigs),
    pbgene = as.character(pbgene),
    CellLineName = as.character(cell),
    stringsAsFactors = FALSE
  )
  
  db <- list(
    h5file = h5file,
    genes = as.character(genes),
    sigs = as.character(sigs),
    G = length(genes),
    N = length(sigs),
    gene2idx = setNames(seq_along(genes), as.character(genes)),
    sig2idx = setNames(seq_along(sigs), as.character(sigs)),
    meta = meta
  )
  if (!is.null(Z_mem)) db$Z_mem <- Z_mem
  
  if (verbose) {
    zd <- tryCatch(rhdf5::h5read(h5file, "Z_definition"), error = function(e) NA_character_)
    message(sprintf("Loaded DB: %d genes x %d signatures; Z_definition=%s", db$G, db$N, zd))
    if (!.gpsa_h5exists(h5file, "logFC")) {
      message("NOTE: /logFC not found. Cannot extract logFC geneList for GSEA plots.")
    }
  }
  
  class(db) <- c("gpsa_db", "list")
  db
}

#' Print method for gpsa_db
#' @param x A gpsa_db object
#' @param ... Ignored
#' @export
print.gpsa_db <- function(x, ...) {
  cat("GPSA Database\n")
  cat(sprintf("  Genes: %d\n", x$G))
  cat(sprintf("  Signatures: %d\n", x$N))
  cat(sprintf("  H5 file: %s\n", x$h5file))
  cat(sprintf("  Memory preloaded: %s\n", ifelse(is.null(x$Z_mem), "No", "Yes")))
  invisible(x)
}

#' Extract a signature vector from the database
#'
#' Fast extraction of a single signature's logFC or Z values.
#'
#' @param db A gpsa_db object
#' @param signature_id Signature ID to extract
#' @param dataset Which dataset: "logFC" or "Z"
#' @param sort_decreasing Sort by value decreasing?
#'
#' @return Named numeric vector (gene symbols as names)
#' @export
#'
#' @examples
#' \dontrun{
#' db <- gpsa_load_db()
#' gl <- gpsa_get_signature(db, "D21455", dataset = "logFC")
#' head(sort(gl, decreasing = TRUE))
#' }
gpsa_get_signature <- function(db, signature_id, dataset = c("logFC", "Z"), sort_decreasing = FALSE) {
  dataset <- match.arg(dataset)
  if (!.gpsa_h5exists(db$h5file, dataset)) stop("Dataset /", dataset, " not found in HDF5.")
  
  j <- unname(db$sig2idx[signature_id])
  if (is.na(j)) stop("signature_id not found: ", signature_id)
  
  if (!is.null(db$Z_mem) && dataset == "Z") {
    v <- as.numeric(db$Z_mem[, j])
  } else {
    v <- rhdf5::h5read(db$h5file, dataset, index = list(1:db$G, j), drop = FALSE, read.attributes = FALSE)
    if (is.matrix(v)) v <- as.numeric(v[, 1]) else v <- as.numeric(v)
  }
  
  if (length(v) != length(db$genes)) stop("Length mismatch: signature vector != genes length.")
  v <- setNames(v, db$genes)
  
  if (sort_decreasing) v <- sort(v, decreasing = TRUE)
  v
}

