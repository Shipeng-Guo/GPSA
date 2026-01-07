#' Build GPSA HDF5 Database
#'
#' Constructs an HDF5 database from a diffTable RDS file containing
#' genes x signatures logFC matrix. The database stores rankZ-transformed
#' values for fast similarity search.
#'
#' @param rds_file Path to RDS file containing diffTable (genes x signatures)
#' @param h5file Output HDF5 file path
#' @param meta_file Optional path to metadata RDS (must contain signature_id, pbgene, CellLineName)
#' @param gene_col Column name for gene symbols (default: first column)
#' @param na_fill_rank Value to fill NA before ranking (default: 0)
#' @param ties_method Ties method for ranking (default: "average")
#' @param block_s Block size for writing (default: 128)
#' @param chunk_g_Z Chunk size for genes in Z matrix (default: 1024)
#' @param chunk_s_Z Chunk size for signatures in Z matrix (default: 512)
#' @param chunk_g_logFC Chunk size for genes in logFC matrix (default: 2048)
#' @param compress_level Compression level 0-9 (default: 1, lower = faster read)
#' @param store_logFC Store raw logFC values? (default: TRUE)
#' @param use_float32 Use float32 storage? (default: TRUE)
#' @param verbose Print progress? (default: TRUE)
#'
#' @return TRUE invisibly on success
#' @export
#'
#' @examples
#' \dontrun{
#' gpsa_build_db(
#'   rds_file = "diffTable.rds",
#'   h5file = "gpsa_db.h5",
#'   meta_file = "metadata.rds"
#' )
#' }
gpsa_build_db <- function(
    rds_file,
    h5file,
    meta_file = NULL,
    gene_col = NULL,
    na_fill_rank = 0,
    ties_method = "average",
    block_s = 128,
    chunk_g_Z = 1024,
    chunk_s_Z = 512,
    chunk_g_logFC = 2048,
    compress_level = 1,
    store_logFC = TRUE,
    use_float32 = TRUE,
    verbose = TRUE
) {
  if (!file.exists(rds_file)) stop("RDS not found: ", rds_file)
  dir.create(dirname(h5file), showWarnings = FALSE, recursive = TRUE)
  
  if (verbose) message("Loading diffTable RDS ...")
  dt <- readRDS(rds_file)
  dt <- .gpsa_drop_datatable_class_nocopy(dt)
  if (!inherits(dt, "data.frame")) stop("diffTable must be data.frame-like.")
  
  if (is.null(gene_col)) gene_col <- names(dt)[1]
  if (!gene_col %in% names(dt)) stop("gene_col not found: ", gene_col)
  
  genes <- as.character(dt[[gene_col]])
  sigs  <- setdiff(names(dt), gene_col)
  
  G <- length(genes)
  N <- length(sigs)
  if (verbose) message(sprintf("Input: %d genes x %d signatures", G, N))
  
  if (anyDuplicated(genes)) {
    warning("Duplicate gene symbols detected. Consider upstream deduplication.")
  }
  
  if (file.exists(h5file)) {
    if (verbose) message("Removing existing HDF5: ", h5file)
    file.remove(h5file)
  }
  
  if (verbose) message("Creating HDF5: ", h5file)
  rhdf5::h5createFile(h5file)
  
  H5float <- NULL
  if (use_float32) {
    H5float <- tryCatch(rhdf5::h5const("H5T_NATIVE_FLOAT"), error = function(e) NULL)
    if (is.null(H5float) && verbose) message("NOTE: float32 not available; fallback to double.")
  }
  
  ## /Z dataset
  if (is.null(H5float)) {
    rhdf5::h5createDataset(h5file, "Z", dims = c(G, N),
                           chunk = c(min(chunk_g_Z, G), min(chunk_s_Z, N)),
                           level = compress_level, storage.mode = "double")
  } else {
    rhdf5::h5createDataset(h5file, "Z", dims = c(G, N),
                           chunk = c(min(chunk_g_Z, G), min(chunk_s_Z, N)),
                           level = compress_level, H5type = H5float)
  }
  
  ## /logFC dataset (optional)
  if (store_logFC) {
    if (is.null(H5float)) {
      rhdf5::h5createDataset(h5file, "logFC", dims = c(G, N),
                             chunk = c(min(chunk_g_logFC, G), 1),
                             level = compress_level, storage.mode = "double")
    } else {
      rhdf5::h5createDataset(h5file, "logFC", dims = c(G, N),
                             chunk = c(min(chunk_g_logFC, G), 1),
                             level = compress_level, H5type = H5float)
    }
  }
  
  ## Names & definition
  rhdf5::h5write(genes, h5file, "genes")
  rhdf5::h5write(sigs, h5file, "sigs")
  rhdf5::h5write(as.character(gene_col), h5file, "gene_col")
  rhdf5::h5write("cmap_rankZ_qnorm", h5file, "Z_definition")
  
  ## Meta vectors under /meta
  rhdf5::h5createGroup(h5file, "meta")
  pbgene <- rep(NA_character_, N)
  cell <- rep(NA_character_, N)
  
  if (!is.null(meta_file) && file.exists(meta_file)) {
    if (verbose) message("Loading metadata: ", meta_file)
    m <- readRDS(meta_file)
    m <- .gpsa_drop_datatable_class_nocopy(m)
    if (!inherits(m, "data.frame")) stop("meta_file must be data.frame-like.")
    
    if (!"signature_id" %in% names(m)) {
      if ("Dataset" %in% names(m)) names(m)[names(m) == "Dataset"] <- "signature_id"
    }
    
    need <- c("signature_id", "pbgene", "CellLineName")
    miss <- setdiff(need, names(m))
    if (length(miss) > 0) stop("Metadata missing columns: ", paste(miss, collapse = ", "))
    
    idx <- match(sigs, as.character(m$signature_id))
    pbgene <- as.character(m$pbgene[idx])
    cell <- as.character(m$CellLineName[idx])
  }
  rhdf5::h5write(pbgene, h5file, "meta/pbgene")
  rhdf5::h5write(cell, h5file, "meta/CellLineName")
  
  if (verbose) {
    message("Building /Z (rankZ) and /logFC (optional) ...")
    if (!requireNamespace("matrixStats", quietly = TRUE)) {
      message("NOTE: Install matrixStats for faster building.")
    }
  }
  
  for (j in seq(1, N, by = block_s)) {
    j_end <- min(N, j + block_s - 1)
    cols <- sigs[j:j_end]
    if (verbose) message(sprintf("  signatures %d-%d / %d", j, j_end, N))
    
    block <- as.matrix(dt[, cols, drop = FALSE])
    
    if (store_logFC) {
      rhdf5::h5write(block, h5file, "logFC", index = list(1:G, j:j_end))
    }
    
    if (anyNA(block)) block[is.na(block)] <- na_fill_rank
    
    if (requireNamespace("matrixStats", quietly = TRUE)) {
      rnk <- matrixStats::colRanks(block, ties.method = ties_method, preserveShape = TRUE)
    } else {
      rnk <- apply(block, 2, rank, ties.method = ties_method)
    }
    
    p <- (rnk - 0.5) / nrow(block)
    Zblk <- qnorm(p)
    Zblk[!is.finite(Zblk)] <- 0
    
    rhdf5::h5write(Zblk, h5file, "Z", index = list(1:G, j:j_end))
    
    rm(block, rnk, p, Zblk)
    gc(FALSE)
  }
  
  if (verbose) message("DB build finished: ", h5file)
  invisible(TRUE)
}

#' Install GPSA dependencies
#'
#' Helper to install required and optional packages.
#'
#' @param install_minimal Install minimal dependencies (rhdf5, Matrix)?
#' @param install_optional Install optional packages (matrixStats, msigdbr, fgsea, ggplot2)?
#' @param ask Ask before installing?
#' @param update Update existing packages?
#' @return TRUE invisibly
#' @export
gpsa_install_deps <- function(
    install_minimal = TRUE,
    install_optional = TRUE,
    ask = FALSE,
    update = FALSE
) {
  cran_min <- c("Matrix")
  bioc_min <- c("rhdf5")
  
  cran_opt <- c("matrixStats", "msigdbr", "ggplot2")
  bioc_opt <- c("fgsea")
  
  if (install_minimal) {
    for (p in cran_min) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    for (p in bioc_min) if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = ask, update = update)
  }
  
  if (install_optional) {
    for (p in cran_opt) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    for (p in bioc_opt) if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = ask, update = update)
  }
  
  invisible(TRUE)
}

