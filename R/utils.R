#' @importFrom stats qnorm pnorm p.adjust sd median mad
#' @importFrom graphics plot abline rug mtext
#' @importFrom grDevices dev.new
#' @importFrom utils head tail
NULL

## Internal utilities (not exported)

.gpsa_as_matrix2d <- function(x, nrow, ncol) {
  if (is.matrix(x) && all(dim(x) == c(nrow, ncol))) return(x)
  matrix(as.numeric(x), nrow = nrow, ncol = ncol)
}

.gpsa_pick_block_s <- function(block_s, N, prefer = 1024L, max_block = 4096L) {
  b <- block_s
  if (is.null(b)) b <- prefer
  b <- as.integer(b)
  if (!is.finite(b) || b <= 0) b <- prefer
  b <- max(1L, b)
  b <- min(max_block, b)
  b <- min(b, N)
  b
}

.gpsa_get_Z_block <- function(db, rows, cols) {
  if (!is.null(db$Z_mem)) {
    blk <- db$Z_mem[rows, cols, drop = FALSE]
  } else {
    blk <- rhdf5::h5read(db$h5file, "Z",
                         index = list(rows, cols),
                         drop = FALSE,
                         read.attributes = FALSE)
  }
  if (!is.matrix(blk)) {
    blk <- matrix(as.numeric(blk), nrow = length(rows), ncol = length(cols))
  } else {
    storage.mode(blk) <- "double"
  }
  blk
}

.gpsa_drop_datatable_class_nocopy <- function(x) {
  if (!is.null(class(x)) && any(class(x) == "data.table")) class(x) <- "data.frame"
  x
}

.gpsa_h5exists <- function(h5file, name) {
  ok <- tryCatch({
    ls0 <- rhdf5::h5ls(h5file, recursive = FALSE)
    name %in% ls0$name
  }, error = function(e) FALSE)
  
  if (ok) return(TRUE)
  
  tryCatch({
    rhdf5::h5read(h5file, name, index = list(1, 1), read.attributes = FALSE)
    TRUE
  }, error = function(e2) FALSE)
}

#' Collapse duplicate gene names by max absolute value
#'
#' @param stats Named numeric vector
#' @return Named numeric vector with unique names
#' @export
gpsa_collapse_duplicates_maxabs <- function(stats) {
  stats <- stats[!is.na(stats)]
  nm <- names(stats)
  if (is.null(nm) || length(nm) == 0) stop("stats must be a named numeric vector.")
  if (!anyDuplicated(nm)) return(stats)
  
  o <- order(abs(stats), decreasing = TRUE)
  s2 <- stats[o]
  keep <- !duplicated(names(s2))
  s2[keep]
}

#' Prepare a geneList for GSEA/search
#'
#' @param x Named numeric vector (gene symbols as names)
#' @param decreasing Sort decreasing?
#' @param collapse_duplicates Collapse duplicates by max|value|?
#' @param na.rm Remove NA values?
#' @return Sorted named numeric vector
#' @export
gpsa_prepare_geneList <- function(x, decreasing = TRUE, collapse_duplicates = TRUE, na.rm = TRUE) {
  nm <- names(x)
  x <- as.numeric(x)
  if (!is.null(nm) && length(nm) == length(x)) names(x) <- nm
  
  if (na.rm) x <- x[!is.na(x)]
  if (is.null(names(x)) || any(names(x) == "")) stop("geneList must have gene symbols as names.")
  if (collapse_duplicates) x <- gpsa_collapse_duplicates_maxabs(x)
  sort(x, decreasing = decreasing)
}

