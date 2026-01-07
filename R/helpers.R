#' Build Up/Down query from a signature geneList
#'
#' Extracts top N and bottom N genes from a geneList to create
#' an Up/Down bidirectional query.
#'
#' @param geneList Named numeric vector (e.g., logFC values)
#' @param n_each Number of genes for each direction (default: 500)
#'
#' @return List with Up and Down gene sets ready for gpsa_search
#' @export
#'
#' @examples
#' \dontrun{
#' db <- gpsa_load_db()
#' gl <- gpsa_get_signature(db, "D21455", dataset = "logFC")
#' query <- gpsa_make_updown_query(gl, n_each = 300)
#' out <- gpsa_search(db, query)
#' }
gpsa_make_updown_query <- function(geneList, n_each = 500) {
  geneList <- geneList[!is.na(geneList)]
  if (is.null(names(geneList))) stop("geneList must be named.")
  if (anyDuplicated(names(geneList))) geneList <- gpsa_collapse_duplicates_maxabs(geneList)
  if (length(geneList) < 2 * n_each) stop("geneList non-NA entries < 2*n_each")
  
  ord <- order(geneList, decreasing = TRUE)
  up_genes <- names(geneList)[ord][seq_len(n_each)]
  dn_genes <- names(geneList)[ord][(length(geneList) - n_each + 1):length(geneList)]
  
  list(
    Up = list(genes = up_genes, sign = +1, weight = 1),
    Down = list(genes = dn_genes, sign = -1, weight = 1)
  )
}

#' Adaptive Up/Down gene selection by energy
#'
#' Selects up and down genes adaptively based on energy capture threshold,
#' avoiding arbitrary fixed cutoffs like 500.
#'
#' @param geneList Named numeric vector
#' @param eta Energy capture threshold (default: 0.90)
#' @param min_n Minimum genes per direction
#' @param max_n Maximum genes per direction
#'
#' @return List with Up, Down gene sets and .info
#' @export
gpsa_make_updown_query_adaptive <- function(geneList, eta = 0.90, min_n = 200, max_n = 2000) {
  if (is.null(names(geneList))) stop("geneList must be named.")
  geneList <- geneList[!is.na(geneList)]
  if (anyDuplicated(names(geneList))) geneList <- gpsa_collapse_duplicates_maxabs(geneList)
  
  up <- geneList[geneList > 0]
  dn <- geneList[geneList < 0]
  
  if (length(up) < min_n || length(dn) < min_n) {
    stop("Not enough positive/negative genes for adaptive up/down (need at least min_n each).")
  }
  
  pick_by_energy <- function(v, eta, min_n, max_n) {
    v2 <- v^2
    ord <- order(v2, decreasing = TRUE)
    cum <- cumsum(v2[ord]) / sum(v2)
    k0 <- which(cum >= eta)[1]
    if (is.na(k0)) k0 <- length(v)
    
    k <- max(k0, min_n)
    k <- min(k, max_n, length(v))
    names(v)[ord][seq_len(k)]
  }
  
  up_genes <- pick_by_energy(up, eta, min_n, max_n)
  dn_genes <- pick_by_energy(abs(dn), eta, min_n, max_n)
  
  list(
    Up = list(genes = up_genes, sign = +1, weight = 1),
    Down = list(genes = dn_genes, sign = -1, weight = 1),
    .info = list(n_up = length(up_genes), n_down = length(dn_genes), eta = eta)
  )
}

#' Get rank of a signature in a result data.frame
#'
#' @param df Result data.frame with signature_id column
#' @param signature_id Signature ID to find
#'
#' @return Integer rank (row number) or NA
#' @export
gpsa_rank_in <- function(df, signature_id) {
  if (is.null(df) || nrow(df) == 0) return(NA_integer_)
  w <- which(df$signature_id == signature_id)
  if (length(w) == 0) return(NA_integer_)
  w[1]
}

#' Pick a signature by gene and cell line
#'
#' @param db A gpsa_db object
#' @param pbgene Perturbation gene name
#' @param cellline Cell line name (optional)
#' @param fallback Fallback signature ID if not found
#'
#' @return Signature ID or NULL
#' @export
gpsa_pick_signature <- function(db, pbgene, cellline = NULL, fallback = NULL) {
  m <- db$meta
  keep <- which(!is.na(m$pbgene) & m$pbgene == pbgene)
  if (!is.null(cellline)) keep <- keep[m$CellLineName[keep] == cellline]
  if (length(keep) > 0) return(m$signature_id[keep][1])
  if (!is.null(fallback) && fallback %in% db$sigs) return(fallback)
  NULL
}

