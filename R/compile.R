#' Compile gene-set query
#'
#' Compiles one or more gene sets into query weights for search.
#'
#' @param db A gpsa_db object
#' @param query_sets Named list of gene sets, each with $genes (character vector),
#'   $sign (+1 or -1), and optional $weight (default 1)
#' @param min_size Minimum genes after intersection with DB (default: 20)
#'
#' @return Compiled query object
#' @export
#'
#' @examples
#' \dontrun{
#' db <- gpsa_load_db()
#' query <- list(
#'   Up = list(genes = up_genes, sign = +1),
#'   Down = list(genes = down_genes, sign = -1)
#' )
#' cq <- gpsa_compile_query(db, query)
#' }
gpsa_compile_query <- function(db, query_sets, min_size = 20) {
  if (!is.list(query_sets) || length(query_sets) < 1) stop("query_sets must be a non-empty list.")
  if (is.null(names(query_sets)) || any(names(query_sets) == "")) {
    names(query_sets) <- paste0("S", seq_along(query_sets))
  }
  
  set_names <- names(query_sets)
  idx_list <- vector("list", length(query_sets))
  k <- integer(length(query_sets))
  sign <- numeric(length(query_sets))
  weight <- numeric(length(query_sets))
  
  for (j in seq_along(query_sets)) {
    s <- query_sets[[j]]
    if (is.null(s$genes) || is.null(s$sign)) stop("Each set needs $genes and $sign.")
    wj <- if (is.null(s$weight)) 1 else as.numeric(s$weight)
    sj <- ifelse(as.numeric(s$sign) >= 0, 1, -1)
    
    idx <- unname(db$gene2idx[as.character(s$genes)])
    idx <- idx[!is.na(idx)]
    idx <- unique(idx)
    
    if (length(idx) < min_size) {
      stop("Gene set too small after intersection: ", set_names[j], " (", length(idx), ")")
    }
    
    idx_list[[j]] <- idx
    k[j] <- length(idx)
    sign[j] <- sj
    weight[j] <- wj
  }
  
  union_idx <- sort(unique(unlist(idx_list, use.names = FALSE)))
  
  ## Compiled weights for total score
  w_union <- numeric(length(union_idx))
  names(w_union) <- as.character(union_idx)
  for (j in seq_along(idx_list)) {
    coef <- weight[j] * sign[j] / sqrt(k[j])
    w_union[as.character(idx_list[[j]])] <- w_union[as.character(idx_list[[j]])] + coef
  }
  w_union <- as.numeric(w_union)
  
  list(
    union_idx = union_idx,
    w_union = w_union,
    idx_list = idx_list,
    k = k,
    sign = sign,
    weight = weight,
    set_names = set_names
  )
}

#' Compile full geneList to vector query weights
#'
#' Converts a full geneList (named numeric vector) into query weights
#' using rankZ transformation.
#'
#' @param db A gpsa_db object
#' @param geneList Named numeric vector (gene symbols as names)
#' @param transform Transformation: "rankZ" or "zscore"
#' @param ties_method Ties method for ranking
#' @param normalize L2-normalize weights?
#' @param min_size Minimum overlap with DB
#'
#' @return Compiled vector query object
#' @export
gpsa_compile_signature_vector <- function(
    db,
    geneList,
    transform = c("rankZ", "zscore"),
    ties_method = "average",
    normalize = TRUE,
    min_size = 200
) {
  transform <- match.arg(transform)
  
  if (is.null(names(geneList)) || any(names(geneList) == "")) {
    stop("geneList must be a named numeric vector (names=gene symbols).")
  }
  geneList <- geneList[!is.na(geneList)]
  
  if (anyDuplicated(names(geneList))) geneList <- gpsa_collapse_duplicates_maxabs(geneList)
  
  idx <- unname(db$gene2idx[names(geneList)])
  keep <- !is.na(idx)
  idx <- idx[keep]
  v <- as.numeric(geneList[keep])
  
  M <- length(v)
  if (M < min_size) stop("vector query overlap too small: ", M, " (<", min_size, ")")
  
  if (transform == "rankZ") {
    r <- rank(v, ties.method = ties_method)
    p <- (r - 0.5) / M
    w <- qnorm(p)
    w[!is.finite(w)] <- 0
  } else {
    mu <- mean(v)
    sdv <- sd(v)
    if (!is.finite(sdv) || sdv <= 0) sdv <- 1
    w <- (v - mu) / sdv
    w[!is.finite(w)] <- 0
  }
  
  if (normalize) {
    s <- sqrt(sum(w^2))
    if (is.finite(s) && s > 0) w <- w / s
  }
  
  o <- order(idx)
  list(
    union_idx = idx[o],
    w_union = w[o],
    M = M,
    transform = transform,
    normalize = normalize
  )
}

#' Adaptive sparsification by energy capture
#'
#' Reduces query size while preserving a specified fraction of signal energy.
#'
#' @param union_idx Gene indices
#' @param w_union Weights
#' @param eta Energy capture threshold (default: 0.90)
#' @param min_k Minimum genes to keep
#' @param max_k Maximum genes to keep
#'
#' @return Sparsified query object
#' @export
gpsa_sparsify_query_energy <- function(
    union_idx,
    w_union,
    eta = 0.90,
    min_k = 200,
    max_k = NULL
) {
  if (!is.finite(eta) || eta <= 0 || eta > 1) stop("eta must be in (0,1].")
  w2 <- w_union^2
  tot <- sum(w2)
  if (!is.finite(tot) || tot <= 0) stop("Invalid weights: sum(w^2) <= 0.")
  
  ord <- order(w2, decreasing = TRUE)
  cum <- cumsum(w2[ord]) / tot
  k0 <- which(cum >= eta)[1]
  if (is.na(k0)) k0 <- length(w_union)
  
  k <- max(k0, min_k)
  if (!is.null(max_k)) k <- min(k, as.integer(max_k))
  k <- min(k, length(w_union))
  
  pick <- ord[seq_len(k)]
  capture <- sum(w2[pick]) / tot
  
  o <- order(union_idx[pick])
  list(
    union_idx = union_idx[pick][o],
    w_union = w_union[pick][o],
    k = k,
    eta = eta,
    capture = capture,
    total_k = length(w_union)
  )
}

