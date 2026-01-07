#' Compute score statistics
#'
#' @param score_vec Numeric vector of scores
#' @param robust Use robust statistics (median/MAD)?
#' @return List with tau, p_rank, z_ref, p_ref
#' @export
gpsa_score_stats <- function(score_vec, robust = TRUE) {
  n <- length(score_vec)
  
  r <- rank(score_vec, ties.method = "average")
  pct <- (r - 0.5) / n
  tau <- 100 * (2 * pct - 1)
  
  p_rank <- 2 * pmin(pct, 1 - pct)
  p_rank[p_rank < (0.5 / n)] <- (0.5 / n)
  
  if (robust) {
    mu <- median(score_vec)
    s <- mad(score_vec, constant = 1.4826)
    if (!is.finite(s) || s <= 0) s <- sd(score_vec)
  } else {
    mu <- mean(score_vec)
    s <- sd(score_vec)
  }
  if (!is.finite(s) || s <= 0) s <- 1
  
  z_ref <- (score_vec - mu) / s
  p_ref <- 2 * pnorm(-abs(z_ref))
  
  list(tau = tau, p_rank = p_rank, z_ref = z_ref, p_ref = p_ref)
}

#' Search with gene-set query
#'
#' Performs similarity search using one or more gene sets.
#'
#' @param db A gpsa_db object
#' @param query_sets Named list of gene sets (see gpsa_compile_query)
#' @param return_components Return per-channel NES scores?
#' @param block_s Block size for reading (NULL = auto)
#' @param robust Use robust statistics?
#'
#' @return List with res (data.frame), compile info, and optionally NES_raw/NES_signed
#' @export
#'
#' @examples
#' \dontrun{
#' db <- gpsa_load_db(preload_Z = TRUE)
#' 
#' # Single gene set query
#' query <- list(MySet = list(genes = my_genes, sign = +1))
#' out <- gpsa_search(db, query)
#' head(out$res)
#' 
#' # Up/Down query
#' query <- list(
#'   Up = list(genes = up_genes, sign = +1),
#'   Down = list(genes = down_genes, sign = -1)
#' )
#' out <- gpsa_search(db, query, return_components = TRUE)
#' }
gpsa_search <- function(
    db,
    query_sets,
    return_components = TRUE,
    block_s = NULL,
    robust = TRUE
) {
  cq <- gpsa_compile_query(db, query_sets)
  
  N <- db$N
  k_union <- length(cq$union_idx)
  Score <- numeric(N)
  
  if (return_components) {
    m <- length(cq$idx_list)
    Wch <- matrix(0, nrow = k_union, ncol = m)
    colnames(Wch) <- cq$set_names
    
    for (i in seq_len(m)) {
      idx <- cq$idx_list[[i]]
      pos <- match(idx, cq$union_idx)
      Wch[pos, i] <- cq$weight[i] / sqrt(length(idx))
    }
    
    NES_raw <- matrix(NA_real_, nrow = m, ncol = N)
    rownames(NES_raw) <- cq$set_names
    colnames(NES_raw) <- db$sigs
  }
  
  bsz <- .gpsa_pick_block_s(block_s, N, prefer = 1024L, max_block = 4096L)
  
  t0 <- system.time({
    for (j in seq.int(1L, N, by = bsz)) {
      j_end <- min(N, j + bsz - 1L)
      
      blk <- .gpsa_get_Z_block(db, cq$union_idx, j:j_end)
      
      Score[j:j_end] <- as.numeric(crossprod(cq$w_union, blk))
      
      if (return_components) {
        NES_raw[, j:j_end] <- crossprod(Wch, blk)
      }
    }
  })
  
  st <- gpsa_score_stats(Score, robust = robust)
  FDR_rank <- p.adjust(st$p_rank, "BH")
  FDR_ref <- p.adjust(st$p_ref, "BH")
  
  res <- data.frame(
    signature_id = db$sigs,
    Score = Score,
    tau = st$tau,
    p_rank = st$p_rank,
    FDR_rank = FDR_rank,
    z_ref = st$z_ref,
    p_ref = st$p_ref,
    FDR_ref = FDR_ref,
    stringsAsFactors = FALSE
  )
  
  res <- merge(res, db$meta, by = "signature_id", all.x = TRUE, sort = FALSE)
  res <- res[order(res$Score, decreasing = TRUE), ]
  
  out <- list(res = res, compile = cq)
  attr(out, "timing") <- list(search_sec = unname(t0["elapsed"]))
  
  if (return_components) {
    NES_signed <- NES_raw
    for (i in seq_len(nrow(NES_signed))) {
      NES_signed[i, ] <- NES_signed[i, ] * cq$sign[i]
    }
    out$NES_raw <- NES_raw
    out$NES_signed <- NES_signed
  }
  
  out
}

#' Search with vector query (full geneList)
#'
#' Performs similarity search using a full gene expression profile.
#'
#' @param db A gpsa_db object
#' @param geneList Named numeric vector (gene symbols as names)
#' @param eta Energy capture threshold for sparsification (NULL = no sparsification)
#' @param min_k Minimum genes to use
#' @param max_k Maximum genes to use
#' @param transform Transformation: "rankZ" or "zscore"
#' @param ties_method Ties method for ranking
#' @param normalize L2-normalize weights?
#' @param block_s Block size for reading
#' @param robust Use robust statistics?
#'
#' @return List with res (data.frame) and compile info
#' @export
gpsa_search_vector <- function(
    db,
    geneList,
    eta = 0.90,
    min_k = 300,
    max_k = 5000,
    transform = c("rankZ", "zscore"),
    ties_method = "average",
    normalize = TRUE,
    block_s = NULL,
    robust = TRUE
) {
  transform <- match.arg(transform)
  
  cq0 <- gpsa_compile_signature_vector(
    db = db,
    geneList = geneList,
    transform = transform,
    ties_method = ties_method,
    normalize = normalize,
    min_size = min_k
  )
  
  if (!is.null(eta)) {
    cq <- gpsa_sparsify_query_energy(
      union_idx = cq0$union_idx,
      w_union = cq0$w_union,
      eta = eta,
      min_k = min_k,
      max_k = max_k
    )
  } else {
    cq <- list(
      union_idx = cq0$union_idx,
      w_union = cq0$w_union,
      k = length(cq0$w_union),
      eta = NA_real_,
      capture = 1,
      total_k = length(cq0$w_union)
    )
  }
  
  N <- db$N
  Score <- numeric(N)
  
  bsz <- .gpsa_pick_block_s(block_s, N, prefer = 1024L, max_block = 4096L)
  
  t0 <- system.time({
    for (j in seq.int(1L, N, by = bsz)) {
      j_end <- min(N, j + bsz - 1L)
      
      blk <- .gpsa_get_Z_block(db, cq$union_idx, j:j_end)
      Score[j:j_end] <- as.numeric(crossprod(cq$w_union, blk))
    }
  })
  
  st <- gpsa_score_stats(Score, robust = robust)
  FDR_rank <- p.adjust(st$p_rank, "BH")
  FDR_ref <- p.adjust(st$p_ref, "BH")
  
  res <- data.frame(
    signature_id = db$sigs,
    Score = Score,
    tau = st$tau,
    p_rank = st$p_rank,
    FDR_rank = FDR_rank,
    z_ref = st$z_ref,
    p_ref = st$p_ref,
    FDR_ref = FDR_ref,
    stringsAsFactors = FALSE
  )
  
  res <- merge(res, db$meta, by = "signature_id", all.x = TRUE, sort = FALSE)
  res <- res[order(res$Score, decreasing = TRUE), ]
  
  info <- list(
    type = "vector",
    transform = cq0$transform,
    normalize = cq0$normalize,
    overlap_M = cq0$M,
    k_used = cq$k,
    total_k = cq$total_k,
    eta = cq$eta,
    energy_capture = cq$capture
  )
  
  out <- list(res = res, compile = info)
  attr(out, "timing") <- list(search_sec = unname(t0["elapsed"]))
  out
}

#' Batch search (multiple queries at once)
#'
#' Performs multiple queries efficiently by sharing I/O.
#'
#' @param db A gpsa_db object
#' @param query_list Named list of queries (each is a query_sets list)
#' @param block_s Block size for reading
#' @param robust Use robust statistics?
#' @param compute_tau Compute tau scores?
#'
#' @return List with Score matrix and statistics matrices
#' @export
#'
#' @examples
#' \dontrun{
#' db <- gpsa_load_db(preload_Z = TRUE)
#' 
#' queries <- list(
#'   Q1 = list(Up = list(genes = genes1, sign = +1)),
#'   Q2 = list(Up = list(genes = genes2, sign = +1))
#' )
#' bat <- gpsa_batch_search(db, queries)
#' }
gpsa_batch_search <- function(
    db,
    query_list,
    block_s = NULL,
    robust = FALSE,
    compute_tau = TRUE
) {
  if (!is.list(query_list) || length(query_list) < 1) stop("query_list must be non-empty list.")
  if (is.null(names(query_list)) || any(names(query_list) == "")) {
    names(query_list) <- paste0("Q", seq_along(query_list))
  }
  
  Q <- length(query_list)
  q_names <- names(query_list)
  
  comp <- lapply(query_list, function(q) gpsa_compile_query(db, q))
  union_all <- sort(unique(unlist(lapply(comp, `[[`, "union_idx"), use.names = FALSE)))
  U <- length(union_all)
  N <- db$N
  
  ii <- integer(0); jj <- integer(0); xx <- numeric(0)
  for (q in seq_len(Q)) {
    cq <- comp[[q]]
    pos <- match(cq$union_idx, union_all)
    ii <- c(ii, pos)
    jj <- c(jj, rep(q, length(pos)))
    xx <- c(xx, cq$w_union)
  }
  W <- Matrix::sparseMatrix(i = ii, j = jj, x = xx, dims = c(U, Q))
  W <- as.matrix(W)
  colnames(W) <- q_names
  
  Score_mat <- matrix(NA_real_, nrow = N, ncol = Q)
  rownames(Score_mat) <- db$sigs
  colnames(Score_mat) <- q_names
  
  bsz <- .gpsa_pick_block_s(block_s, N, prefer = 1024L, max_block = 4096L)
  
  t0 <- system.time({
    for (j in seq.int(1L, N, by = bsz)) {
      j_end <- min(N, j + bsz - 1L)
      
      blk <- .gpsa_get_Z_block(db, union_all, j:j_end)
      
      Score_mat[j:j_end, ] <- as.matrix(crossprod(blk, W))
    }
  })
  
  z_ref <- Score_mat
  p_ref <- Score_mat
  tau <- if (compute_tau) Score_mat else NULL
  p_rank <- if (compute_tau) Score_mat else NULL
  
  for (q in seq_len(Q)) {
    sc <- Score_mat[, q]
    st <- gpsa_score_stats(sc, robust = robust)
    z_ref[, q] <- st$z_ref
    p_ref[, q] <- st$p_ref
    if (compute_tau) {
      tau[, q] <- st$tau
      p_rank[, q] <- st$p_rank
    }
  }
  
  FDR_ref <- apply(p_ref, 2, p.adjust, method = "BH")
  FDR_rank <- if (compute_tau) apply(p_rank, 2, p.adjust, method = "BH") else NULL
  
  list(
    Score = Score_mat,
    z_ref = z_ref,
    p_ref = p_ref,
    FDR_ref = FDR_ref,
    tau = tau,
    p_rank = p_rank,
    FDR_rank = FDR_rank,
    union_all = union_all,
    query_info = comp,
    timing = list(elapsed_sec = unname(t0["elapsed"]), union_size = U, block_s = block_s)
  )
}

#' Meta-analysis across queries (Stouffer method)
#'
#' Combines z-scores from multiple queries using Stouffer's method.
#'
#' @param z_mat Matrix of z-scores (signatures x queries)
#' @param weights Optional weights for each query
#'
#' @return data.frame with meta z-scores, p-values, and FDR
#' @export
gpsa_meta_stouffer <- function(z_mat, weights = NULL) {
  z_mat <- as.matrix(z_mat)
  Q <- ncol(z_mat)
  if (is.null(weights)) weights <- rep(1, Q)
  weights <- as.numeric(weights)
  if (length(weights) != Q) stop("weights length must match ncol(z_mat).")
  
  denom <- sqrt(sum(weights^2))
  z_meta <- as.numeric(z_mat %*% weights) / denom
  p_meta <- 2 * pnorm(-abs(z_meta))
  FDR_meta <- p.adjust(p_meta, "BH")
  
  data.frame(
    signature_id = rownames(z_mat),
    z_meta = z_meta,
    p_meta = p_meta,
    FDR_meta = FDR_meta,
    stringsAsFactors = FALSE
  )[order(z_meta, decreasing = TRUE), ]
}

