#' Attach per-channel scores to result data.frame
#'
#' Adds raw and signed NES scores for each channel to the result table.
#'
#' @param res Result data.frame from gpsa_search
#' @param NES_raw Raw NES matrix (channels x signatures)
#' @param NES_signed Signed NES matrix
#' @param suffix_raw Suffix for raw score columns
#' @param suffix_signed Suffix for signed score columns
#'
#' @return Updated data.frame with channel columns
#' @export
gpsa_attach_channel_scores <- function(res, NES_raw, NES_signed,
                                       suffix_raw = "_raw", suffix_signed = "_signed") {
  if (is.null(NES_raw) || is.null(NES_signed)) stop("NES_raw / NES_signed missing.")
  sig_ids <- colnames(NES_raw)
  if (!identical(sig_ids, colnames(NES_signed))) stop("NES_raw and NES_signed column names mismatch.")
  idx <- match(res$signature_id, sig_ids)
  
  for (ch in rownames(NES_raw)) {
    res[[paste0(ch, suffix_raw)]] <- as.numeric(NES_raw[ch, idx])
    res[[paste0(ch, suffix_signed)]] <- as.numeric(NES_signed[ch, idx])
  }
  res
}

#' Compute dual-mode scores (Phenocopy vs Discovery)
#'
#' Computes Score_bi (strict phenocopy) and Score_uni (discovery) from
#' signed NES matrix.
#'
#' @param NES_signed Signed NES matrix (positive = desired direction)
#' @param eps Small value to avoid division by zero
#' @param p_mode P-value mode: "none" or "norm"
#' @param adjust_method P-value adjustment method
#'
#' @return List with Score_bi, Score_uni, Balance, Driver, and optional p/FDR
#' @export
gpsa_compute_mode_scores <- function(NES_signed,
                                     eps = 1e-9,
                                     p_mode = c("none", "norm"),
                                     adjust_method = "BH") {
  p_mode <- match.arg(p_mode)
  m <- nrow(NES_signed)
  N <- ncol(NES_signed)
  if (is.null(rownames(NES_signed))) rownames(NES_signed) <- paste0("C", seq_len(m))
  
  Score_bi <- apply(NES_signed, 2, min)
  Score_uni <- apply(NES_signed, 2, max)
  Driver_idx <- apply(NES_signed, 2, which.max)
  Driver <- rownames(NES_signed)[Driver_idx]
  Balance <- Score_bi / (Score_uni + eps)
  
  p_bi <- p_uni <- FDR_bi <- FDR_uni <- rep(NA_real_, N)
  
  if (p_mode == "norm") {
    p_ch <- pnorm(-NES_signed)
    p_bi <- apply(p_ch, 2, max)
    p_uni <- pmin(1, m * apply(p_ch, 2, min))
    FDR_bi <- p.adjust(p_bi, method = adjust_method)
    FDR_uni <- p.adjust(p_uni, method = adjust_method)
  }
  
  list(
    Score_bi = Score_bi,
    Score_uni = Score_uni,
    Balance = Balance,
    Driver = Driver,
    p_bi = p_bi,
    FDR_bi = FDR_bi,
    p_uni = p_uni,
    FDR_uni = FDR_uni
  )
}

#' Postprocess search results with dual-mode scoring
#'
#' Adds phenocopy (Score_bi) and discovery (Score_uni) scores to search results.
#' Creates two ranked tables: res_similarity and res_discovery.
#'
#' @param out Output from gpsa_search with return_components=TRUE
#' @param attach_channel_cols Attach per-channel score columns?
#' @param p_mode P-value computation mode
#' @param eps Epsilon for Balance calculation
#' @param filter_positive Filter to positive scores only?
#'
#' @return Updated output with res, res_similarity, res_discovery
#' @export
#'
#' @examples
#' \dontrun{
#' db <- gpsa_load_db(preload_Z = TRUE)
#' gl <- gpsa_get_signature(db, "D21455", dataset = "logFC")
#' query <- gpsa_make_updown_query(gl)
#' out <- gpsa_search(db, query, return_components = TRUE)
#' out <- gpsa_postprocess_modes(out, p_mode = "norm")
#' 
#' # Phenocopy hits (all channels positive)
#' head(out$res_similarity)
#' 
#' # Discovery hits (at least one channel strong)
#' head(out$res_discovery)
#' }
gpsa_postprocess_modes <- function(out,
                                   attach_channel_cols = TRUE,
                                   p_mode = c("none", "norm"),
                                   eps = 1e-9,
                                   filter_positive = TRUE) {
  p_mode <- match.arg(p_mode)
  if (is.null(out$NES_signed)) {
    warning("No NES_signed found. Run gpsa_search(..., return_components=TRUE) first.")
    return(out)
  }
  
  res <- out$res
  NES_raw <- out$NES_raw
  NES_signed <- out$NES_signed
  sig_ids <- colnames(NES_signed)
  idx <- match(res$signature_id, sig_ids)
  
  if (attach_channel_cols) {
    res <- gpsa_attach_channel_scores(res, NES_raw, NES_signed,
                                      suffix_raw = "_raw", suffix_signed = "_signed")
  }
  
  ms <- gpsa_compute_mode_scores(NES_signed, eps = eps, p_mode = p_mode, adjust_method = "BH")
  
  res$Score_bi <- as.numeric(ms$Score_bi[idx])
  res$Score_uni <- as.numeric(ms$Score_uni[idx])
  res$Balance <- as.numeric(ms$Balance[idx])
  res$Driver <- as.character(ms$Driver[idx])
  
  res$p_bi <- as.numeric(ms$p_bi[idx])
  res$FDR_bi <- as.numeric(ms$FDR_bi[idx])
  res$p_uni <- as.numeric(ms$p_uni[idx])
  res$FDR_uni <- as.numeric(ms$FDR_uni[idx])
  
  out$res <- res
  out$mode_info <- list(
    definition = list(
      Score_bi = "min over NES_signed channels (phenocopy strict)",
      Score_uni = "max over NES_signed channels (one-sided discovery)",
      Balance = "Score_bi / (Score_uni + eps)",
      Driver = "channel name that maximizes NES_signed"
    ),
    p_mode = p_mode
  )
  
  ## Create two entrance tables
  res_similarity <- res
  if (filter_positive) res_similarity <- res_similarity[res_similarity$Score_bi > 0, , drop = FALSE]
  res_similarity <- res_similarity[order(res_similarity$Score_bi, decreasing = TRUE), , drop = FALSE]
  
  res_discovery <- res
  if (filter_positive) res_discovery <- res_discovery[res_discovery$Score_uni > 0, , drop = FALSE]
  res_discovery <- res_discovery[order(res_discovery$Score_uni, decreasing = TRUE), , drop = FALSE]
  
  out$res_similarity <- res_similarity
  out$res_discovery <- res_discovery
  
  out
}

