#' Compute GSEA running enrichment score
#'
#' @param geneList Named numeric vector (sorted descending)
#' @param geneSet Character vector of gene symbols
#' @param gsea_param GSEA weight parameter (default: 1)
#'
#' @return List with ES, peak, running scores, hit positions
#' @export
gpsa_gsea_running <- function(geneList, geneSet, gsea_param = 1) {
  geneList <- gpsa_prepare_geneList(geneList, decreasing = TRUE, collapse_duplicates = TRUE, na.rm = TRUE)
  geneSet <- unique(as.character(geneSet))
  overlap <- intersect(geneSet, names(geneList))
  if (length(overlap) < 10) stop("GeneSet overlap too small (<10).")
  
  hits <- names(geneList) %in% overlap
  N <- length(geneList)
  Nh <- sum(hits)
  
  w <- abs(geneList)^gsea_param
  norm_hit <- sum(w[hits])
  norm_miss <- N - Nh
  
  step <- numeric(N)
  step[hits] <- w[hits] / norm_hit
  step[!hits] <- -1 / norm_miss
  
  running <- cumsum(step)
  ES_pos <- max(running)
  ES_neg <- min(running)
  if (abs(ES_pos) >= abs(ES_neg)) {
    ES <- ES_pos
    peak <- which.max(running)
  } else {
    ES <- ES_neg
    peak <- which.min(running)
  }
  
  list(
    ES = ES,
    peak = peak,
    running = running,
    hit_positions = which(hits),
    overlap = overlap,
    geneList = geneList
  )
}

#' Plot GSEA running enrichment score
#'
#' @param geneList Named numeric vector
#' @param geneSet Character vector of gene symbols
#' @param main Plot title
#' @param method "base" (default) or "fgsea" if available
#' @param gsea_param GSEA weight parameter
#'
#' @return Invisibly returns the running score object (base) or ggplot (fgsea)
#' @export
#'
#' @examples
#' \dontrun{
#' db <- gpsa_load_db()
#' gl <- gpsa_get_signature(db, "D21455", dataset = "logFC")
#' gpsa_plot_gsea(gl, my_gene_set, main = "My Gene Set")
#' }
gpsa_plot_gsea <- function(geneList, geneSet, main = NULL, method = c("base", "fgsea"), gsea_param = 1) {
  method <- match.arg(method)
  
  if (method == "fgsea" && requireNamespace("fgsea", quietly = TRUE)) {
    geneList <- gpsa_prepare_geneList(geneList, decreasing = TRUE, collapse_duplicates = TRUE, na.rm = TRUE)
    geneSet <- unique(as.character(geneSet))
    overlap <- intersect(geneSet, names(geneList))
    if (length(overlap) < 10) stop("GeneSet overlap too small (<10).")
    
    p <- fgsea::plotEnrichment(overlap, geneList)
    if (!is.null(main) && requireNamespace("ggplot2", quietly = TRUE)) {
      p <- p + ggplot2::ggtitle(main)
    }
    print(p)
    return(invisible(p))
  }
  
  out <- gpsa_gsea_running(geneList, geneSet, gsea_param = gsea_param)
  if (is.null(main)) main <- "GSEA running enrichment score"
  
  plot(out$running, type = "l",
       xlab = "Rank in geneList (descending)",
       ylab = "Running ES",
       main = main)
  abline(h = 0, lty = 2, col = "grey60")
  abline(v = out$peak, lty = 3, col = "grey60")
  rug(out$hit_positions, ticksize = 0.02)
  mtext(sprintf("ES = %.3f | overlap = %d | peak@%d", out$ES, length(out$overlap), out$peak),
        side = 3, line = 0.2, cex = 0.85)
  
  invisible(out)
}

#' GSEA-style validation with fgsea
#'
#' Runs fgsea on up and down gene sets for validation.
#'
#' @param geneList Named numeric vector
#' @param up_genes Up-regulated gene set
#' @param down_genes Down-regulated gene set
#' @param nperm Number of permutations
#' @param minSize Minimum gene set size
#' @param maxSize Maximum gene set size
#'
#' @return fgsea result data.frame
#' @export
gpsa_fgsea_updown <- function(geneList, up_genes, down_genes,
                              nperm = 2000,
                              minSize = 10,
                              maxSize = 5000) {
  if (!requireNamespace("fgsea", quietly = TRUE)) stop("fgsea not installed.")
  gl <- gpsa_prepare_geneList(geneList, decreasing = TRUE, collapse_duplicates = TRUE, na.rm = TRUE)
  pathways <- list(upGenes = unique(as.character(up_genes)),
                   downGenes = unique(as.character(down_genes)))
  fg <- fgsea::fgsea(pathways = pathways, stats = gl,
                     nperm = nperm, minSize = minSize, maxSize = maxSize)
  fg[order(fg$pathway), , drop = FALSE]
}

