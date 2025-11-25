#' Volcano Plot of Gene Effect Sizes vs. p-values
#'
#' Produces volcano plot of gene effect sizes vs. -log10(p-values).
#' Labels significant genes.
#' Meant to take in TWAS results from PrediXcan which contain a gene, zscore
#' (standardized effect size), and a pvalue column, summarizing the individual
#' associations each gene has with your phenotype of interest.
#'
#' @param betas_p_vals Data frame with effect sizes and p-values.
#' @param p_thresh Numeric. P-value threshold for significance (default = 0.05).
#' @param vline Numeric vector. Vertical reference lines (default uses ±1.96).
#' @param add_sig_gene_labels Logical. Whether to label significant genes.
#'
#' @return Returns a ggplot volcano plot.
#'
#' @examples
#' volcano_plot(betas_p_vals = data.frame(zscore = c(1, 0, 3, 4, 0.5),
#'                                  pvalue = c(0.1, 0.9, 0.02, 0.11, 0.3)))
#'
#' @export
#' @import ggplot2
#' @import ggrepel

volcano_plot <- function(betas_p_vals, p_thresh = 0.05,
                         vline = c(qnorm(0.025), qnorm(0.975)),
                         add_sig_gene_labels = TRUE){
  # Input validation
  if (!is.data.frame(betas_p_vals)) {
    stop("betas_p_vals must be a data.frame.")
  }

  # Check for p-value column
  if (!"pvalue" %in% colnames(betas_p_vals)) {
    stop("betas_p_vals must contain a column named 'pvalue'.")
  }

  # Check for zscore column
  if (!"zscore" %in% colnames(betas_p_vals)) {
    stop("betas_p_vals must contain a column named 'zscore'.")
  }

  betas_p_vals$pvalue_significant <- betas_p_vals$pvalue < p_thresh
  if(!("gene" %in% colnames(betas_p_vals))){
    betas_p_vals$gene <- rownames(betas_p_vals)
  }

  # Base ggplot
  plt <- ggplot2::ggplot(
    betas_p_vals,
    ggplot2::aes(
      x = zscore,
      y = -log10(pvalue),
      color = pvalue_significant,
      label = gene)) +
    ggplot2::geom_point() +
    ggplot2::theme_minimal() +
    ggplot2::scale_color_manual(
      values = c("FALSE" = "black", "TRUE" = "blue"),
      name = paste0("p < ", p_thresh)
    ) +
    ggplot2::geom_vline(xintercept = vline, color = "orange") +
    ggplot2::geom_hline(yintercept = -log10(p_thresh), color = "red") +
    ggplot2::labs(
      title = "Volcano Plot",
      x = "Effect size (z-score)",
      y = expression(-log[10](pvalue)))

  # Add gene labels if add_sig_gene_labels TRUE
  if (add_sig_gene_labels == TRUE) {
    plt <- plt +
      ggrepel::geom_text_repel(
        data = betas_p_vals[betas_p_vals$pvalue_significant, ],
        max.overlaps = Inf)
  }

  return(plt)
}

# [END]
