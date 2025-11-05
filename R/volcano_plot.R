#' Volcano Plot of Gene Effect Sizes vs. p-values
#'
#' Produces volcano plot of gene effect sizes vs. -log10(p-values).
#' Labels significant genes.
#'
#' @param betas_p_vals Dataframe with columns gene_name (to label points with),
#' betas (containing gene effect sizes), pvalue (significance of effect sizes)
#' @param p_thresh Max value of p value to be considered significant.
#' @param vline Scalar or vector of values at which to draw vertical lines
#' on volcano plot.
#' @param add_sig_gene_labels Boolean. If TRUE, significant points will be
#' labeled with their gene names. If FALSE, no points will be labeled.
#'
#' @return Returns a ggplot volcano plot.
#'
#' @examples
#' volcano_plot(betas_p_vals = data.frame(beta = c(1, 0, 3, 4, 0.5),
#'                                  pvalue = c(0.1, 0.9, 0.02, 0.11, 0.3)))
#'
#' @export
#' @import ggplot2
#' @import ggrepel

volcano_plot <- function(betas_p_vals, p_thresh = 0.05,
                         vline = c(qnorm(0.025), qnorm(0.975)),
                         add_sig_gene_labels = TRUE){
  betas_p_vals$pvalue_significant <- betas_p_vals$pvalue < p_thresh
  if(!("gene" %in% colnames(betas_p_vals))){
    betas_p_vals$gene <- rownames(betas_p_vals)
  }
  if(add_sig_gene_labels == TRUE){
    return(ggplot2::ggplot(betas_p_vals,
                           ggplot2::aes(x=zscore, y=-log10(pvalue),
                                        col=pvalue_significant, label=gene)) +
             ggplot2::geom_point() +
             ggplot2::theme_minimal() +
             ggrepel::geom_text_repel() +
             ggplot2::scale_color_manual(values = c("blue", "black")) +
             ggplot2::geom_vline(xintercept = vline, col="orange") +
             ggplot2::geom_hline(yintercept = -log10(p_thresh), col="red"))
  } else{
    return(return(ggplot2::ggplot(betas_p_vals,
                                  ggplot2::aes(x=zscore, y=-log10(pvalue),
                                               col=pvalue_significant, label=gene)) +
                    ggplot2::geom_point() +
                    ggplot2::theme_minimal() +
                    ggplot2::scale_color_manual(values = c("blue", "black")) +
                    ggplot2::geom_vline(xintercept = vline, col="orange") +
                    ggplot2::geom_hline(yintercept = -log10(p_thresh), col="red")))
  }


}

# [END]
