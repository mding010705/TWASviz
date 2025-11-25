#' Plot Pairwise Correlations and Overlaps of Multiple Gene Sets
#'
#' Plots heatmap for pairwise complete correlations for effect sizes of
#' multiple gene sets in lower triangle, number of pairwise complete effect
#' sizes for multiple gene sets in upper triangle. Diagonal displays total
#' number of effect sizes for a single gene set.
#'
#' @param betas A dataframe with rows of genes and columns of tissue/cell types.
#' Entry (i, j) is the TWAS effect size of gene i in tissue j, or NULL if gene i
#' is not expressed/significant in tissue j.
#' @param cor_method A string specifying the correlation method
#' (pearson, kendall, spearman).
#' @param tissue_names A vector of strings with length of the number of columns
#' in betas. These are the tissue names in the order that they appear in betas.
#'
#' @return Returns a heatmap plot of correlations and overlaps
#'
#' @examples
#' correlation_overlap_heatmap(betas = data.frame(x = c(1, 7, 3), y = c(1, -9, 3),
#'                             z = c(-1, 0, NA)), tissue_names = c("A", "B", "C"))
#'
#' @export
#' @import corrplot
#' @import psych


correlation_overlap_heatmap <- function(betas, cor_method = "pearson",
                                        tissue_names = 1:ncol(betas)){
  if (ncol(betas) == 0){
    return(NULL)
  }

  test_cor <- psych::corr.test(betas, method = cor_method, adjust = "fdr")
  corr_coef <- test_cor$r
  colnames(corr_coef) <- tissue_names
  row.names(corr_coef) <- tissue_names
  n_mat <- test_cor$n
  colnames(n_mat) <- tissue_names
  row.names(n_mat) <- tissue_names
  mat <- matrix(0, nrow = ncol(betas), ncol = ncol(betas))
  mat[upper.tri(mat, diag = TRUE)] <- test_cor$p.adj
  p_adj <- t(mat)
  p_adj <- rbind(0, p_adj)
  p_adj <- cbind(p_adj, 0)
  p_adj[is.na(p_adj)] <- 1
  upper_plot <- corrplot::corrplot(n_mat, method = "color", addCoef.col = "grey50",
                                   type='upper', number.cex=0.5,
                                   col = corrplot::COL1("YlGn"), is.corr = F,
                                   number.digits = 0, tl.pos = "lt",
                                   cl.ratio = 0.1, cl.pos = "r", na.label = "NA")
  return(tryCatch(corrplot::corrplot(corr_coef, p.mat = p_adj,
                            sig.level=c(0.01, 0.05, 0.1), pch.cex = 0.9,
                            type='lower', insig = 'label_sig', pch.col = 'green',
                            tl.pos = "n", is.cor = T, method = "square",
                            diag = FALSE, cl.ratio = 0.1, na.label = "NA",
                            add = T), error=function(e){}))

}

# [END]
