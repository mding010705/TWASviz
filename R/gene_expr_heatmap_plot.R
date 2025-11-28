#' Reorder Correlation Matrix using Complete Hierarchical Clustering
#'
#' Produces an a matrix reordered based on complete hierarchical clustering.
#'
#' @param corr_mat Correlation matrix.
#'
#' @return Returns a reordered correlation matrix.

reorder_corr_hclust <- function(corr_mat){
  dist_mat <- as.dist((1-corr_mat)/2)
  hclust_mat <- hclust(dist_mat)
  corr_mat <- corr_mat[hclust_mat$order, hclust_mat$order]
  return(corr_mat)
}

#' Down Sample Matrix
#'
#' Produces an a matrix of smaller dimensions from the input matrix, but with
#' similar patterns. Adapted from
#' https://gdevailly.netlify.app/post/plotting-big-matrices-in-r/ .
#'
#' @param mat Gene expression correlation matrix.
#' @param cl Parallel cluster.
#' @param target_height Max number of rows that the output matrix should have.
#' @param target_width Max number of columns that the output matrix should have.
#' @param summary_func Function with which to down sample the matrix.
#'
#' @return Returns a heatmap image of gene expression correlations.
#'
#' @import parallel
#'
#' @references
#' Devailly, G. (2021, October 14). Plotting heatmaps from big matrices
#' in R. <https://gdevailly.netlify.app/post/plotting-big-matrices-in-r/>.
#'
#' R Core Team (2024). R: A Language and Environment for Statistical Computing.
#' R Foundation for Statistical Computing, Vienna, Austria. https://www.R-project.org/.

redim_matrix <- function(
    mat, cl = NULL,
    target_height = 300,
    target_width = 300,
    summary_func = function(x) mean(x, na.rm = TRUE)) {

  if(target_height > nrow(mat) | target_width > ncol(mat)) {
    return(reorder_corr_hclust(mat))
  }
  mat <- reorder_corr_hclust(mat)
  seq_height <- round(seq(1, nrow(mat), length.out = target_height + 1))
  seq_width  <- round(seq(1, ncol(mat), length.out = target_width  + 1))
  # complicated way to write a double for loop
  parallel::clusterExport(cl, varlist = c("target_width", "target_height",
                                          "summary_func", "mat", "seq_height",
                                          "seq_width"),
                          envir = environment())


  return(do.call(rbind, parallel::parLapply(cl, seq_len(target_height), function(i) { # i is row
    vapply(seq_len(target_width), function(j) { # j is column
      summary_func(
        mat[
          seq(seq_height[i], seq_height[i + 1]),
          seq(seq_width[j] , seq_width[j + 1] )
        ]
      )
    }, 0.0)
  })))
}


#' Heatmap of all Pairwise Correlations of Expression of Many Genes
#'
#' Produces an image of a heatmap for all of the pairwise correlations of
#' many gene expressions. Intended for use for expression of > 500 genes.
#'
#' @param gene_expr Gene expression matrix. Rows are individual samples.
#' Columns are genes.
#' @param target_height Max number of rows that the heatmap should have.
#' @param target_width Max number of columns that the heatmap should have.
#' @param summary_func Function with which to down sample a large
#' correlation matrix.
#' @param n_core Number of cores for parallel processing
#' @param col_pal Vector of hex colors, in order of corresponding colours
#' from -1 to 1.
#' @param title Title for the heatmap.
#'
#' @return Returns a heatmap image of gene correlations.
#'
#' @examples
#' \dontrun{
#'   corr_heatmap(gene_expr = matrix(rnorm(500^2), nrow = 500))
#' }
#'
#' @export
#' @import grDevices
#' @import parallel
#'
#' @references R Core Team (2024). R: A Language and Environment for Statistical
#' Computing. R Foundation for Statistical Computing, Vienna, Austria.
#' https://www.R-project.org/.

corr_heatmap <- function(gene_expr, target_height = 300,
                         target_width = 300,
                summary_func = function(x) x[which.max(abs(x))],
                n_core = parallel::detectCores() - 1,
                col_pal =
                  grDevices::colorRampPalette(c("red", "white", "blue"))(50),
                title = "Gene expression correlation matrix"){
  # Input validation
  if (!is.matrix(gene_expr)) stop("gene_expr must be a numeric matrix.")
  if (nrow(gene_expr) < 2 || ncol(gene_expr) < 2)
    stop("gene_expr must have at least 2 rows and 2 columns.")

  # Compute correlation and downsample
  cl <- parallel::makeCluster(n_core, type="PSOCK")
  corr_mat <- redim_matrix(cor(gene_expr), cl = cl,
                           summary_func = summary_func,
                           target_height = target_height,
                           target_width = target_width)
  parallel::stopCluster(cl)

  # Plot heatmap
  return(image(
    corr_mat[,seq(ncol(corr_mat), 1)],
    axes = FALSE,
    col = col_pal,
    breaks = seq(-1, 1, length.out = length(col_pal) + 1),
    main = title, asp = 1
  ))
}

# [END]
