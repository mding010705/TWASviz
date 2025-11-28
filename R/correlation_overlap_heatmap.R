#' Plot Pairwise Correlations and Overlaps of Multiple Gene Sets
#'
#' Produces a two-layer heatmap summarizing:
#'
#' Lower triangle:
#'  Pairwise complete correlations between gene-level effect sizes
#'  Significance stars computed from FDR-adjusted p-values
#'
#' Upper triangle:
#'  Number of shared non-missing effect sizes (pairwise N)
#'  Text labels show N
#'
#' Diagonal:
#'  Number of non-missing values per tissue
#'
#' Internally uses `psych::corr.test()` to compute correlations, p-values,
#' and pairwise complete sample sizes.
#'
#' @references
#' Revelle, W. (2024). *psych: Procedures for Psychological, Psychometric,
#' and Personality Research.* Northwestern University.
#' Wickham, H. (2016). *ggplot2: Elegant Graphics for Data Analysis.* Springer.
#'
#' @param betas A numeric data frame or matrix where rows are genes and
#'   columns are tissues/cell types. Missing values (NA) indicate gene absence.
#' @param cor_method Correlation method passed to `psych::corr.test()`.
#'   One of `"pearson"`, `"kendall"`, `"spearman"`.
#' @param tissue_names Optional vector of tissue names. Must match
#'   `ncol(betas)`. Defaults to column names or sequential numbering.
#' @param main_title Optional main title for plot.
#'
#' @return A ggplot2 heatmap object, or the number of non-NA entries in `betas`
#' if it only contains one column.
#'
#' @examples
#' correlation_overlap_heatmap(
#'   betas = data.frame(
#'     A = c(1, 7, 3),
#'     B = c(1, 6, 3),
#'     C = c(10, -9, NA)
#'   ),
#'   tissue_names = c("A", "B", "C"))
#'
#' @export
#' @importFrom psych corr.test
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @import ggnewscale
#' @import magrittr

correlation_overlap_heatmap <- function(
    betas,
    tissue_names = colnames(betas),
    main_title = "Correlation (lower) and number of shared features (upper)",
    cor_method = "pearson") {

  # Input validation
  if (!is.numeric(as.matrix(betas))) {
    stop("`betas` must contain numeric values.")
  }

  if (is.vector(betas) || ncol(betas) == 1){
    message("Only one column was detected in `betas`. Returning the number of
            non-NA entries.")
    return(sum(!is.na(betas)))
  }

  if (!is.data.frame(betas) && !is.matrix(betas)) {
    stop("`betas` must be a data frame or matrix.")
  }

  if (ncol(betas) == 0) {
    stop("`betas` must have at least one column.")
  }

  valid_methods <- c("pearson", "spearman", "kendall")
  if (!(cor_method %in% valid_methods)) {
    stop("`cor_method` must be one of: 'pearson', 'spearman', 'kendall'.")
  }

  if (length(tissue_names) != ncol(betas)) {
    stop("`tissue_names` must have length equal to ncol(betas).")
  }

  # betas: samples × tissues matrix
  test_cor <- psych::corr.test(betas, method = cor_method, adjust = "fdr")

  corr_coef <- test_cor$r
  # Fill matrix with NA then insert p-values only in lower triangle
  p_adj <- matrix(1, nrow = ncol(betas), ncol = ncol(betas))
  p_adj[lower.tri(p_adj, diag = FALSE)] <- test_cor$p.adj
  diag(p_adj) <- 0
  p_adj[is.na(p_adj)] <- 1
  n_mat <- test_cor$n

  # when no NA, test_cor$n is an integer, not matrix
  if (!is.matrix(n_mat)){
    n_mat <- matrix(rep(n_mat, nrow(corr_coef) * ncol(corr_coef)),
                    nrow = nrow(corr_coef), ncol = ncol(corr_coef))
  }

  colnames(corr_coef) <- tissue_names
  rownames(corr_coef) <- tissue_names
  colnames(p_adj) <- tissue_names
  rownames(p_adj) <- tissue_names
  colnames(n_mat) <- tissue_names
  rownames(n_mat) <- tissue_names

  # Significance labels
  sig_label <- function(p) {
    if (is.na(p)) return("")
    if (p < 0.01) return("***")
    if (p < 0.05) return("**")
    if (p < 0.10) return("*")
    return("")
  }

  sig_mat <- matrix(sapply(p_adj, sig_label),
                    nrow = nrow(p_adj),
                    dimnames = dimnames(p_adj))

  # Build long-format plot data
  `%>%` <- magrittr::`%>%`
  df_corr <- as.data.frame(corr_coef) %>%
    dplyr::mutate(Row = rownames(.)) %>%
    tidyr::pivot_longer(-Row, names_to = "Col", values_to = "corr")

  df_p <- as.data.frame(p_adj) %>%
    dplyr::mutate(Row = rownames(.)) %>%
    tidyr::pivot_longer(-Row, names_to = "Col", values_to = "p_adj")

  df_sig <- as.data.frame(sig_mat) %>%
    dplyr::mutate(Row = rownames(.)) %>%
    tidyr::pivot_longer(-Row, names_to = "Col", values_to = "sig")

  df_n <- as.data.frame(n_mat) %>%
    dplyr::mutate(Row = rownames(.)) %>%
    tidyr::pivot_longer(-Row, names_to = "Col", values_to = "n")

  plot_df <- df_corr %>%
    dplyr::left_join(df_p,  by = c("Row", "Col")) %>%
    dplyr::left_join(df_sig, by = c("Row", "Col")) %>%
    dplyr::left_join(df_n,  by = c("Row", "Col"))

  # Add triangle indicator
  plot_df$tri <-
    ifelse(match(plot_df$Row, tissue_names) > match(plot_df$Col, tissue_names),
           "lower", "upper")
  plot_df$fill <-ifelse(plot_df$tri == "lower", plot_df$corr, plot_df$n)
  plot_df$lab <-ifelse(plot_df$tri == "lower", plot_df$sig, plot_df$n)

  # Make sure tiles are ordered the same as specified in tissue_names
  plot_df$Col <- factor(as.character(plot_df$Col), levels = tissue_names)
  plot_df$Row <- factor(as.character(plot_df$Row), levels = rev(tissue_names))

    # Combined heatmap
  return(ggplot2::ggplot(plot_df) +
    # Lower triangle: correlation heatmap
      ggplot2::geom_tile(ggplot2::aes(Col, Row, fill = fill),
                         dplyr::filter(plot_df, tri == "lower"),
              color = "white", na.rm = TRUE) +

    # Significance stars on lower triangle
      ggplot2::geom_text(ggplot2::aes(Col, Row, label = lab),
                         dplyr::filter(plot_df, tri == "lower"),
              color = "green", size = 4) +

    # Color scales
      ggplot2::scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
      na.value = "grey80",
      name = "Correlation"
    ) +
    ggnewscale::new_scale_fill() +

    # Upper triangle: sample size n
      ggplot2::geom_tile(ggplot2::aes(Col, Row, fill = fill),
                         dplyr::filter(plot_df, tri == "upper"),
              color = "white", na.rm = TRUE) +
      ggplot2::scale_fill_gradient(
      low = "yellow", high = "forestgreen",
      na.value = "grey80",
      name = "Number of shared features"
    ) +
    # Text on lower triangle
      ggplot2::geom_text(ggplot2::aes(Col, Row, label = lab),
                         dplyr::filter(plot_df, tri == "upper"),
              color = "grey30", size = 3) +

      ggplot2::labs(x = "", y = "",
         title = main_title) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank()))
}

# [END]
