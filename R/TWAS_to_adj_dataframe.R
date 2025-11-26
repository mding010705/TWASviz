#' Format Sparse Group Lasso Coefficients into Gene/Pathway Adjacency Matrices
#'
#' Converts coefficient output from multiple Sparse Group Lasso (SGL) models
#' into two adjacency-matrix-like data frames:
#' 1. genes × tissues and
#' 2. pathways × tissues.
#'
#' For each SGL model, nonzero coefficients are parsed into gene and pathway
#' components. Coefficients belonging to the same gene or pathway are summed.
#' Intercepts are removed.
#'
#' **References**
#' - Simon, N., Friedman, J., Hastie, T., & Tibshirani, R. (2013).
#'   *A Sparse-Group Lasso.* Journal of Computational and Graphical Statistics, 22(2), 231–245.
#' - `sparsegl` package documentation: https://cran.r-project.org/package=sparsegl
#'
#' @param gene_coefs A list of coefficient vectors (e.g., from `sparsegl::cv.sparsegl()` or sgl_TWAS).
#'   Each element should be a named numeric vector with coefficients for one tissue/condition.
#' @param tissue_names Character vector giving column names for each coefficient set.
#'   Defaults to sequence along `gene_coefs`.
#'
#' @return A list with two data frames:
#' \describe{
#'   \item{gene}{Gene × tissue/cell type adjacency matrix of summed effect sizes.}
#'   \item{pathway}{Pathway × tissue/cell type adjacency matrix of summed effect sizes.}
#' }
#'
#' @examples
#' \dontrun{
#'   # Example coefficient vectors (normally from cv.sparsegl$beta or sgl_TWAS$lambda1se_coef)
#'   coef1 <- c("(Intercept)" = 0, "GENE1_PATHA" = 0.5, "GENE2_PATHB" = -0.2)
#'   coef2 <- c("(Intercept)" = 0, "GENE1_PATHA" = 0.1, "GENE3_PATHC" = 0.4)
#'
#'   out <- sgl2adj_df(
#'     gene_coefs = list(coef1, coef2),
#'     tissue_names = c("TissueA", "TissueB")
#'   )
#'
#'   out$gene
#'   out$pathway
#' }
#'
#' @export
#' @import stringr

sgl2adj_df <- function(gene_coefs,
                       tissue_names = seq_along(gene_coefs)){
  # Input Validation
  if (!is.list(gene_coefs))
    stop("gene_coefs must be a list of named numeric vectors.")

  if (length(gene_coefs) == 0)
    stop("gene_coefs is empty. Provide coefficient vectors.")

  # Validate tissue_names
  if (length(tissue_names) != length(gene_coefs))
    stop("Length of tissue_names must match length of gene_coefs.")

  if (any(duplicated(tissue_names)))
    stop("tissue_names must contain unique values.")

  # Initialize with placeholder NA to enable iterative merging
  gene_adj_df <- data.frame(gene = c(NA))
  pathway_adj_df <- data.frame(pathway = c(NA))

  # Loop through each set of coefficients (each tissue/cell type)
  for(coef in gene_coefs){

    # Remove intercept term
    coef <- coef[(rownames(coef) != "(Intercept)"), ]

    # Handle case where all coefficients are zero
    if (all(coef == 0)){
      gene_adj_df <- cbind(gene_adj_df, NA)
      pathway_adj_df <- cbind(pathway_adj_df, NA)
      message("A set of coefficients is all 0.")
      next
    }

    # Keep only non-zero coefficients
    coef <- coef[coef != 0]

    # Split coefficient names into gene and pathway components
    gene_pathway <- stringr::str_split_fixed(names(coef), "_", n = 2)

    # If pathway is missing, treat gene name as the pathway name
    gene_pathway[gene_pathway[, 2] == "", 2] <-
      gene_pathway[gene_pathway[, 2] == "", 1]

    # Aggregate Gene Effects
    gene_df <- data.frame(gene = gene_pathway[, 1], betas = coef)
    gene_df <- aggregate(. ~ gene, data = gene_df, FUN = sum)

    # Aggregate Pathway Effects
    pathway_df <- data.frame(pathway = gene_pathway[, 2], betas = coef)
    pathway_df <- aggregate(. ~ pathway, data = pathway_df, FUN = sum)

    # Merge these columns into growing adjacency matrices
    gene_adj_df <- merge(gene_adj_df, gene_df, by = "gene", all = TRUE)
    pathway_adj_df <- merge(pathway_adj_df, pathway_df, by = "pathway", all = TRUE)
  }

  # Assign column names using provided tissue names
  colnames(gene_adj_df) <- c("gene", tissue_names)
  colnames(pathway_adj_df) <- c("pathway", tissue_names)

  # Remove initial NA placeholder row
  gene_adj_df <- gene_adj_df[!is.na(gene_adj_df$gene), ]
  pathway_adj_df <- pathway_adj_df[!is.na(pathway_adj_df$pathway), ]

  # Convert rownames to gene/pathway names
  if (!is.null(dim(gene_adj_df))){
    rownames(gene_adj_df) <- gene_adj_df$gene
    rownames(pathway_adj_df) <- pathway_adj_df$pathway

    # Drop gene/pathway columns. They are now rownames
    gene_adj_df <- gene_adj_df[, -1]
    pathway_adj_df <- pathway_adj_df[, -1]
  }

  return(list(gene = gene_adj_df, pathway = pathway_adj_df))
}





#' Format PrediXcan Association Output into Gene Adjacency Matrix
#'
#' Reads multiple PrediXcan association result files and constructs an
#' adjacency-matrix-like data frame where rows are genes and columns are
#' tissues/conditions, storing the effect size (standard error or beta).
#'
#' Filters genes using either a p-value threshold or FDR-based q-value
#' threshold (Storey & Tibshirani 2003).
#'
#' **References**
#' - Gamazon, E. R., et al. (2015).
#'   *PrediXcan: Trait prediction from transcriptome data integrates genetic regulation.*
#'   Nature Genetics, 47, 1091–1098.
#' - Storey, J. D., & Tibshirani, R. (2003).
#'   *Statistical significance for genomewide studies.* PNAS, 100(16), 9440–9445.
#'
#' @param predixcan_assoc_filenames Vector of file paths to PrediXcan output files.
#'   Each must contain columns: \code{gene}, \code{se}, and \code{pvalue}.
#' @param tissue_names Names for tissues/conditions. Defaults to index sequence.
#' @param use_fdr If TRUE, convert p-values to q-values.
#' @param pvalue_thresh Threshold for p-value or q-value significance.
#'
#' @return A gene × tissue adjacency matrix of PrediXcan effect sizes.
#'
#' @examples
#' \dontrun{
#'   files <- c("muscle_assoc.txt", "blood_assoc.txt")
#'
#'   adj <- predixcan2adj_df(
#'     predixcan_assoc_filenames = files,
#'     tissue_names = c("Cell1", "Cell2"),
#'     use_fdr = TRUE,
#'     pvalue_thresh = 0.10
#'   )
#'
#'   head(adj)
#' }
#'
#' @export
#' @import data.table
#' @import qvalue


predixcan2adj_df <- function(predixcan_assoc_filenames,
                             tissue_names = seq_along(predixcan_assoc_filenames),
                             use_fdr = FALSE,
                             pvalue_thresh = 0.05){
  # Input validation
  if (length(predixcan_assoc_filenames) == 0)
    stop("No PrediXcan filenames provided.")

  # Check files exist
  missing_files <- predixcan_assoc_filenames[!file.exists(predixcan_assoc_filenames)]
  if (length(missing_files) > 0)
    stop("The following PrediXcan files do not exist:\n",
         paste0("  - ", missing_files, collapse = "\n"))

  # Validate tissue_names
  if (length(tissue_names) != length(predixcan_assoc_filenames))
    stop("Length of tissue_names must match number of filenames.")

  if (any(duplicated(tissue_names)))
    stop("tissue_names must contain unique values.")


  # Initialize with placeholder NA to enable iterative merging
  gene_adj_df <- data.frame(gene = c(NA))

  # Read each PrediXcan file
  for (f in predixcan_assoc_filenames){

    # Efficient loading using data.table
    assoc <- as.data.frame(data.table::fread(f, header = TRUE))

    # Use q-value FDR correction if use_fdr is TRUE
    if (use_fdr == TRUE){
      assoc$qvalue <- qvalue::qvalue(assoc$pvalue)
      filtered <- assoc[assoc$qvalue < pvalue_thresh, c("gene", "se")]
    } else {
      filtered <- assoc[assoc$pvalue < pvalue_thresh, c("gene", "se")]
    }

    # Iteratively merge into adjacency matrix
    gene_adj_df <- merge(gene_adj_df, filtered,
                         by = "gene", sort = FALSE, all = TRUE)
  }

  # Remove placeholder NA row
  gene_adj_df <- gene_adj_df[!is.na(gene_adj_df$gene), ]

  # Assign user-specified tissue names
  colnames(gene_adj_df) <- c("gene", tissue_names)

  # Move gene column into rownames
  rownames(gene_adj_df) <- gene_adj_df$gene
  gene_adj_df <- gene_adj_df[, -1]

  return(gene_adj_df)
}

# [END]
