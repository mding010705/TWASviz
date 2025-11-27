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
#' @importFrom stringr str_split_fixed

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


#' Export Sparse Group Lasso (SGL) Gene and Pathway Coefficients to Text Files
#'
#' Takes a vector of gene-level coefficients from an SGL TWAS model, removes
#' intercept terms, splits combined `gene_pathway` names, aggregates gene-level
#' and pathway-level effects, and writes two text files:
#'   - `gene<suffix>.txt`
#'   - `pathway<suffix>.txt`
#'
#' @param gene_coef Named numeric vector of SGL coefficients.
#'   Names must follow the format `GENE_PATHWAY`, e.g., `"BRCA1_DNArepair"`.
#'   If a pathway is missing (e.g., `"BRCA1_"`), the gene name is used as the pathway name.
#'
#' @param filename_prefix Character prefix to prepend to the output files.
#'
#' @return A list containing:
#' \describe{
#'   \item{gene}{Data frame of aggregated gene-level effects}
#'   \item{pathway}{Data frame of aggregated pathway-level effects}
#' }
#'
#'
#' @examples
#' # Example SGL coefficient vector
#' coef_vec <- c(
#'   gene1_p1 = 0.5,
#'   gene2_p1 = -0.3,
#'   gene3_p2 = 0.2,
#'   "(Intercept)" = 1.2
#' )
#'
#' sgl2txt_file(coef_vec, filename_prefix = "demo_")
#'
#' @export
#' @import stringr

sgl2txt_file <- function(gene_coef,
                         filename_prefix = "sgl_TWAS_") {

  # Input Validation
  if (missing(gene_coef) || length(gene_coef) == 0)
    stop("gene_coef is empty. Provide a named numeric vector of coefficients.")

  if (is.null(rownames(gene_coef)))
    stop("gene_coef must be a row named numeric matrix.")

  gene_coef <- as.matrix(gene_coef)

  # Remove intercept term
  if ("(Intercept)" %in% rownames(gene_coef)) {
    gene_coef <- gene_coef[rownames(gene_coef) != "(Intercept)",]
  }

  # Keep only non-zero coefficients
  gene_coef <- gene_coef[(gene_coef != 0)]


  # All coefficients are zero
  if (length(gene_coef) == 0) {
    message("All coefficients are 0. No output files created.")
    return(NULL)
  }

  # Split coefficient names into gene and pathway
  # Expects: GENE_PATHWAY
  gene_pathway <- stringr::str_split_fixed(names(gene_coef), "_", n = 2)

  # If pathway is missing, use the gene name as pathway
  missing_pathway <- gene_pathway[, 2] == "" | gene_pathway[, 2] == " "
  gene_pathway[missing_pathway, 2] <- gene_pathway[missing_pathway, 1]

  # Aggregate at the gene level
  gene_df <- data.frame(
    gene = gene_pathway[, 1],
    betas = gene_coef,
    stringsAsFactors = FALSE
  )
  gene_df <- aggregate(betas ~ gene, data = gene_df, sum)

  # Aggregate at the pathway level
  pathway_df <- data.frame(
    pathway = gene_pathway[, 2],
    betas = gene_coef,
    stringsAsFactors = FALSE
  )
  pathway_df <- aggregate(betas ~ pathway, data = pathway_df, sum)

  # Write output files
  write.table(gene_df,
                     paste0(filename_prefix, "gene.txt"),
                     row.names = FALSE,
                     quote = FALSE,
                     sep = "\t")

  write.table(pathway_df,
                     paste0(filename_prefix, "pathway.txt"),
                     row.names = FALSE,
                     quote = FALSE,
                     sep = "\t")

  #Return structured output
  return(list(gene = gene_df, pathway = pathway_df))
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
#' @param gene_colname Name of column in predixcan_assoc_filenames that contains
#' the gene name.
#' @param effect_size_colname Name of column in predixcan_assoc_filenames that
#' contains the effect size.
#' @param pvalue_colname Name of column in predixcan_assoc_filenames that
#' contains the p-value.
#' @param tissue_names Names for tissues/conditions. Defaults to index sequence.
#' @param use_fdr If TRUE, convert p-values to q-values.
#' @param pvalue_thresh Threshold for p-value or q-value significance.
#'
#' @return A gene × tissue adjacency matrix of PrediXcan effect sizes.
#'
#' @examples
#' \dontrun{
#'   files <- c("inst/extdata/predixcan_twas1.txt",
#'   "inst/extdata/predixcan_twas2.txt")
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
                             gene_colname = "gene",
                             effect_size_colname = "zscore",
                             pvalue_colname = "pvalue",
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
  colnames(gene_adj_df) <- gene_colname

  # Read each PrediXcan file
  for (f in predixcan_assoc_filenames){

    # Loading using data.table
    assoc <- as.data.frame(data.table::fread(f, header = TRUE))

    # Check for p-value column, use q-value FDR correction if use_fdr is TRUE
    if (is.null(pvalue_colname)){
      filtered <- assoc[, c(gene_colname, effect_size_colname)]
    } else if (use_fdr == TRUE){
      assoc$qvalue <- qvalue::qvalue(assoc[, pvalue_colname])
      filtered <- assoc[assoc$qvalue < pvalue_thresh, c(gene_colname,
                                                        effect_size_colname)]
    } else {
      filtered <- assoc[assoc[, pvalue_colname] < pvalue_thresh, c(gene_colname,
                                                        effect_size_colname)]
    }

    # Iteratively merge into adjacency matrix
    gene_adj_df <- merge(gene_adj_df, filtered,
                         by = gene_colname, sort = FALSE, all = TRUE)
  }

  # Remove placeholder NA row
  gene_adj_df <- gene_adj_df[!is.na(gene_adj_df[, gene_colname]), ]

  # Assign user-specified tissue names
  colnames(gene_adj_df) <- c("gene", tissue_names)

  # Move gene column into rownames
  rownames(gene_adj_df) <- gene_adj_df$gene
  gene_adj_df <- gene_adj_df[, -1]

  return(gene_adj_df)
}

# [END]
