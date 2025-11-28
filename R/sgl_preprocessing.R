#' Formatting Gene Expression and Pathway Data for Sparse Group Lasso
#'
#' This function prepares imputed gene expression data, phenotype data, and
#' pathway annotations for downstream sparse group lasso analysis.
#'
#' It performs the following:
#'   1. Sorts gene expression data by IID and removes FID/IID columns.
#'   2. Subsets and aligns phenotype and covariate files to the same individuals
#'      and in the same order as `gene_expr`.
#'   3. If covariates are provided:
#'        • Fits `glm(phenotype ~ covariates, family = family_func)`
#'        • Returns offsets = predicted values (for non-Gaussian families)
#'        • Returns residualized phenotype (for Gaussian families only)
#'   4. Expands pathways so genes appearing in multiple pathways are duplicated,
#'      renamed as `gene_pathway`.
#'   5. Constructs pathway group indices for sparse group lasso.
#'
#' @references
#' Friedman, J., Hastie, T., Tibshirani, R. (2010). Regularization Paths for
#' Generalized Linear Models via Coordinate Descent. *Journal of Statistical
#' Software.*
#' Meier, L., van de Geer, S., Bühlmann, P. (2008). The Group Lasso for Logistic
#' Regression. *J. Royal Statistical Society B.*
#'
#' @param gene_expr Data frame containing imputed expression data.
#'   Must contain first two columns named **FID** and **IID** and remaining
#'   columns representing gene expression.
#' @param all_pathways A named list where each element is a character vector
#'   of genes belonging to a pathway. Names represent pathway names.
#' @param phenotype_filename Path to phenotype text file containing at least
#'   the column `IID` and the phenotype column specified by `phenotype_colname`.
#' @param phenotype_colname Name of the phenotype column in phenotype file.
#' @param covariates_filename (Optional) path to covariates file that must
#'   contain `IID` and columns listed in `covariates_colnames`.
#' @param covariates_colnames Vector of covariate column names to extract.
#' @param family_func Model family specification (i.e., `"gaussian"`, `"binomial"`).
#'   Passed to `glm()`. See `?stats::family` for details.
#'
#' @return A list containing:
#'   \describe{
#'     \item{X}{Processed gene expression matrix with duplicated pathway genes.}
#'     \item{y}{Vector of phenotype values (residualized if Gaussian).}
#'     \item{groups}{Integer vector defining pathway grouping for group lasso.}
#'     \item{gene2pathway}{List mapping genes to the pathways they belong to.}
#'     \item{IID}{Participant IDs used (sorted).}
#'     \item{offsets}{NULL (Gaussian) or vector of glm predicted offsets
#'                   (non-Gaussian).}
#'     \item{processed_pathways}{Named list of pathway-specific duplicated genes.}
#'   }
#'
#' @examples
#' \dontrun{
#' preprocess_expressions_pathways(
#'   gene_expr = data.frame(
#'     FID = 1:3, IID = 1:3,
#'     G1 = rnorm(3), G2 = rnorm(3), G3 = rnorm(3)
#'   ),
#'   all_pathways = list(
#'     PW1 = c("G1", "G2"),
#'     PW2 = c("G2", "G3")
#'   ),
#'   phenotype_filename = "pheno.txt",
#'   phenotype_colname = "BMI",
#'   covariates_filename = "covars.txt",
#'   covariates_colnames = c("age", "sex"),
#'   family_func = "gaussian"
#' )
#' }
#'
#' @export
#' @import data.table

preprocess_expressions_pathways <- function(
    gene_expr,
    all_pathways,
    phenotype_filename = NULL,
    phenotype_colname = NULL,
    covariates_filename = NULL,
    covariates_colnames = NULL,
    family_func = "gaussian") {

  # Input validation
  if (missing(gene_expr)) stop("gene_expr must be provided.")
  if (!is.data.frame(gene_expr)) stop("gene_expr must be a data frame.")
  if (!identical(colnames(gene_expr)[1:2], c("FID", "IID")))
    stop("First 2 columns of gene_expr must be FID and IID columns.")

  if (!is.list(all_pathways) || is.null(names(all_pathways)))
    stop("all_pathways must be a named list of gene vectors.")

  if (!is.null(phenotype_filename) && is.null(phenotype_colname))
    stop("phenotype_colname must be provided when phenotype_filename is given.")

  valid_families <- c("gaussian", "binomial")
  if (!family_func %in% valid_families)
    stop("family_func must be one of: ",
         paste(valid_families, collapse = ", "), ".")


  # Process gene expression matrix (sort, strip FID/IID)
  gene_expr <- as.data.frame(gene_expr)
  gene_expr <- gene_expr[order(gene_expr$IID), ]
  rownames(gene_expr) <- gene_expr$IID
  gene_expr <- gene_expr[, -(1:2)]  # remove FID/IID columns

  # Load and align phenotype (if provided)
  if (is.null(phenotype_filename)) {
    phenotype <- NULL
  } else {
    phenotype <- data.table::fread(phenotype_filename)
    phenotype <- phenotype[!duplicated(phenotype$IID), ]
    phenotype <- as.data.frame(phenotype)

    # align to gene_expr individuals
    phenotype <- phenotype[match(rownames(gene_expr), phenotype$IID), ]
    phenotype <- phenotype[!is.na(phenotype$IID), ]
  }

  # Load covariates + compute glm offsets or residuals
  offsets <- NULL

  if (!is.null(covariates_filename)) {

    covars <- data.table::fread(covariates_filename)
    covars <- covars[!duplicated(covars$IID), ]
    covars <- as.data.frame(covars)

    # align covars, phenotype, gene_expr
    covars <- covars[match(phenotype$IID, covars$IID), ]
    covars <- covars[!is.na(covars$IID), ]

    # re-align phenotype
    phenotype <- phenotype[match(covars$IID, phenotype$IID), ]

    covars <- covars[, covariates_colnames, drop = FALSE]

    # Fit GLM
    fit <- glm(
      formula = phenotype[, phenotype_colname] ~ .,
      data    = covars,
      family  = family_func
    )

    offsets <- predict(fit, covars)

    # For Gaussian, residualize the phenotype instead of storing offsets
    if (family_func == "gaussian") {
      phenotype[, phenotype_colname] <- phenotype[, phenotype_colname] - offsets
      offsets <- NULL
    }

    # Reorder gene expression rows accordingly
    gene_expr <- gene_expr[match(phenotype$IID, rownames(gene_expr)), ]
  }

  # Extract phenotype vector
  if (!is.null(phenotype)) {
    y <- as.vector(phenotype[, phenotype_colname])
    if(family_func == "binomial"){
      y <- as.factor(y)
    }
  } else {
    y <- NULL
  }

  # Process pathways: duplicate genes, rename gene_pathway
  all_genes <- colnames(gene_expr)
  gene_pathways <- list()
  proc_pathways <- list()

  for (p in seq_along(all_pathways)) {

    pathway_name <- names(all_pathways)[p]
    genes <- intersect(all_pathways[[p]], all_genes)

    if (length(genes) == 0) next

    duplicated_genes <- character(length(genes))

    for (g in seq_along(genes)) {
      gene <- genes[g]

      # record membership
      if (gene %in% names(gene_pathways)) {
        gene_pathways[[gene]] <- c(gene_pathways[[gene]], length(genes))
        names(gene_pathways[[gene]])[length(gene_pathways[[gene]])] <- pathway_name
      } else {
        gene_pathways[[gene]] <- c(length(genes))
        names(gene_pathways[[gene]]) <- pathway_name
      }

      # duplicate + label as gene_pathway
      duplicated_genes[g] <- paste0(gene, "_", pathway_name)
    }

    proc_pathways[[pathway_name]] <- duplicated_genes
  }

  # Add singleton “pathways” for genes in no pathway
  for (g in all_genes) {
    if (!g %in% names(gene_pathways)) {
      gene_pathways[[g]] <- 1
      names(gene_pathways[[g]]) <- g
      proc_pathways[[g]] <- paste0(g, "_", g)
    }
  }

  # Build processed gene expression matrix with duplicated columns
  flat_path_genes <- unname(unlist(proc_pathways))

  gene_expr_proc <- gene_expr[, match(sub("_.*", "", flat_path_genes),
                                      colnames(gene_expr))]

  colnames(gene_expr_proc) <- flat_path_genes

  # pathway group indices
  groups <- rep(
    seq_along(proc_pathways),
    lengths(proc_pathways)
  )

  # Final output
  return(list(
    X = as.matrix(gene_expr_proc),
    y = y,
    groups = groups,
    gene2pathway = gene_pathways,
    IID = rownames(gene_expr_proc),
    offsets = offsets,
    processed_pathways = proc_pathways))
}

# [END]
