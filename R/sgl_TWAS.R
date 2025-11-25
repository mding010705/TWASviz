#' Perform Sparse Group Lasso TWAS
#'
#' Performs a Transcriptome-Wide Association Study (TWAS) using the
#' Sparse Group Lasso (SGL), given imputed gene expression data, phenotype
#' values, pathway-based group assignments, and optional offsets.
#'
#' Internally, this function fits a cross-validated Sparse Group Lasso model
#' using `sparsegl::cv.sparsegl()`, then extracts the coefficients
#' at the *lambda.1se* value (the most regularized model within one standard
#' error of the minimum CV error), typically used for stability and
#' interpretability in TWAS.
#'
#' **References**
#' - Liang, X., Cohen, A., Sólon Heinsfeld, A., Pestilli, F., & McDonald, D.J. (2024).
#'   *sparsegl: An R Package for Estimating Sparse Group Lasso.*
#'   Journal of Statistical Software, 110(6), 1–23.
#'   doi:10.18637/jss.v110.i06
#'
#' - Simon, N., Friedman, J., Hastie, T., & Tibshirani, R. (2013).
#'   *A Sparse-Group Lasso.* Journal of Computational and Graphical Statistics, 22(2), 231–245.
#'
#'
#' @param X Numeric matrix of imputed gene expression values,
#'   where rows are samples and columns are genes.
#' @param y Numeric vector of phenotype values (continuous or binary).
#' @param grouping Integer vector assigning each gene to a pathway group.
#'   Must be consecutive integers (1, 2, 3, …).
#'   If `NULL`, the sparse group lasso reduces to standard lasso.
#' @param family_func Model family specification (e.g., `"gaussian"`, `"binomial"`).
#'   Passed to `cv.sparsegl()`. See `?stats::family` for details.
#' @param pred_loss Loss function for cross-validation. Valid options:
#'   `"default"`, `"mse"`, `"deviance"`, `"mae"`, `"misclass"`.
#' @param offsets Optional numeric vector of offsets. Useful for models with
#'   baseline risk or exposure variables.
#' @param nfolds Number of cross-validation folds. Must be ≥3.
#'
#' @return A list containing:
#' \describe{
#'   \item{sgl_fit}{The full `cv.sparsegl` model object.}
#'   \item{lambda1se_coef}{A sparse coefficient matrix (`dgCMatrix`)
#'         containing gene names and coefficients at the `lambda.1se` value.}
#' }
#'
#' @examples
#' \dontrun{
#'   # Simulated data for demonstration
#'   set.seed(123)
#'   X <- matrix(rnorm(200 * 50), nrow = 200, ncol = 50)
#'   y <- rnorm(200)
#'
#'   # Group 50 genes into 10 pathways
#'   grouping <- rep(1:10, each = 5)
#'
#'   # Run TWAS using Sparse Group Lasso
#'   twas_res <- sgl_TWAS(
#'     X = X,
#'     y = y,
#'     grouping = grouping,
#'     family_func = "gaussian",
#'     pred_loss = "mse",
#'     nfolds = 5
#'   )
#'
#'   # Inspect fitted object
#'   twas_res$sgl_fit
#'
#'   # Extract gene coefficients at lambda.1se
#'   coef_1se <- twas_res$lambda1se_coef
#'   head(coef_1se)
#' }
#'
#' @export
#' @import sparsegl


sgl_TWAS <- function(X, y, grouping = NULL,
                     family_func = "gaussian",
                     pred_loss = "mse", offsets = NULL, nfolds = 10){
  # Input Validation
  # Validate X
  if (!is.matrix(X))
    stop("X must be a numeric matrix.")

  if (!is.numeric(X))
    stop("X must contain numeric values.")

  # Validate y
  if (!is.numeric(y))
    stop("y must be a numeric vector.")

  if (length(y) != nrow(X))
    stop("Length of y (", length(y),") must match number of rows in X
         (", nrow(X), ").")

  # Validate grouping
  if (!is.null(grouping)) {
    if (!is.numeric(grouping) || !is.integer(grouping))
      stop("grouping must be an integer vector.")

    if (length(grouping) != ncol(X))
      stop("Length of grouping must match number of columns in X (genes).")

    # Check consecutive group IDs (1, 2, 3, ...)
    if (!identical(sort(unique(grouping)), seq_along(unique(grouping))))
      stop("grouping must consist of consecutive integers (e.g., 1, 2, 3, ...).")
  }

  # Validate family_func
  valid_families <- c("gaussian", "binomial", "poisson")
  if (!family_func %in% valid_families)
    stop("family_func must be one of: ",
         paste(valid_families, collapse = ", "), ".")

  # Validate pred_loss
  valid_losses <- c("default", "mse", "deviance", "mae", "misclass")
  if (!pred_loss %in% valid_losses)
    stop("pred_loss must be one of: ",
         paste(valid_losses, collapse = ", "), ".")

  # Validate offsets
  if (!is.null(offsets)) {
    if (!is.numeric(offsets))
      stop("offsets must be numeric.")

    if (length(offsets) != length(y))
      stop("offsets must have the same length as y.")
  }

  # Validate nfolds
  if (!is.numeric(nfolds) || nfolds < 3)
    stop("nfolds must be a numeric value >= 3.")

  # Fit Sparse Group Lasso with cross-validation.
  sgl_fit <- sparsegl::cv.sparsegl(
    X, y,
    group      = grouping,
    family     = family_func,
    pred.loss  = pred_loss,
    offset     = offsets,
    nfolds     = nfolds,
    intercept  = TRUE
  )

  # Extract coefficients at lambda.1se (one-standard-error rule)
  lambda1se_coef <- coef(sgl_fit, s = "lambda.1se")

  # Return list containing both the full model and selected coefficients
  return(list(
    sgl_fit        = sgl_fit,
    lambda1se_coef = lambda1se_coef))
}


# [END]
