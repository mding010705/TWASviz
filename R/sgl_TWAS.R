#' Perform Sparse Group Lasso TWAS
#'
#' Performs sparse group lasso with the formatted imputed gene expression data,
#' phenotype data, groupings, and offsets from sgl_preprocessing().
#'
#' @param X Imputed gene expression matrix.
#' @param y Vector of phenotype of interest.
#' @param grouping Vector of consecutive integers describing the pathway
#' grouping of the genes.
#' @param family_func Family object specifying the type of model you
#' want to use in the sparse group lasso function. See ?stats::family for more
#' details.
#' @param pred_loss Loss to use for cross-validation error. Valid options are:
#' "default" the same as deviance (mse for regression and deviance otherwise)
#' "mse" mean square error
#' "deviance" the default (mse for Gaussian regression, and negative
#' log-likelihood otherwise)
#' "mae" mean absolute error, can apply to any family
#' "misclass" for classification only, misclassification error.
#' @param offsets Optional offset (constant predictor without a corresponding
#' coefficient)
#' @param nfolds Number of folds, ranging from 3 to
#' sample size (leave-one-out CV).
#'
#' @return Returns a list(
#' sgl_fit = cv.sparsegl object from sparsegl::cv.sparsegl(),
#' lambda1se_coef = dgCMatrix object containing gene names and coefficients
#' estimated from the sparse group lasso with lambda = largest value of lambda
#' such that error is within 1 standard error of the minimum cross-validated
#' errors).
#'
#' @references
#' Liang, X., Cohen, A., Sólon Heinsfeld, A., Pestilli, F., and McDonald, D.J.
#' 2024. sparsegl: An R Package for Estimating Sparse Group Lasso. Journal of
#' Statistical Software, Vol. 110(6): 1–23. doi:10.18637/jss.v110.i06.
#'
#' @export
#' @import sparsegl

sgl_TWAS <- function(X, y, grouping = NULL,
                     family_func = "gaussian",
                     pred_loss = "mse", offsets = NULL, nfolds = 10){
  sgl_fit <- sparsegl::cv.sparsegl(X, y, group = grouping, family = family_func,
                                   pred.loss = pred_loss, offset = offsets,
                                   nfolds = nfolds, intercept = TRUE)

  return(list(sgl_fit = sgl_fit,
              lambda1se_coef = coef(sgl_fit, s = "lambda.1se")))
}

# [END]
