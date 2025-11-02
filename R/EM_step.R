
#' Expectation-Maximization (EM) Algorithm for TIPS-TWAS
#'
#' Executes the EM algorithm for TIPS. Estimates the effect sizes of
#' genes on a phenotype of interest in a pathway-aware way, using block-wise
#' gradient descent.
#'
#' @param old Initial parameter estimates (sigma1, sigma2, sigmau, alpha)
#' @param a_para Regularization parameter
#' @param wg1 Genotype matrix for group 1 (reference set)
#' @param wg2 Genotype matrix for group 2 (outcome set)
#' @param z Phenotype vector
#' @param y Gene expression matrix
#' @param lambda Penalty strength
#' @param k Vector of group sizes for each gene
#' @param m1 Number of genes
#' @param n1 Number of individuals in group 1
#' @param n2 Number of individuals in group 2
#' @param p1 Number of SNPs
#' @param max_i Max number of EM iterations
#' @param tol Convergence tolerance
#'
#' @return The final parameters (final estimates for sigma1, sigma2, sigmau,
#' alpha, lambda, log likelihood)


EM_updated <- function(old, a_para, wg1, wg2, z, y, lambda, k,
                       m1, n1, n2, p1, max_i = 10, tol = 0.001) {
  # init param containers
  theta_old <- vector("list", max_i)
  theta_history <- vector("list", max_i)
  current_diff <- rep(NA, max_i)
  prev_diff <- rep(Inf, max_i)

  theta_old[[1]] <- old[1:(3 + m1)]

  # EM iters
  for (i in seq_len(max_i)) {
    cat(sprintf("Start iteration %d in EM_updated.\n", i))

    # call M-step function
    step_results <- mstep_updated(theta_old[[i]], a_para, wg1, wg2, z, y,
                                  lambda, k, m1, n1, n2, p1)

    # extract params and likelihood
    theta_history[[i]] <- step_results[1:(3 + m1)]
    like <- step_results[(4 + m1):length(step_results)]

    # avg absolute difference between old and new parameter sets
    current_diff[i] <- mean(abs(theta_history[[i]] - theta_old[[i]]))

    # convergence criterion
    if (abs(current_diff[i] - prev_diff[i]) < tol) {
      cat(sprintf("Converged at iteration %d (|Δdiff| < tol).\n", i))
      final_parameters <- c(theta_history[[i]], lambda, like)
      cat("Finished EM_updated.\n")
      return(final_parameters)
    } else {
      theta_old[[i + 1]] <- theta_history[[i]]
      prev_diff[i + 1] <- current_diff[i]
      cat(sprintf("Finished iteration %d in EM_updated.\n", i))
    }
  }


  # if model doesn't converge
  cat("Reached maximum iterations without convergence.\n")
  final_parameters <- c(theta_history[[max_i]], lambda, like)
  cat("Finished EM_updated function.\n")

  return(final_parameters)
}
