



mstep_updated <- function(old, a_para, wg1, wg2, z, y,
                          lambda, k, m1, n1, n2, p1,
                          index_filename,
                          n_core = parallel::detectCores() - 1){
  # init parameters
  sigma1 <- old[1]
  sigma2 <- old[2]
  sigmau <- old[3]
  alpha_g <- old[4:length(old)]

  # read gene block indices
  index_gene_geno <- data.table::fread(index_filename, header = TRUE)
  index_gene_array <- as.numeric(index_gene_geno[[1]])
  start_indices <- c(1, cumsum(index_gene_array[-length(index_gene_array)]) + 1)
  end_indices <- cumsum(index_gene_array)

  # setup parallel proc
  cl <- parallel::makeCluster(n_core, type="PSOCK")
  parallel::clusterExport(cl, varlist = c("m1", "start_indices",
                                          "end_indices", "wg1", "wg2",
                                          "sigma1", "alpha_g", "sigma2",
                                          "sigmau", "z", "a_para",
                                          "lambda"),
                          envir = environment())

  # parallelized m-step, per gene
  results <- do.call(rbind, parallel::parSapply(cl, seq_len(m1), function(i) {
    cat(sprintf("Index: %d in m1 using mstep_updated function.\n", i))

    block_start <- start_indices[i]
    block_end <- end_indices[i]
    block_indices <- block_start:block_end

    wg1_block <- wg1[, block_indices, drop = FALSE]
    wg2_block <- wg2[, block_indices, drop = FALSE]

    wg1_product <- t(wg1_block) %*% wg1_block
    wg2_product <- t(wg2_block) %*% wg2_block

    col <- ncol(wg1_block)


    # make A matrix
    A <- (1 / sigma1) * wg1_product +
      (alpha_g[i]^2 / sigma2) * wg2_product +
      (1 / sigmau) * diag(col)

    # regularization
    A <- A + 1e-4 * diag(col)

    # solve for sigma_ui_i and mu_ui_i
    sigma_ui_i <- solve(A)
    mu_ui_i <- sigma_ui_i %*% ((1 / sigma1) * t(wg1_block) %*% y[, i] +
                                 (alpha_g[i] / sigma2) * t(wg2_block) %*% z)


    # alpha est
    si <- 2 * t(z) %*% wg2_block %*% mu_ui_i - a_para * lambda * sign(alpha_g[i])
    abc_i <- (1 - a_para) * lambda * k[i]

    if (abs(si) <= abc_i || abs(si) <= (a_para * lambda)) {
      alpha_est <- 0
    } else {
      denom <- (2 * t(mu_ui_i) %*% wg2_product %*% mu_ui_i +
                  2 * sum(diag(wg2_block %*% sigma_ui_i %*% t(wg2_block))))
      alpha_est <- (1 - abc_i / abs(si)) * as.numeric(si / denom)
    }

    # compute E1, E2, E3
    E1 <- as.numeric(t(y[, i]) %*% y[, i] -
                       2 * t(y[, i]) %*% wg1_block %*% mu_ui_i +
                       t(mu_ui_i) %*% wg1_product %*% mu_ui_i +
                       sum(diag(wg1_product %*% sigma_ui_i)))

    E2 <- as.numeric(t(z) %*% z -
                       2 * alpha_est * t(z) %*% wg2_block %*% mu_ui_i +
                       alpha_est^2 * (t(mu_ui_i) %*% wg2_product %*% mu_ui_i) +
                       sum(diag((alpha_est^2) * wg2_product %*% sigma_ui_i)))

    E3 <- as.numeric(t(mu_ui_i) %*% mu_ui_i + sum(diag(sigma_ui_i)))

    return(c(E1 = E1, E2 = E2, E3 = E3, alpha_est = alpha_est))
  }))

  stopCluster(cl)

  # final ests
  sigma1_est <- (1 / (m1 * n1)) * sum(results[, "E1"])
  sigma2_est <- (1 / (m1 * n2)) * sum(results[, "E2"])
  sigmau_est <- (1 / p1) * sum(results[, "E3"])

  log_likelihood_y <- -0.5 * n1 * m1 * log(2 * pi * sigma1_est) -
    0.5 * sum(results[, "E1"]) / sigma1_est
  log_likelihood_z <- -0.5 * n2 * log(2 * pi * sigma2_est) -
    0.5 * sum(results[, "E2"]) / sigma2_est
  log_likelihood_u <- -0.5 * p1 * m1 * log(2 * pi * sigmau_est) -
    0.5 * sum(results[, "E3"]) / sigmau_est

  log_likelihood <- log_likelihood_y + log_likelihood_z + log_likelihood_u

  all_E2_values <- results[, "E2"]
  E2_null <- as.numeric(t(z) %*% z)

  output <- c(sigma1_est, sigma2_est, sigmau_est,
              results[, "alpha_est"],
              all_E2_values,
              E2_null, log_likelihood)

  return(output)
}
