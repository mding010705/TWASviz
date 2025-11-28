library(parallel)

test_that("No error on good input.", {
  expect_no_error(TWASviz::corr_heatmap(gene_expr =
                                  matrix(rnorm(500^2), nrow = 500)))
})


test_that("Error on empty input.", {
  expect_error(TWASviz::corr_heatmap(gene_expr = matrix()))
})

# Helper cluster

make_test_cluster <- function() {
  parallel::makeCluster(1, type = "PSOCK")
}

test_that("reorder_corr_hclust returns reordered matrix with same
          dims and names.", {
  mat <- matrix(c(1,0.2,0.2,
                  0.2,1,0.9,
                  0.2,0.9,1), nrow=3)
  rownames(mat) <- colnames(mat) <- c("A","B","C")

  out <- reorder_corr_hclust(mat)

  expect_true(is.matrix(out))
  expect_equal(dim(out), dim(mat))
  expect_setequal(rownames(out), rownames(mat))
  expect_setequal(colnames(out), colnames(mat))
})

test_that("redim_matrix returns reordered full matrix when target
          larger than original.", {
  mat <- matrix(runif(9), nrow = 3)
  rownames(mat) <- colnames(mat) <- LETTERS[1:3]

  cl <- make_test_cluster()

  out <- redim_matrix(mat, cl = cl, target_height = 10, target_width = 10)

  parallel::stopCluster(cl)

  # Should simply return reordered matrix of same size
  expect_equal(dim(out), c(3,3))
})

test_that("redim_matrix downsamples correctly.", {
  set.seed(1)
  mat <- matrix(runif(100), nrow=10)

  cl <- make_test_cluster()

  out <- redim_matrix(
    mat, cl = cl,
    target_height = 4,
    target_width = 5,
    summary_func = function(x) mean(x)
  )

  parallel::stopCluster(cl)

  expect_equal(dim(out), c(4,5))
})

test_that("redim_matrix uses supplied summary_func.", {
  mat <- matrix(1:16, nrow=4)

  cl <- make_test_cluster()

  # summary_func selects max element in each block
  out <- redim_matrix(
    mat,
    cl = cl,
    target_height = 2,
    target_width = 2,
    summary_func = function(x) max(x)
  )

  parallel::stopCluster(cl)

  expect_equal(out[1,1], max(mat[1:2,1:2]))
})

# main corr_heatmap function

test_that("corr_heatmap errors on invalid gene_expr inputs.", {
  expect_error(corr_heatmap(1), "gene_expr must be a numeric matrix")
  expect_error(corr_heatmap(matrix(1,1,1)), "at least 2 rows and 2 columns")
})


test_that("corr_heatmap uses custom summary_func.", {
  gene_expr <- matrix(rnorm(100), nrow=10)

  # summary_func returns maximum of block
  res <- corr_heatmap(
    gene_expr,
    summary_func = function(x) max(x),
    n_core = 1,
    target_height = 4,
    target_width = 4
  )

  # z contains block summaries
  expect_true(max(res$z) <= 1 && min(res$z) >= -1)
})

