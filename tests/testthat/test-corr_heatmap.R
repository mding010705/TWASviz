test_that("no error on good input", {
  expect_no_error(TWASviz::corr_heatmap(gene_expr =
                                  matrix(rnorm(500^2), nrow = 500)))
})


test_that("error on empty input", {
  expect_error(TWASviz::corr_heatmap(gene_expr = matrix()))
})
