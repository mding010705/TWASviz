
test_that("no error on good input", {
  expect_no_error(TWASviz::correlation_overlap_heatmap(
    betas = data.frame(x = c(1, 7, 3), y = c(1, -9, 3), z = c(-1, 0, NA)),
                                           tissue_names = c("A", "B", "C")))
})


test_that("error on no", {
  expect_error(TWASviz::correlation_overlap_heatmap())
})


test_that("null output on too few observations input", {
  expect_null(TWASviz::correlation_overlap_heatmap(
    betas = data.frame(x = c(NA, 7, NA), y = c(1, NA, 3), z = c(-1, 0, NA)),
                                          tissue_names = c("A", "B", "C")))
})
