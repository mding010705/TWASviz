test_that("no error on good input", {
  expect_no_error(TWASviz::volcano_plot(data.frame(beta = c(1, 2, 3),
                                                   p_value = c(0.5, 0.001, 1)),
                                        add_sig_gene_labels = FALSE))
})


test_that("error on empty input", {
  expect_error(TWASviz::volcano_plot())
})
