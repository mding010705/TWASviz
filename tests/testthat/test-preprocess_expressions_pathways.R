test_that("no error on good input", {
  expect_no_error(TWASviz::preprocess_expressions_pathways(
    gene_expr = data.frame(FID = c(1,2,3), IID = c(1,1,1), gene1 = 1:3,
                                 gene2 = 1:3),
    all_pathways = list(p1 = c("gene1", "gene2"), p2 = c("gene4"))))
})


test_that("error on empty input", {
  expect_error(TWASviz::preprocess_expressions_pathways(
    gene_expr = data.frame(FID = c(), IID = c()),
    all_pathways = list(p1 = c(), p2 = c())))
})
