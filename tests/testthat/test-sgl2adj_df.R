test_that("sgl2adj_df fails on invalid input", {
  expect_error(sgl2adj_df(gene_coefs = NULL))
  expect_error(sgl2adj_df(gene_coefs = list(), tissue_names = "A"))
  expect_error(sgl2adj_df(gene_coefs = list(1, 2), tissue_names = "A"))
})

test_that("sgl2adj_df computes correct aggregation", {
  coef1 <- c("(Intercept)" = 0, "G1_P1" = 1, "G2_P2" = 2)
  coef2 <- c("(Intercept)" = 0, "G1_P1" = 3, "G3_P3" = 4)

  res <- sgl2adj_df(
    gene_coefs = list(coef1, coef2),
    tissue_names = c("T1", "T2")
  )

  # Correct genes
  expect_true(all(c("G1", "G2", "G3") %in% rownames(res$gene)))

  # Correct values
  expect_equal(as.numeric(res$gene["G1", ]), c(1, 3))
  expect_equal(res$gene["G2", "T2"], NA_real_)

  # Pathways
  expect_true("P3" %in% rownames(res$pathway))
})

test_that("sgl2adj_df handles all zero coefficients.", {
  coef1 <- c("(Intercept)" = 1, "G1_P1" = 0)
  coef2 <- c("(Intercept)" = 0, "G2_P2" = 5)

  res <- sgl2adj_df(
    gene_coefs = list(coef1, coef2),
    tissue_names = c("A", "B")
  )

  expect_true(is.na(res$gene["G1", "A"]))
  expect_equal(res$gene["G2", "B"], 5)
})
