test_that("sgl2txt_file fails on invalid input.", {
  expect_error(sgl2txt_file(NULL))
  expect_error(sgl2txt_file(numeric(0)))
})

test_that("sgl2txt_file handles all zero coefficients.", {
  coef <- as.matrix(c("(Intercept)" = 1, "G1_P1" = 0), nrow = 2)

  expect_message(result <- sgl2txt_file(coef, filename_prefix = "tmp_"))
  expect_null(result)
})

test_that("sgl2txt_file aggregates correctly.", {
  coef <- as.matrix(c("(Intercept)" = 1, "G1_P1" = 1, "G1_P1" = 2,
                      "G2_P2" = -1), nrow = 4)

  res <- sgl2txt_file(coef, filename_prefix = "tmp_")

  expect_true("G1" %in% res$gene$gene)
  expect_equal(res$gene$betas[res$gene$gene == "G1"], 3)

  expect_true("P2" %in% res$pathway$pathway)
  expect_equal(res$pathway$betas[res$pathway$pathway == "P2"], -1)
})

test_that("sgl2txt_file writes output files.", {
  coef <- as.matrix(c("G1_P1" = 1))
  prefix <- tempfile(pattern = "sgltest_")

  res <- sgl2txt_file(coef, filename_prefix = prefix)

  expect_true(file.exists(paste0(prefix, "gene.txt")))
  expect_true(file.exists(paste0(prefix, "pathway.txt")))
})
