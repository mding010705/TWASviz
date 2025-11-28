library(data.table)

test_that("predixcan2adj_df rejects missing files.", {
  expect_error(predixcan2adj_df("blbla.txt"))
})

test_that("predixcan2adj_df processes single file.", {
  tmp <- tempfile(fileext = ".txt")
  dt <- data.frame(
    gene = c("G1", "G2"),
    zscore = c(2, 3),
    pvalue = c(0.01, 0.2)
  )
  fwrite(dt, tmp)

  res <- predixcan2adj_df(
    predixcan_assoc_filenames = tmp,
    tissue_names = "T1",
    pvalue_thresh = 0.05
  )

  expect_true("G1" %in% rownames(res))
  expect_false("G2" %in% rownames(res))
  expect_equal(res["G1", "T1"], 2)
})

test_that("predixcan2adj_df applies fdr correction.", {
  tmp <- tempfile(fileext = ".txt")

  dt <- data.frame(
    gene = paste0("G", 1:10),
    zscore = rnorm(10),
    pvalue = rep(c(0.2,0.1), 5) # should not pass fdr 0.05
  )
  fwrite(dt, tmp)

  res <- predixcan2adj_df(
    predixcan_assoc_filenames = tmp,
    tissue_names = "T",
    use_fdr = TRUE,
    pvalue_thresh = 0.05
  )

  expect_equal(nrow(res), 0)  # no gene passes fdr
})

test_that("predixcan2adj_df merges multiple files correctly", {
  tmp1 <- tempfile(fileext = ".txt")
  tmp2 <- tempfile(fileext = ".txt")

  dt1 <- data.frame(gene = "G1", zscore = 1, pvalue = 0.01)
  dt2 <- data.frame(gene = "G2", zscore = 2, pvalue = 0.01)

  fwrite(dt1, tmp1)
  fwrite(dt2, tmp2)

  res <- predixcan2adj_df(
    predixcan_assoc_filenames = c(tmp1, tmp2),
    tissue_names = c("T1", "T2"),
    pvalue_thresh = 0.05
  )

  expect_equal(res["G1", "T1"], 1)
  expect_true(is.na(res["G1", "T2"]))
  expect_equal(res["G2", "T2"], 2)
})
