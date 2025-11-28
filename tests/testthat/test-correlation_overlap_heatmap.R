

test_that("Non-numeric betas throw an error.", {
  betas <- data.frame(A = c("x", "y"), B = c("a", "b"))
  expect_error(
    correlation_overlap_heatmap(betas),
    "`betas` must contain numeric values."
  )
})

test_that("Non data frame / matrix input throws an error.", {
  expect_error(
    correlation_overlap_heatmap(list(a = 1:3))
  )
})

test_that("Empty column matrix throws an error.", {
  betas <- matrix(nrow = 3, ncol = 0)
  expect_error(
    correlation_overlap_heatmap(betas)
  )
})

test_that("Invalid correlation method throws an error.", {
  betas <- data.frame(A = 1:3, B = 3:1)
  expect_error(
    correlation_overlap_heatmap(betas, cor_method = "foo"),
    "`cor_method` must be one of"
  )
})

test_that("Mismatched tissue_names length throws error.", {
  betas <- data.frame(A = 1:3, B = 3:1)
  expect_error(
    correlation_overlap_heatmap(betas, tissue_names = c("A")),
    "`tissue_names` must have length equal to ncol"
  )
})

test_that("One-column input returns count of non-NA elements.", {
  betas <- data.frame(A = c(1, NA, 3, NA))
  result <- correlation_overlap_heatmap(betas)
  expect_equal(result, 2)
})

test_that("Function runs and returns a ggplot object for valid input.", {
  betas <- data.frame(
    A = c(1, 7, 3),
    B = c(1, 6, 3),
    C = c(10, -9, NA)
  )

  p <- correlation_overlap_heatmap(betas)

  expect_s3_class(p, "ggplot")
})

test_that("Tissue names override column names in output plot.", {
  betas <- data.frame(
    X = c(1, 2, 3),
    Y = c(2, 3, 4)
  )

  custom_names <- c("T1", "T2")
  p <- correlation_overlap_heatmap(betas, tissue_names = custom_names)

  # Extract factor levels from the plot
  df <- p$data
  expect_true(all(levels(df$Col) == custom_names))
  expect_true(all(levels(df$Row) == rev(custom_names)))
})

test_that("Lower vs upper triangle labeling works.", {
  betas <- data.frame(
    A = c(1, 2, 3),
    B = c(2, 4, 6)
  )

  p <- correlation_overlap_heatmap(betas)
  df <- p$data

  # lower triangle = correlation values
  lower <- subset(df, tri == "lower")
  upper <- subset(df, tri == "upper")

  expect_true(all(lower$fill == lower$corr))
  expect_true(all(upper$fill == upper$n))
})
