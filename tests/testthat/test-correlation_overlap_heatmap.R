
test_that("no error on good input", {
  expect_no_error(TWASviz::correlation_overlap_heatmap(
    betas = data.frame(x = c(1, 7, 3), y = c(1, -9, 3), z = c(-1, 0, NA)),
                                           tissue_names = c("A", "B", "C")))
})


test_that("error on no input", {
  expect_error(TWASviz::correlation_overlap_heatmap())
})


# Mock psych::corr.test() so test suite does not require psych

fake_corr_test <- function(betas, method, adjust) {

  p <- ncol(betas)

  # Simple artificial correlation matrix
  r <- diag(1, p)
  r[r == 0] <- 0.5

  # Fake p-values
  p_adj <- matrix(0.05, nrow = p, ncol = p)
  diag(p_adj) <- 0

  # Fake sample size matrix
  n_mat <- matrix(10, nrow = p, ncol = p)

  structure(
    list(
      r = r,
      p.adj = p_adj[lower.tri(p_adj)],
      n = n_mat
    ),
    class = "psych_corr_test"
  )
}

mock_psych <- function(expr) {
  with_mock(
    `psych::corr.test` = fake_corr_test,
    expr
  )
}


test_that("single column returns count of non-NA values.", {
  betas <- data.frame(A = c(1, NA, 3, 4))
  expect_message(res <- correlation_overlap_heatmap(betas))
  expect_equal(res, 3)
})

test_that("input must be numeric.", {
  betas <- data.frame(A = c("x", "y"))
  expect_error(correlation_overlap_heatmap(betas))
})

test_that("cor_method must be valid.", {
  betas <- data.frame(A = 1:3, B = 2:4)
  expect_error(correlation_overlap_heatmap(betas, cor_method = "wrong"))
})

test_that("tissue_names must match number of columns.", {
  betas <- data.frame(A = 1:3, B = 2:4)
  expect_error(correlation_overlap_heatmap(betas, tissue_names = c("A")))
})

test_that("function returns ggplot object when >= 2 columns.", {
  betas <- data.frame(
    A = c(1, 2, 3),
    B = c(1, 2, 3),
    C = c(1, 2, 3)
  )

  res <- mock_psych(
    correlation_overlap_heatmap(betas)
  )

  expect_true(inherits(res, "ggplot"))
})

test_that("mocked psych::corr.test called and combined output structured properly.", {
  betas <- data.frame(
    T1 = c(1, 2, NA),
    T2 = c(2, 3, 1),
    T3 = c(5, 6, 7)
  )

  res <- mock_psych(
    correlation_overlap_heatmap(betas)
  )

  expect_true(inherits(res, "ggplot"))

  # Extract plot data
  plot_data <- res$data

  # Expect rows = tissues^2
  expect_equal(nrow(plot_data), 3 * 3)

  # All expected fields exist
  expect_true(all(c("Row", "Col", "corr", "p_adj", "sig", "n", "tri") %in%
                    names(plot_data)))
})

test_that("significance stars computed correctly.", {
  # Build fake data where we can control p-values
  fake_corr_test2 <- function(betas, method, adjust) {

    p <- ncol(betas)

    r <- diag(1, p)
    p_vals <- c(0.005, 0.03, 0.08)  # upper triangle

    structure(
      list(
        r = r,
        p.adj = p_vals,
        n = matrix(10, p, p)
      ),
      class = "psych_corr_test"
    )
  }

  res <- with_mock(
    `psych::corr.test` = fake_corr_test2,
    correlation_overlap_heatmap(data.frame(A=1:3, B=1:3, C=1:3))
  )

  df <- res$data

  # sig stars only in lower triangle
  sig_values <- unique(df$sig)

  expect_true("***" %in% sig_values) # p=0.005
  expect_true("**"  %in% sig_values) # p=0.03
  expect_true("*"   %in% sig_values) # p=0.08
})

test_that("upper triangle uses n matrix.", {
  betas <- data.frame(
    A = c(1, 2, 3),
    B = c(4, 5, 6),
    C = c(7, 8, 9)
  )

  res <- mock_psych(
    correlation_overlap_heatmap(betas)
  )

  df <- res$data

  upp <- df[df$tri == "upper", ]

  # All upper triangle fill values should equal fake n (10)
  expect_true(all(upp$fill == 10))
})

test_that("lower triangle uses correlation matrix.", {
  betas <- data.frame(
    A = c(1, 2, 3),
    B = c(4, 5, 6),
    C = c(7, 8, 9)
  )

  res <- mock_psych(
    correlation_overlap_heatmap(betas)
  )

  df <- res$data

  low <- df[df$tri == "lower", ]

  expect_true(all(low$fill == 1 | low$fill == 0.5))
})
