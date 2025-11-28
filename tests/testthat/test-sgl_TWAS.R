
# Mock sparsegl::cv.sparsegl() and coef()

fake_cv_sparsegl <- function(X, y, group, family, pred.loss,
                             offset, nfolds, intercept) {

  # Return object mimicking cv.sparsegl
  structure(
    list(
      lambda = c(0.5, 0.2, 0.1),
      lambda.1se = 0.2,
      fit = list(dummy = TRUE)
    ),
    class = "cv.sparsegl"
  )
}

fake_coef <- function(object, s) {
  # Return a dgCMatrix like sparse matrix structure
  out <- matrix(c(0.1, -0.2, 0), ncol = 1)
  rownames(out) <- c("gene1", "gene2", "(Intercept)")
  return(out)
}

# Use with_mock to patch sparsegl functions
mock_sgl <- function(expr) {
  with_mock(
    `sparsegl::cv.sparsegl` = fake_cv_sparsegl,
    `coef` = fake_coef,
    expr
  )
}

# -------------------------------------------------------------------
# Begin test suite
# -------------------------------------------------------------------

test_that("sgl_TWAS validates input matrix X.", {
  expect_error(sgl_TWAS(X = "not matrix", y = 1:10))
  expect_error(sgl_TWAS(X = matrix(letters, 5, 2), y = rnorm(5)))
})

test_that("sgl_TWAS validates y.", {
  X <- matrix(rnorm(20), nrow = 10)
  expect_error(sgl_TWAS(X, y = c("a","b"), family_func = "gaussian"))
  expect_error(sgl_TWAS(X, y = rnorm(5))) # y length mismatch
})

test_that("sgl_TWAS validates grouping.", {
  X <- matrix(rnorm(20), nrow = 10, ncol = 2)
  y <- rnorm(10)

  expect_error(sgl_TWAS(X, y, grouping = "not integer"))
  expect_error(sgl_TWAS(X, y, grouping = c(1,2,3)))
  expect_error(sgl_TWAS(X, y, grouping = c(1,1), family_func = "gaussian"))
})

test_that("sgl_TWAS validates family and pred_loss.", {
  X <- matrix(rnorm(20), nrow = 10, ncol = 2)
  y <- rnorm(10)

  expect_error(sgl_TWAS(X, y, family_func = "wrong"))
  expect_error(sgl_TWAS(X, y, pred_loss = "badloss"))
})

test_that("sgl_TWAS validates offsets.", {
  X <- matrix(rnorm(20), nrow = 10, ncol = 2)
  y <- rnorm(10)

  expect_error(sgl_TWAS(X, y, offsets = "bad"))
  expect_error(sgl_TWAS(X, y, offsets = 1:5))
})

test_that("sgl_TWAS validates nfolds.", {
  X <- matrix(rnorm(20), nrow = 10, ncol = 2)
  y <- rnorm(10)

  expect_error(sgl_TWAS(X, y, nfolds = 2))
  expect_error(sgl_TWAS(X, y, nfolds = "five"))
})

# Successful mocked model fit

test_that("sgl_TWAS runs successfully.", {
  X <- matrix(rnorm(50), nrow = 10, ncol = 5)
  y <- rnorm(10)
  grouping <- rep(1:2, length.out = 5)

  res <- mock_sgl(
    sgl_TWAS(X, y, grouping = grouping, nfolds = 5)
  )

  expect_true(is.list(res))
  expect_true("sgl_fit" %in% names(res))
  expect_true("lambda1se_coef" %in% names(res))

  expect_true(inherits(res$sgl_fit, "cv.sparsegl"))
  expect_true(is.matrix(res$lambda1se_coef))
})

test_that("sgl_TWAS handles binomial family.", {
  X <- matrix(rnorm(50), nrow = 10, ncol = 5)
  y <- rbinom(10, 1, 0.5)

  res <- mock_sgl(
    sgl_TWAS(X, y, family_func = "binomial", nfolds = 5)
  )

  expect_true(inherits(res$sgl_fit, "cv.sparsegl"))
  expect_true(is.matrix(res$lambda1se_coef))
})

# grouping NULL -> lasso

test_that("sgl_TWAS works when grouping is NULL (lasso).", {
  X <- matrix(rnorm(50), nrow = 10, ncol = 5)
  y <- rnorm(10)

  res <- mock_sgl(
    sgl_TWAS(X, y, grouping = NULL, nfolds = 5)
  )

  expect_true(is.list(res))
  expect_true(is.matrix(res$lambda1se_coef))
})

