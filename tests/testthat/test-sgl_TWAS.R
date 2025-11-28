
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


