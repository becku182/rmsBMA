test_that("coef_to_full maps intercept and slopes into full-length vector", {
  # K = 4 candidate regressors
  model_row <- c(1, 0, 1, 0)  # select regressors 1 and 3

  # model_coefs: intercept + slopes for selected regressors (in that order)
  model_coefs <- c(10, 0.5, -2)

  out <- coef_to_full(model_coefs, model_row)

  # length = K + 1 (intercept + K slopes)
  expect_true(is.numeric(out))
  expect_length(out, 5)

  # intercept
  expect_equal(out[1], 10)

  # slopes placed in positions (idx + 1)
  # idx = c(1, 3) => positions 2 and 4
  expect_equal(out, c(10, 0.5, 0, -2, 0))
})

test_that("coef_to_full preserves slope order implied by which(model_row==1)", {
  model_row <- c(0, 1, 1, 0, 1) # idx = c(2,3,5)

  # slopes must correspond to regressors 2, 3, 5 in that order
  model_coefs <- c(1.2, 10, 20, 30)

  out <- coef_to_full(model_coefs, model_row)

  expect_equal(out, c(1.2, 0, 10, 20, 0, 30))
})

test_that("coef_to_full works for intercept-only model (no selected regressors)", {
  model_row <- c(0, 0, 0)
  model_coefs <- 3.14

  out <- coef_to_full(model_coefs, model_row)

  expect_equal(out, c(3.14, 0, 0, 0))
})

test_that("coef_to_full coerces inputs (matrix coefs, logical model_row)", {
  model_row <- c(TRUE, FALSE, TRUE)  # coerced to 1,0,1
  model_coefs <- matrix(c(5, -1, 2), ncol = 1) # intercept + 2 slopes

  out <- coef_to_full(model_coefs, model_row)

  # idx = c(1,3) -> positions 2 and 4
  expect_equal(out, c(5, -1, 0, 2))
})

test_that("coef_to_full errors if model_coefs has no intercept", {
  model_row <- c(1, 0, 1)

  expect_error(
    coef_to_full(numeric(0), model_row),
    "must contain at least an intercept"
  )
})

test_that("coef_to_full errors on mismatch between selected regressors and slope count", {
  model_row <- c(1, 0, 1, 1) # k = 3 selected

  # intercept + only 2 slopes -> mismatch
  model_coefs <- c(1, 0.1, 0.2)

  expect_error(
    coef_to_full(model_coefs, model_row),
    "Mismatch: model_row selects 3 regressors"
  )
})
