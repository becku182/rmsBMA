test_that("fast_ols_const matches lm(y ~ 1) on beta, SE, and df", {
  set.seed(123)

  m <- 60
  y <- 2 + rnorm(m, sd = 0.5)

  out <- fast_ols_const(y)

  betas    <- out[[1]]
  se_B     <- out[[2]]
  log_like <- out[[3]]
  R2       <- out[[4]]
  df       <- out[[5]]
  Dil      <- out[[6]]

  fit <- lm(y ~ 1)

  # beta0 = mean(y)
  expect_equal(as.numeric(betas), as.numeric(coef(fit)), tolerance = 1e-12)

  # df = m - 1
  expect_equal(as.numeric(df), m - 1)

  # SE for intercept in intercept-only model:
  # Var(beta0) = sigma^2 / sum(1^2) = sigma^2 / m
  # sigma^2 estimated with df = m-1
  SSR <- sum(resid(fit)^2)
  sigma2 <- SSR / (m - 1)
  se_ref <- sqrt(sigma2 / m)

  expect_equal(as.numeric(se_B), as.numeric(se_ref), tolerance = 1e-12)

  # dilution is fixed to 1
  expect_equal(as.numeric(Dil), 1)

  # log_like as implemented: (-m/2) * log(SSR)
  log_like_ref <- (-m / 2) * log(SSR)
  expect_equal(as.numeric(log_like), as.numeric(log_like_ref), tolerance = 1e-12)

  # R^2 for intercept-only model is not well-defined in the usual way.
  # In your code, SST is sum((y - mean(y))^2) which equals SSR here,
  # so R2 becomes 0 (up to numerical tolerance), unless SST=0.
  expect_equal(as.numeric(R2), 0, tolerance = 1e-12)
})

test_that("fast_ols_const returns expected structure and types", {
  set.seed(1)

  m <- 10
  y <- rnorm(m)

  out <- fast_ols_const(y)

  expect_type(out, "list")
  expect_length(out, 6)

  # betas: 1x1 matrix
  expect_true(is.matrix(out[[1]]))
  expect_equal(dim(out[[1]]), c(1, 1))

  # se: length 1 numeric
  expect_true(is.numeric(out[[2]]))
  expect_length(out[[2]], 1)

  # scalars are numeric length-1
  expect_true(is.numeric(out[[3]])); expect_length(out[[3]], 1) # log_like
  expect_true(is.numeric(out[[4]])); expect_length(out[[4]], 1) # R2
  expect_true(is.numeric(out[[5]])); expect_length(out[[5]], 1) # df
  expect_true(is.numeric(out[[6]])); expect_length(out[[6]], 1) # dilution
})

test_that("fast_ols_const handles constant y (SST = 0) by returning NaN/Inf R2", {
  # This documents current behavior rather than enforcing a change.
  # If y is constant, SST=0 and SSR=0, so R2 = 1 - 0/0 = NaN.
  y <- rep(5, 12)

  out <- fast_ols_const(y)

  expect_true(is.nan(out[[4]]) || is.infinite(out[[4]]))
  # SSR=0 implies log(SSR) is -Inf, so log_like should be Inf (because -m/2 * -Inf)
  expect_true(is.infinite(out[[3]]))
})

test_that("fast_ols_const is deterministic for a constant y on any BLAS", {
  # CRAN's BLIS machine failed the old form of this test: with a constant y the
  # fit is exact, so the residual sum of squares is zero in exact arithmetic,
  # but a BLAS returning an intercept off by 1e-15 gives a residual around
  # 1e-29 instead, and the log likelihood comes out as a large finite number
  # rather than Inf. The estimator now snaps rounding-level sums of squares to
  # zero, so the answer is the same everywhere.
  for (val in c(5, 0, -3, 1e6)) {
    y <- rep(val, 12)
    out <- fast_ols_const(y)
    expect_true(is.infinite(out[[3]]), info = paste("y =", val))
    expect_true(is.nan(out[[4]]), info = paste("y =", val))
  }

  # ordinary data must be untouched: the tolerance scales with the data, so a
  # genuine residual sum of squares can never reach it
  set.seed(12)
  y <- stats::rnorm(40, mean = 3)
  out <- fast_ols_const(y)
  expect_true(is.finite(out[[3]]))
  expect_true(is.finite(out[[4]]))
})
