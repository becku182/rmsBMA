test_that("fast_ols matches lm() for coefficients, SEs, and R2", {
  set.seed(123)

  n <- 50
  x1 <- rnorm(n)
  x2 <- rnorm(n, sd = 2)
  e  <- rnorm(n, sd = 0.5)
  y  <- 2 + 1 * x1 + 2 * x2 + e
  x  <- cbind(x1, x2)

  out <- fast_ols(y, x)

  # unpack
  betas    <- out[[1]]
  se_B     <- out[[2]]
  log_like <- out[[3]]
  R2       <- out[[4]]
  df       <- out[[5]]
  Dil      <- out[[6]]

  # reference via lm (with intercept)
  dat <- data.frame(y = y, x1 = x1, x2 = x2)
  fit <- lm(y ~ x1 + x2, data = dat)

  # coefficients (including intercept)
  expect_equal(as.numeric(betas), as.numeric(coef(fit)), tolerance = 1e-10)

  # classical (OLS) standard errors
  expect_equal(as.numeric(se_B), as.numeric(summary(fit)$coef[, 2]), tolerance = 1e-10)

  # R^2
  expect_equal(as.numeric(R2), as.numeric(summary(fit)$r.squared), tolerance = 1e-10)

  # df = n - k, where k = p + 1 for intercept
  expect_equal(as.numeric(df), n - (ncol(x) + 1))

  # dilution = det(cor(x))
  expect_equal(as.numeric(Dil), as.numeric(det(stats::cor(x))), tolerance = 1e-12)

  # log_like term as implemented in fast_ols()
  # log_like <- (-r/2) * log(m) + (-m/2) * log(SSR)
  m <- n
  r <- ncol(x) # number of regressors without intercept
  SSR <- sum(resid(fit)^2)
  log_like_ref <- (-r / 2) * log(m) + (-m / 2) * log(SSR)
  expect_equal(as.numeric(log_like), as.numeric(log_like_ref), tolerance = 1e-10)
})

test_that("fast_ols returns expected structure and types", {
  set.seed(1)

  n <- 20
  x <- cbind(rnorm(n), rnorm(n))
  y <- rnorm(n)

  out <- fast_ols(y, x)

  expect_type(out, "list")
  expect_length(out, 6)

  # betas: (p+1) x 1 matrix
  expect_true(is.matrix(out[[1]]))
  expect_equal(dim(out[[1]]), c(ncol(x) + 1, 1))

  # se: length p+1
  expect_true(is.numeric(out[[2]]))
  expect_length(out[[2]], ncol(x) + 1)

  # scalars are numeric length-1
  expect_true(is.numeric(out[[3]])); expect_length(out[[3]], 1) # log_like
  expect_true(is.numeric(out[[4]])); expect_length(out[[4]], 1) # R2
  expect_true(is.numeric(out[[5]])); expect_length(out[[5]], 1) # df
  expect_true(is.numeric(out[[6]])); expect_length(out[[6]], 1) # dilution
})

test_that("fast_ols errors on singular X'X (non-invertible design)", {
  set.seed(42)

  n <- 30
  x1 <- rnorm(n)
  x2 <- 2 * x1  # perfectly collinear
  x  <- cbind(x1, x2)
  y  <- rnorm(n)

  expect_error(fast_ols(y, x))
})
