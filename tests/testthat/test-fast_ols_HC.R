test_that("fast_ols_HC matches lm() for betas and R2; matches manual HC1 SEs", {
  set.seed(123)

  n <- 80
  x1 <- rnorm(n)
  x2 <- rnorm(n, sd = 2)
  # add heteroskedastic errors to make robust SEs meaningful
  e  <- rnorm(n, sd = 0.2 + 0.5 * abs(x1))
  y  <- 2 + 1 * x1 + 2 * x2 + e
  x  <- cbind(x1, x2)

  out <- fast_ols_HC(y, x)

  betas    <- out[[1]]
  se_HC1   <- out[[2]]
  log_like <- out[[3]]
  R2       <- out[[4]]
  df       <- out[[5]]
  Dil      <- out[[6]]

  # reference: lm with intercept
  dat <- data.frame(y = y, x1 = x1, x2 = x2)
  fit <- lm(y ~ x1 + x2, data = dat)

  # betas equal
  expect_equal(as.numeric(betas), as.numeric(coef(fit)), tolerance = 1e-10)

  # R2 equal
  expect_equal(as.numeric(R2), as.numeric(summary(fit)$r.squared), tolerance = 1e-10)

  # df = n - (r + 1)
  expect_equal(as.numeric(df), n - (ncol(x) + 1))

  # dilution = det(cor(x))
  expect_equal(as.numeric(Dil), as.numeric(det(stats::cor(x))), tolerance = 1e-12)

  # manual HC1 SEs (MacKinnon & White style via (n/df) scaling)
  X <- cbind(1, x)
  u <- as.numeric(resid(fit))
  XtX_inv <- solve(crossprod(X))
  meat <- crossprod(X, X * (u^2))        # X' diag(u^2) X
  V_HC1 <- (n / df) * XtX_inv %*% meat %*% XtX_inv
  se_ref <- sqrt(diag(V_HC1))

  expect_equal(as.numeric(se_HC1), as.numeric(se_ref), tolerance = 1e-10)

  # log_like term as implemented:
  # (-r/2) * log(n) + (-n/2) * log(SSR)
  r <- ncol(x)
  SSR <- sum(u^2)
  log_like_ref <- (-r / 2) * log(n) + (-n / 2) * log(SSR)
  expect_equal(as.numeric(log_like), as.numeric(log_like_ref), tolerance = 1e-10)
})

test_that("fast_ols_HC returns expected structure and types", {
  set.seed(1)

  n <- 25
  x <- cbind(rnorm(n), rnorm(n))
  y <- rnorm(n)

  out <- fast_ols_HC(y, x)

  expect_type(out, "list")
  expect_length(out, 6)

  # betas: (p+1) x 1 matrix
  expect_true(is.matrix(out[[1]]))
  expect_equal(dim(out[[1]]), c(ncol(x) + 1, 1))

  # se: length p+1 numeric
  expect_true(is.numeric(out[[2]]))
  expect_length(out[[2]], ncol(x) + 1)

  # scalars are numeric length-1
  expect_true(is.numeric(out[[3]])); expect_length(out[[3]], 1) # log_like
  expect_true(is.numeric(out[[4]])); expect_length(out[[4]], 1) # R2
  expect_true(is.numeric(out[[5]])); expect_length(out[[5]], 1) # df
  expect_true(is.numeric(out[[6]])); expect_length(out[[6]], 1) # dilution
})

test_that("fast_ols_HC errors on singular design matrix", {
  set.seed(42)

  n <- 30
  x1 <- rnorm(n)
  x2 <- 2 * x1  # perfect collinearity
  x  <- cbind(x1, x2)
  y  <- rnorm(n)

  expect_error(fast_ols_HC(y, x))
})
