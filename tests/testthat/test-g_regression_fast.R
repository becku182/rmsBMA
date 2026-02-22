test_that("g_regression_fast matches closed-form quantities and lm() R2", {
  set.seed(123)

  m <- 120
  x1 <- rnorm(m)
  x2 <- rnorm(m, sd = 2)
  e  <- rnorm(m, sd = 1)
  y  <- 2 + 1.5 * x1 - 0.7 * x2 + e

  dat <- cbind(y, x1, x2)
  g <- 0.99

  out <- g_regression_fast(dat, g = g)

  betas    <- out[[1]]
  post_se  <- out[[2]]
  log_like <- out[[3]]
  R2       <- out[[4]]
  df       <- out[[5]]
  Dil      <- out[[6]]

  # --- basic shape checks ---
  expect_true(is.matrix(betas))
  expect_true(is.numeric(post_se))
  expect_equal(dim(betas), c(ncol(dat), 1))   # intercept + 2 regressors = 3
  expect_length(post_se, ncol(dat))

  # --- construct reference objects (independent of function internals) ---
  Y <- matrix(dat[, 1], ncol = 1)
  X <- dat[, -1, drop = FALSE]
  Z <- cbind(1, X)

  k <- ncol(Z)
  r <- ncol(X)
  expect_equal(df, m - k)

  # dilution
  expect_equal(as.numeric(Dil), as.numeric(det(stats::cor(X))), tolerance = 1e-12)

  # OLS
  ZtZ <- crossprod(Z)
  Zty <- crossprod(Z, Y)
  ZtZ_inv <- solve(ZtZ)
  beta_ols <- ZtZ_inv %*% Zty

  # posterior mean in your code: beta = (Z'Z)^-1 Z'y / (1+g)
  beta_ref <- beta_ols / (1 + g)
  expect_equal(as.numeric(betas), as.numeric(beta_ref), tolerance = 1e-10)

  # quantities for EX_1 and posterior SE (exactly as your code)
  yty <- drop(crossprod(Y))
  yPzy <- yty - drop(t(Zty) %*% ZtZ_inv %*% Zty)

  y_m <- mean(Y)
  y_center_ss <- yty - m * (y_m^2)

  EX_1_ref <- (1 / (g + 1)) * yPzy + (g / (g + 1)) * y_center_ss

  se_2_ref <- EX_1_ref / m
  V <- ZtZ_inv / (1 + g)
  co_var_ref <- V * ((m * se_2_ref) / (m - 2))
  post_se_ref <- sqrt(diag(co_var_ref))

  expect_equal(as.numeric(post_se), as.numeric(post_se_ref), tolerance = 1e-10)

  # log_like in your code
  log_like_ref <- (r / 2) * log(g / (g + 1)) - ((m - 1) / 2) * log(EX_1_ref)
  expect_equal(as.numeric(log_like), as.numeric(log_like_ref), tolerance = 1e-10)

  # R2 in your code equals OLS R2 with intercept
  fit <- lm(y ~ x1 + x2)
  expect_equal(as.numeric(R2), as.numeric(summary(fit)$r.squared), tolerance = 1e-10)
})

test_that("g_regression_fast returns expected list structure and scalar types", {
  set.seed(1)

  m <- 25
  x <- cbind(rnorm(m), rnorm(m))
  y <- rnorm(m)
  dat <- cbind(y, x)

  out <- g_regression_fast(dat, g = 0.5)

  expect_type(out, "list")
  expect_length(out, 6)

  expect_true(is.matrix(out[[1]]))                 # betas
  expect_true(is.numeric(out[[2]]))                # post_se
  expect_true(is.numeric(out[[3]])); expect_length(out[[3]], 1)  # log_like
  expect_true(is.numeric(out[[4]])); expect_length(out[[4]], 1)  # R2
  expect_true(is.numeric(out[[5]])); expect_length(out[[5]], 1)  # df
  expect_true(is.numeric(out[[6]])); expect_length(out[[6]], 1)  # Dilution
})

test_that("g_regression_fast errors on singular Z'Z", {
  set.seed(42)

  m <- 40
  x1 <- rnorm(m)
  x2 <- 2 * x1          # perfect collinearity => singular Z'Z
  y  <- rnorm(m)

  dat <- cbind(y, x1, x2)

  expect_error(g_regression_fast(dat, g = 0.5), "Determinant of Z'Z=0")
})
