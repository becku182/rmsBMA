# tests/testthat/test-g_regression.R

test_that("g_regression: g keyword mappings behave as implemented", {
  set.seed(123)

  m <- 60
  x1 <- rnorm(m)
  x2 <- rnorm(m, sd = 2)
  y  <- 2 + x1 + 2 * x2 + rnorm(m)

  dat <- cbind(y, x1, x2)
  r <- ncol(dat) - 1

  # helper to pull betas for comparing implied g choices via closed form
  ref_beta <- function(g_num) {
    Y <- matrix(dat[, 1], ncol = 1)
    X <- dat[, -1, drop = FALSE]
    Z <- cbind(1, X)
    beta_ols <- solve(crossprod(Z), crossprod(Z, Y))
    as.numeric(beta_ols / (1 + g_num))
  }

  # UIP: g = 1/m (also NULL and unknown)
  out_uip  <- g_regression(dat, g = "UIP")
  out_null <- g_regression(dat, g = NULL)
  out_unk  <- g_regression(dat, g = "something_else")

  g_uip <- 1 / m
  expect_equal(as.numeric(out_uip[[1]]),  ref_beta(g_uip), tolerance = 1e-10)
  expect_equal(as.numeric(out_null[[1]]), ref_beta(g_uip), tolerance = 1e-10)
  expect_equal(as.numeric(out_unk[[1]]),  ref_beta(g_uip), tolerance = 1e-10)

  # RIC: g = 1/(r^2)
  g_ric <- 1 / (r^2)
  out_ric <- g_regression(dat, g = "RIC")
  expect_equal(as.numeric(out_ric[[1]]), ref_beta(g_ric), tolerance = 1e-10)

  # Benchmark (note: your code expects exactly "Benchmark" with capital B)
  g_bench <- 1 / max(m, (r^2))
  out_bench <- g_regression(dat, g = "Benchmark")
  expect_equal(as.numeric(out_bench[[1]]), ref_beta(g_bench), tolerance = 1e-10)

  # HQ: g = 1/(log(m)^3)
  g_hq <- 1 / (log(m)^3)
  out_hq <- g_regression(dat, g = "HQ")
  expect_equal(as.numeric(out_hq[[1]]), ref_beta(g_hq), tolerance = 1e-10)

  # rootUIP: g = sqrt(1/m)
  g_root <- sqrt(1 / m)
  out_root <- g_regression(dat, g = "rootUIP")
  expect_equal(as.numeric(out_root[[1]]), ref_beta(g_root), tolerance = 1e-10)
})

test_that("g_regression matches closed-form quantities for a numeric g", {
  set.seed(321)

  m <- 120
  x1 <- rnorm(m)
  x2 <- rnorm(m, sd = 2)
  e  <- rnorm(m, sd = 1)
  y  <- 2 + 1.5 * x1 - 0.7 * x2 + e

  dat <- cbind(y, x1, x2)
  g <- 0.8

  out <- g_regression(dat, g = g)

  betas    <- out[[1]]
  post_se  <- out[[2]]
  log_like <- out[[3]]
  R2       <- out[[4]]
  df       <- out[[5]]
  Dil      <- out[[6]]

  # independent reference calculations
  Y <- matrix(dat[, 1], ncol = 1)
  X <- dat[, -1, drop = FALSE]
  Z <- cbind(1, X)

  r <- ncol(X)
  k <- ncol(Z)
  df_ref <- m - k

  expect_equal(as.numeric(df), df_ref)

  # dilution
  expect_equal(as.numeric(Dil), as.numeric(det(stats::cor(X))), tolerance = 1e-12)

  ZtZ <- crossprod(Z)
  ZtZ_inv <- solve(ZtZ)
  Zty <- crossprod(Z, Y)
  yty <- drop(crossprod(Y))

  beta_ols <- ZtZ_inv %*% Zty
  beta_ref <- beta_ols / (1 + g)
  expect_equal(as.numeric(betas), as.numeric(beta_ref), tolerance = 1e-10)

  yPzy <- yty - drop(t(Zty) %*% ZtZ_inv %*% Zty)
  ybar <- mean(Y)
  SST <- yty - m * (ybar^2)

  EX_1_ref <- (1/(g+1)) * yPzy + (g/(g+1)) * SST

  V <- ZtZ_inv / (1 + g)
  se_2_ref <- EX_1_ref / m
  co_var_ref <- V * ((m * se_2_ref) / (m - 2))
  post_se_ref <- sqrt(diag(co_var_ref))

  expect_equal(as.numeric(post_se), as.numeric(post_se_ref), tolerance = 1e-10)

  log_like_ref <- (r/2) * log(g/(g+1)) - ((m-1)/2) * log(EX_1_ref)
  expect_equal(as.numeric(log_like), as.numeric(log_like_ref), tolerance = 1e-10)

  R2_ref <- 1 - (yPzy / SST)
  expect_equal(as.numeric(R2), as.numeric(R2_ref), tolerance = 1e-10)
})

test_that("g_regression errors on invalid numeric g", {
  set.seed(1)
  dat <- cbind(rnorm(10), rnorm(10), rnorm(10))

  expect_error(g_regression(dat, g = 0), "g must be strictly positive")
  expect_error(g_regression(dat, g = -1), "g must be strictly positive")
  expect_error(g_regression(dat, g = NaN), "g must be strictly positive")
  expect_error(g_regression(dat, g = c(0.1, 0.2)), "g must be strictly positive")
})

test_that("g_regression errors on singular Z'Z", {
  set.seed(42)

  m <- 40
  x1 <- rnorm(m)
  x2 <- 2 * x1   # perfect collinearity
  y  <- rnorm(m)

  dat <- cbind(y, x1, x2)

  expect_error(g_regression(dat, g = 0.5), "Determinant of Z'Z=0")
})
