test_that("g_regression_fast_HC: betas match beta_ols/(1+g); SEs match manual HC1", {
  set.seed(123)

  m <- 120
  x1 <- rnorm(m)
  x2 <- rnorm(m, sd = 2)

  # heteroskedastic errors so HC SEs differ from classical
  e <- rnorm(m, sd = 0.3 + 0.7 * abs(x1))
  y <- 2 + 1.2 * x1 - 0.8 * x2 + e

  dat <- cbind(y, x1, x2)
  g <- 0.9

  out <- g_regression_fast_HC(dat, g = g)

  # structure
  expect_type(out, "list")
  expect_true(all(c("betas","post_se","log_marg_like","R2","df","Dilution") %in% names(out)))

  betas    <- out$betas
  se_HC1   <- out$post_se
  log_like <- out$log_marg_like
  R2       <- out$R2
  df       <- out$df
  Dil      <- out$Dilution

  # --- independent reference calculations ---
  Y <- matrix(dat[, 1], ncol = 1)
  X <- dat[, -1, drop = FALSE]
  Z <- cbind(1, X)

  r <- ncol(X)
  k <- ncol(Z)
  df_ref <- m - k
  expect_equal(as.numeric(df), df_ref)

  # dilution
  expect_equal(as.numeric(Dil), as.numeric(det(stats::cor(X))), tolerance = 1e-12)

  # OLS
  ZtZ <- crossprod(Z)
  ZtZ_inv <- solve(ZtZ)
  Zty <- crossprod(Z, Y)
  beta_ols <- ZtZ_inv %*% Zty

  # g-prior posterior mean in your code
  beta_ref <- beta_ols / (1 + g)
  expect_equal(as.numeric(betas), as.numeric(beta_ref), tolerance = 1e-10)

  # R2 as in your function: based on yPzy and centered SST
  yty <- drop(crossprod(Y))
  yPzy <- yty - drop(t(Zty) %*% ZtZ_inv %*% Zty)  # SSR from OLS residuals
  ybar <- mean(Y)
  SST <- yty - m * (ybar^2)
  R2_ref <- 1 - (yPzy / SST)
  expect_equal(as.numeric(R2), as.numeric(R2_ref), tolerance = 1e-10)

  # log marginal likelihood as in your function
  EX_1 <- (1/(g+1)) * yPzy + (g/(g+1)) * SST
  log_like_ref <- (r/2) * log(g/(g+1)) - ((m-1)/2) * log(EX_1)
  expect_equal(as.numeric(log_like), as.numeric(log_like_ref), tolerance = 1e-10)

  # HC1 SEs (manual)
  u <- as.numeric(Y - Z %*% beta_ols)
  meat <- crossprod(Z, Z * (u^2))
  vcov_HC1 <- (m / df_ref) * ZtZ_inv %*% meat %*% ZtZ_inv
  se_ref <- sqrt(diag(vcov_HC1))

  expect_equal(as.numeric(se_HC1), as.numeric(se_ref), tolerance = 1e-10)
})

test_that("g_regression_fast_HC errors when df <= 0", {
  # m = 3, k = 1 + 2 = 3 => df = 0
  dat <- matrix(c(
    1, 1, 2,
    2, 2, 4,
    3, 3, 6
  ), ncol = 3, byrow = TRUE) # y, x1, x2

  expect_error(
    g_regression_fast_HC(dat, g = 0.5),
    "Non-positive degrees of freedom"
  )
})

test_that("g_regression_fast_HC errors on singular Z'Z (perfect multicollinearity)", {
  set.seed(42)

  m <- 50
  x1 <- rnorm(m)
  x2 <- 2 * x1
  y  <- rnorm(m)

  dat <- cbind(y, x1, x2)

  expect_error(
    g_regression_fast_HC(dat, g = 0.5),
    "Determinant of Z'Z=0"
  )
})

test_that("g_regression_fast_HC returns expected dimensions", {
  set.seed(1)

  m <- 40
  x <- cbind(rnorm(m), rnorm(m), rnorm(m))
  y <- rnorm(m)
  dat <- cbind(y, x)

  out <- g_regression_fast_HC(dat, g = 1.2)

  # betas: (p+1) x 1
  expect_true(is.matrix(out$betas))
  expect_equal(dim(out$betas), c(ncol(x) + 1, 1))

  # post_se: length p+1
  expect_true(is.numeric(out$post_se))
  expect_length(out$post_se, ncol(x) + 1)

  # scalars
  expect_true(is.numeric(out$log_marg_like)); expect_length(out$log_marg_like, 1)
  expect_true(is.numeric(out$R2));            expect_length(out$R2, 1)
  expect_true(is.numeric(out$df));            expect_length(out$df, 1)
  expect_true(is.numeric(out$Dilution));      expect_length(out$Dilution, 1)
})
