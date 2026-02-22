test_that("g_regression_fast_const matches closed-form intercept-only g-prior quantities", {
  set.seed(123)

  m <- 100
  y <- 2 + rnorm(m, sd = 3)
  y <- as.matrix(y)

  g <- 0.75
  out <- g_regression_fast_const(y, g = g)

  betas    <- out[[1]]
  post_se  <- out[[2]]
  log_like <- out[[3]]
  R2       <- out[[4]]
  df       <- out[[5]]
  Dil      <- out[[6]]

  # shapes/types
  expect_true(is.matrix(betas))
  expect_equal(dim(betas), c(1, 1))
  expect_true(is.numeric(post_se))
  expect_length(post_se, 1)

  expect_true(is.numeric(log_like)); expect_length(log_like, 1)
  expect_true(is.numeric(R2));       expect_length(R2, 1)
  expect_true(is.numeric(df));       expect_length(df, 1)
  expect_true(is.numeric(Dil));      expect_length(Dil, 1)

  # constants
  expect_equal(as.numeric(Dil), 1)
  expect_equal(as.numeric(df), m - 1)

  # --- closed-form references for intercept-only case ---
  # z = 1_m, Z'Z = m, (Z'Z)^-1 = 1/m, Z'y = sum(y)
  # V = (Z'Z)^-1 / (1+g) = 1 / (m(1+g))
  # beta = V * Z'y = mean(y)/(1+g)
  y_vec <- as.numeric(y)
  ybar <- mean(y_vec)
  beta_ref <- ybar / (1 + g)
  expect_equal(as.numeric(betas), beta_ref, tolerance = 1e-12)

  # yPzy = y'y - (Z'y)'(Z'Z)^-1(Z'y) = sum((y - ybar)^2) = SST
  SST <- sum((y_vec - ybar)^2)

  # In your code, for intercept-only model, EX_1 reduces to SST (independent of g)
  EX_1_ref <- SST

  # posterior SE per your code:
  # se_2 = EX_1/m
  # co_var = V * ((m*se_2)/(m-2)) = V * (EX_1/(m-2))
  # post_se = sqrt(co_var)
  V <- 1 / (m * (1 + g))
  co_var_ref <- V * (EX_1_ref / (m - 2))
  post_se_ref <- sqrt(co_var_ref)
  expect_equal(as.numeric(post_se), post_se_ref, tolerance = 1e-12)

  # log_like per your code: -((m-1)/2) * log(EX_1)
  log_like_ref <- -((m - 1) / 2) * log(EX_1_ref)
  expect_equal(as.numeric(log_like), log_like_ref, tolerance = 1e-12)

  # R2 in your code: SSR=yPzy, SST=y_center_ss; for intercept-only, SSR==SST => R2=0 (unless SST=0)
  expect_equal(as.numeric(R2), 0, tolerance = 1e-12)
})

test_that("g_regression_fast_const handles constant y (SST=0) with NaN/Inf outputs", {
  y <- as.matrix(rep(5, 12))

  out <- g_regression_fast_const(y, g = 0.5)

  # SST = 0 => R2 = 1 - 0/0 => NaN
  expect_true(is.nan(out[[4]]) || is.infinite(out[[4]]))

  # log(EX_1) with EX_1=0 => -Inf; -((m-1)/2)*(-Inf) => +Inf
  expect_true(is.infinite(out[[3]]))
})

test_that("g_regression_fast_const returns expected list structure", {
  set.seed(1)

  y <- as.matrix(rnorm(30))
  out <- g_regression_fast_const(y, g = 0.5)

  expect_type(out, "list")
  expect_length(out, 6)
})
