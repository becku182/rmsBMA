# tests/testthat/test-model_space.R

test_that("model_space returns correctly sized model space and components (OLS, no HC)", {
  set.seed(1)

  m <- 30
  K <- 3
  x1 <- rnorm(m)
  x2 <- rnorm(m)
  x3 <- rnorm(m)
  y  <- 1 + 2*x1 - 1*x2 + rnorm(m, sd = 0.5)

  dat <- cbind(y, x1, x2, x3)
  colnames(dat) <- c("y", "x1", "x2", "x3")

  M <- 2
  MS_exp <- sum(choose(K, 0:M))

  out <- model_space(dat, M = M, g = "None", HC = FALSE)

  expect_type(out, "list")
  expect_length(out, 6)   # 6 since the info element records how the space was built
  expect_identical(out$info$method, "enumeration")

  x_names <- out[[1]]
  ols_results <- out[[2]]
  MS <- out[[3]]
  M_out <- out[[4]]
  K_out <- out[[5]]

  expect_equal(x_names, c("x1", "x2", "x3"))
  expect_equal(K_out, K)
  expect_equal(M_out, M)
  expect_equal(as.numeric(MS), MS_exp)

  expect_true(is.matrix(ols_results))
  expect_equal(nrow(ols_results), MS_exp)
  # columns: K (id matrix) + (2K+6) = 3K+6
  expect_equal(ncol(ols_results), 3*K + 6)
  expect_false(is.null(colnames(ols_results)))
})

test_that("model_space: intercept-only row matches fast_ols_const (OLS, no HC)", {
  set.seed(2)

  m <- 25
  K <- 3
  x1 <- rnorm(m)
  x2 <- rnorm(m)
  x3 <- rnorm(m)
  y  <- 2 + x1 + rnorm(m)

  dat <- cbind(y, x1, x2, x3)
  colnames(dat) <- c("y", "x1", "x2", "x3")

  out <- model_space(dat, M = 2, g = "None", HC = FALSE)
  ols_results <- out[[2]]

  # By construction: row 1 is the no-regressors model
  y_mat <- as.matrix(dat[, 1])
  ref0 <- fast_ols_const(y_mat)

  K <- ncol(dat) - 1

  # intercept coefficient stored at [1, K+1]
  expect_equal(as.numeric(ols_results[1, K + 1]), as.numeric(ref0[[1]]), tolerance = 1e-10)

  # intercept SE stored at [1, 2K+2]
  expect_equal(as.numeric(ols_results[1, 2*K + 2]), as.numeric(ref0[[2]]), tolerance = 1e-10)

  # log_like, R2, DF, Dilut at the end block
  expect_equal(as.numeric(ols_results[1, 3*K + 3]), as.numeric(ref0[[3]]), tolerance = 1e-10)
  expect_equal(as.numeric(ols_results[1, 3*K + 4]), as.numeric(ref0[[4]]), tolerance = 1e-10)
  expect_equal(as.numeric(ols_results[1, 3*K + 5]), as.numeric(ref0[[5]]), tolerance = 1e-10)
  expect_equal(as.numeric(ols_results[1, 3*K + 6]), as.numeric(ref0[[6]]), tolerance = 1e-10)
})

test_that("model_space: a 1-regressor model row matches fast_ols (OLS, no HC)", {
  set.seed(3)

  m <- 40
  x1 <- rnorm(m)
  x2 <- rnorm(m)
  x3 <- rnorm(m)
  y  <- 1 + 1.5*x1 - 0.3*x3 + rnorm(m, sd = 0.7)

  dat <- cbind(y, x1, x2, x3)
  colnames(dat) <- c("y", "x1", "x2", "x3")

  out <- model_space(dat, M = 1, g = "None", HC = FALSE)
  ols_results <- out[[2]]

  K <- 3
  y_mat <- as.matrix(dat[, 1])
  x_mat <- as.matrix(dat[, 2:4])

  # Find row for model (x1 only): inclusion vector = c(1,0,0)
  target <- c(1, 0, 0)
  ridx <- which(apply(ols_results[, 1:K, drop = FALSE], 1, function(rr) all(rr == target)))
  expect_equal(length(ridx), 1)

  x_ms <- x_mat[, 1, drop = FALSE]
  ref <- fast_ols(y_mat, x_ms)

  # Coefs live in columns (K+1):(2K+1) and are "full" via coef_to_full()
  coef_block <- ols_results[ridx, (K+1):(2*K+1)]
  se_block   <- ols_results[ridx, (2*K+2):(3*K+2)]

  expect_equal(
    as.numeric(coef_block),
    as.numeric(coef_to_full(ref[[1]], target)),
    tolerance = 1e-10
  )

  expect_equal(
    as.numeric(se_block),
    as.numeric(coef_to_full(ref[[2]], target)),
    tolerance = 1e-10
  )

  expect_equal(as.numeric(ols_results[ridx, 3*K+3]), as.numeric(ref[[3]]), tolerance = 1e-10)
  expect_equal(as.numeric(ols_results[ridx, 3*K+4]), as.numeric(ref[[4]]), tolerance = 1e-10)
  expect_equal(as.numeric(ols_results[ridx, 3*K+5]), as.numeric(ref[[5]]), tolerance = 1e-10)
  expect_equal(as.numeric(ols_results[ridx, 3*K+6]), as.numeric(ref[[6]]), tolerance = 1e-10)
})

test_that("model_space: HC=TRUE uses robust SEs (spot-check against fast_ols_HC)", {
  set.seed(4)

  m <- 50
  x1 <- rnorm(m)
  x2 <- rnorm(m)
  x3 <- rnorm(m)
  e  <- rnorm(m, sd = 0.2 + 0.8*abs(x1))
  y  <- 1 + x1 + 0.5*x2 + e

  dat <- cbind(y, x1, x2, x3)
  colnames(dat) <- c("y", "x1", "x2", "x3")

  out <- model_space(dat, M = 1, g = "None", HC = TRUE)
  ols_results <- out[[2]]

  K <- 3
  y_mat <- as.matrix(dat[, 1])
  x_mat <- as.matrix(dat[, 2:4])

  # model (x2 only): c(0,1,0)
  target <- c(0, 1, 0)
  ridx <- which(apply(ols_results[, 1:K, drop = FALSE], 1, function(rr) all(rr == target)))
  expect_equal(length(ridx), 1)

  x_ms <- x_mat[, 2, drop = FALSE]
  ref <- fast_ols_HC(y_mat, x_ms)

  se_block <- ols_results[ridx, (2*K+2):(3*K+2)]
  expect_equal(
    as.numeric(se_block),
    as.numeric(coef_to_full(ref[[2]], target)),
    tolerance = 1e-10
  )
})

test_that("model_space: rejects invalid HC argument", {
  set.seed(5)
  dat <- cbind(rnorm(10), rnorm(10), rnorm(10))
  colnames(dat) <- c("y", "x1", "x2")

  expect_error(
    model_space(dat, M = 1, g = "None", HC = "yes"),
    "HC.*single logical"
  )
})

# UPDATED TEST: new behavior is WARNING + M reset to K (not an error)
test_that("model_space: when M > K it warns and sets M = K", {
  set.seed(6)
  dat <- cbind(rnorm(10), rnorm(10), rnorm(10))
  colnames(dat) <- c("y", "x1", "x2")
  # K = 2, user requests M = 3

  expect_warning(
    out <- model_space(dat, M = 3, g = "None", HC = FALSE),
    "M > K: setting M = K"
  )

  expect_type(out, "list")
  expect_length(out, 6)   # 6 since the info element records how the space was built
  expect_identical(out$info$method, "enumeration")

  M_out <- out[[4]]
  K_out <- out[[5]]
  MS_out <- out[[3]]
  ols_results <- out[[2]]

  expect_equal(K_out, 2)
  expect_equal(M_out, K_out)

  # MS should now reflect M = K
  MS_exp <- sum(choose(K_out, 0:K_out))
  expect_equal(as.numeric(MS_out), MS_exp)

  # And the result matrix should have correct dimensions
  expect_true(is.matrix(ols_results))
  expect_equal(nrow(ols_results), MS_exp)
  expect_equal(ncol(ols_results), 3*K_out + 6)
})

test_that("model_space: g-prior branch works and matches g_regression_fast for a 1-regressor model", {
  set.seed(7)

  m <- 60
  x1 <- rnorm(m)
  x2 <- rnorm(m)
  x3 <- rnorm(m)
  y  <- 2 + 1.2*x1 - 0.4*x2 + rnorm(m)

  dat <- cbind(y, x1, x2, x3)
  colnames(dat) <- c("y", "x1", "x2", "x3")

  g_num <- 0.8
  out <- model_space(dat, M = 1, g = g_num, HC = FALSE)
  ols_results <- out[[2]]

  K <- 3
  y_mat <- as.matrix(dat[, 1])
  x_mat <- as.matrix(dat[, 2:4])

  # model (x3 only): c(0,0,1)
  target <- c(0, 0, 1)
  ridx <- which(apply(ols_results[, 1:K, drop = FALSE], 1, function(rr) all(rr == target)))
  expect_equal(length(ridx), 1)

  x_ms <- x_mat[, 3, drop = FALSE]
  ref <- g_regression_fast(cbind(y_mat, x_ms), g = g_num)

  coef_block <- ols_results[ridx, (K+1):(2*K+1)]
  se_block   <- ols_results[ridx, (2*K+2):(3*K+2)]

  expect_equal(
    as.numeric(coef_block),
    as.numeric(coef_to_full(ref[[1]], target)),
    tolerance = 1e-10
  )

  expect_equal(
    as.numeric(se_block),
    as.numeric(coef_to_full(ref[[2]], target)),
    tolerance = 1e-10
  )
})
