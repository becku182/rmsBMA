# A perfectly fitting model has a residual sum of squares of zero in exact
# arithmetic. Whether it comes out as exactly zero depends on the BLAS, which
# is what made CRAN's tests-BLIS machine fail test-fast_ols_const.R. Every
# estimator now snaps rounding-level sums of squares to zero, so the degenerate
# case is the same everywhere. These tests hold that contract for all of them.

const_inputs <- function(m = 15, seed = 2) {
  set.seed(seed)
  list(x = cbind(stats::rnorm(m), stats::rnorm(m)), m = m)
}

test_that("every estimator is deterministic for a constant y", {
  d <- const_inputs()
  for (val in c(5, 0, -3, 1e6)) {
    y <- rep(val, d$m)
    fits <- list(
      fast_ols                = fast_ols(y, d$x),
      fast_ols_HC             = fast_ols_HC(y, d$x),
      fast_ols_const          = fast_ols_const(y),
      ols                     = ols(y, d$x, const = 1),
      g_regression            = g_regression(cbind(y, d$x), g = 0.5),
      g_regression_fast       = g_regression_fast(cbind(y, d$x), g = 0.5),
      g_regression_fast_HC    = g_regression_fast_HC(cbind(y, d$x), g = 0.5),
      g_regression_fast_const = g_regression_fast_const(y, g = 0.5))
    for (nm in names(fits)) {
      info <- paste(nm, "y =", val)
      # never NaN, which is what log() of a negative would give
      expect_false(is.nan(fits[[nm]][[3]]), info = info)
      expect_true(is.infinite(fits[[nm]][[3]]), info = info)
      expect_true(is.nan(fits[[nm]][[4]]), info = info)
    }
  }
})

test_that("a near-constant y cannot produce a NaN log likelihood", {
  # y_center_ss is built as yty - m * mean(y)^2, a difference of two nearly
  # equal numbers. Cancellation could drive it negative and send log() to NaN.
  d <- const_inputs()
  set.seed(1)
  y <- 5 + stats::rnorm(d$m, sd = 1e-9)
  for (fn in c("g_regression", "g_regression_fast", "g_regression_fast_HC")) {
    out <- do.call(fn, list(cbind(y, d$x), g = 0.5))
    expect_false(is.nan(out[[3]]), info = fn)
  }
})

test_that("the guard is inert on ordinary data", {
  # The tolerance scales with the data, so a genuine sum of squares can never
  # reach it: every estimator must still return finite values.
  set.seed(9)
  m <- 40
  x <- cbind(stats::rnorm(m), stats::rnorm(m))
  y <- 2 + x[, 1] - x[, 2] + stats::rnorm(m)
  fits <- list(
    fast_ols             = fast_ols(y, x),
    fast_ols_HC          = fast_ols_HC(y, x),
    fast_ols_const       = fast_ols_const(y),
    ols                  = ols(y, x, const = 1),
    g_regression         = g_regression(cbind(y, x), g = 0.5),
    g_regression_fast    = g_regression_fast(cbind(y, x), g = 0.5),
    g_regression_fast_HC = g_regression_fast_HC(cbind(y, x), g = 0.5))
  for (nm in names(fits)) {
    expect_true(is.finite(fits[[nm]][[3]]), info = nm)
    expect_true(is.finite(fits[[nm]][[4]]), info = nm)
  }
})

test_that("snap_zero only fires below the rounding level of its data", {
  # 300 is sum(y^2) for y = rep(5, 12); tolerance is eps * m * 300 ~ 8e-13
  expect_identical(rmsBMA:::snap_zero(9.5e-30, 300, 12), 0)  # the BLIS residual
  expect_identical(rmsBMA:::snap_zero(-1e-13, 300, 12), 0)   # cancellation
  expect_identical(rmsBMA:::snap_zero(0, 300, 12), 0)
  expect_equal(rmsBMA:::snap_zero(1e-6, 300, 12), 1e-6)      # a real residual
  expect_equal(rmsBMA:::snap_zero(42, 300, 12), 42)
  expect_true(is.na(rmsBMA:::snap_zero(NA_real_, 300, 12)))  # not swallowed
})
