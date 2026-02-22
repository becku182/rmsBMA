# tests/testthat/test-coef_hist.R

test_that("coef_hist errors on invalid weight", {
  set.seed(1)

  n <- 12
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y  <- 1 + 0.6*x1 + rnorm(n, sd = 0.4)

  data <- cbind(y = y, x1 = x1, x2 = x2)

  ms <- model_space(data, M = 2, g = "None", HC = FALSE)
  b  <- bma(ms, EMS = 1, dilution = 0, Narrative = 0, round = 12)

  expect_error(
    coef_hist(b, weight = "wrong"),
    "weight is wrongly specified"
  )
})

test_that("coef_hist returns K+1 ggplots in kernel mode", {
  skip_if_not_installed("ggplot2")

  set.seed(2)

  n <- 12
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y  <- 1 + 0.5*x1 - 0.2*x2 + rnorm(n, sd = 0.4)

  data <- cbind(y = y, x1 = x1, x2 = x2)

  ms <- model_space(data, M = 2, g = "None", HC = FALSE)
  b  <- bma(ms, EMS = 1, dilution = 0, Narrative = 0, round = 12)

  plots <- coef_hist(b, kernel = 1)

  expect_type(plots, "list")
  expect_equal(length(plots), 3) # Const + 2 regressors

  # names should start with "Const" then regressors
  expect_equal(names(plots), c("Const", "x1", "x2"))

  expect_s3_class(plots[[1]], "ggplot")
  expect_s3_class(plots[[2]], "ggplot")
  expect_s3_class(plots[[3]], "ggplot")

  # light sanity: y axis label is Frequency in your function
  expect_equal(plots[[1]]$labels$y, "Frequency")
})

test_that("coef_hist BN=1 requires num of length K+1", {
  skip_if_not_installed("ggplot2")

  set.seed(3)

  n <- 12
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y  <- 1 + rnorm(n)

  data <- cbind(y = y, x1 = x1, x2 = x2)

  ms <- model_space(data, M = 2, g = "None", HC = FALSE)
  b  <- bma(ms, EMS = 1, dilution = 0, Narrative = 0, round = 12)

  # missing num
  expect_error(
    coef_hist(b, BN = 1, num = NULL),
    "Please provide a vector with number of bins"
  )

  # wrong length: should be K+1 = 3
  expect_error(
    coef_hist(b, BN = 1, num = c(10, 10)),
    "num is missspecified"
  )

  # correct length should work
  plots <- coef_hist(b, BN = 1, num = c(10, 10, 10))
  expect_type(plots, "list")
  expect_equal(length(plots), 3)
})

test_that("coef_hist BW='vec' requires binW of length K+1", {
  skip_if_not_installed("ggplot2")

  set.seed(4)

  n <- 12
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y  <- 1 + 0.3*x1 + rnorm(n)

  data <- cbind(y = y, x1 = x1, x2 = x2)

  ms <- model_space(data, M = 2, g = "None", HC = FALSE)
  b  <- bma(ms, EMS = 1, dilution = 0, Narrative = 0, round = 12)

  # missing binW
  expect_error(
    coef_hist(b, BW = "vec", BN = 0, binW = NULL),
    "Please provide a vector with bin width sizes"
  )

  # wrong length: should be K+1 = 3
  expect_error(
    coef_hist(b, BW = "vec", BN = 0, binW = c(0.1, 0.1)),
    "binW is missspecified"
  )

  # correct length should work
  plots <- coef_hist(b, BW = "vec", BN = 0, binW = c(0.2, 0.2, 0.2))
  expect_type(plots, "list")
  expect_equal(length(plots), 3)
})
