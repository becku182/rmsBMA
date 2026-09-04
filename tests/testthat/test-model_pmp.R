# tests/testthat/test-model_pmp.R

test_that("model_pmp returns three plot objects on a small model space", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("tidyr")
  skip_if_not_installed("ggpubr")

  set.seed(123)

  n <- 12
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y  <- 1 + 0.6*x1 - 0.4*x2 + rnorm(n, sd = 0.3)

  data <- cbind(y = y, x1 = x1, x2 = x2)

  ms <- model_space(data, M = 2, g = "None", HC = FALSE)
  b  <- bma(ms, EMS = 1, dilution = 0, Narrative = 0, round = 12)

  out <- model_pmp(b, top = 3)

  expect_type(out, "list")
  expect_length(out, 3)

  # First two are ggplot objects
  expect_s3_class(out[[1]], "ggplot")
  expect_s3_class(out[[2]], "ggplot")

  # Third is ggarrange output (ggpubr). Class can vary by ggpubr version,
  # but it is at least a "gg" object.
  expect_true(inherits(out[[3]], "gg") || inherits(out[[3]], "ggarrange") || inherits(out[[3]], "gtable"))
})

test_that("model_pmp works when top is NULL (defaults to M)", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("tidyr")
  skip_if_not_installed("ggpubr")

  set.seed(1)

  n <- 10
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y  <- 1 + x1 + rnorm(n, sd = 0.5)

  data <- cbind(y = y, x1 = x1, x2 = x2)

  ms <- model_space(data, M = 2, g = "None", HC = FALSE)
  b  <- bma(ms, EMS = 1, dilution = 0, Narrative = 0, round = 12)

  out <- model_pmp(b, top = NULL)

  expect_type(out, "list")
  expect_length(out, 3)
  expect_s3_class(out[[1]], "ggplot")
  expect_s3_class(out[[2]], "ggplot")
})

test_that("model_pmp handles top > M without error", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("tidyr")
  skip_if_not_installed("ggpubr")

  set.seed(2)

  n <- 10
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y  <- 1 + rnorm(n)

  data <- cbind(y = y, x1 = x1, x2 = x2)

  ms <- model_space(data, M = 2, g = "None", HC = FALSE)
  b  <- bma(ms, EMS = 1, dilution = 0, Narrative = 0, round = 12)

  # Here M = MS = 4. Ask for top bigger than 4.
  expect_message(
    out <- model_pmp(b, top = 999),
    "cannot be higher"
  )

  expect_type(out, "list")
  expect_length(out, 3)
  expect_s3_class(out[[1]], "ggplot")
  expect_s3_class(out[[2]], "ggplot")
})

test_that("model_pmp plots carry expected axis labels (sanity check)", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("tidyr")
  skip_if_not_installed("ggpubr")

  set.seed(3)

  n <- 12
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y  <- 1 + 0.3*x1 + rnorm(n, sd = 0.4)

  data <- cbind(y = y, x1 = x1, x2 = x2)

  ms <- model_space(data, M = 2, g = "None", HC = FALSE)
  b  <- bma(ms, EMS = 1, dilution = 0, Narrative = 0, round = 12)

  out <- model_pmp(b, top = 2)

  # ggplot objects have labels stored in $labels
  expect_true(out[[1]]$labels$x %in% c("Model number in the raniking", "Model number in the ranking"))
  expect_equal(out[[1]]$labels$y, "Prior, Posterior")
})
