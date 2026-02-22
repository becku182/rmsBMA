# tests/testthat/test-posterior_dens.R

test_that("posterior_dens errors on invalid prior", {
  set.seed(1)

  n <- 10
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y  <- 1 + 0.5*x1 + rnorm(n)

  data <- cbind(y = y, x1 = x1, x2 = x2)
  ms <- model_space(data, M = 2, g = "None", HC = FALSE)
  b  <- bma(ms, EMS = 1, dilution = 0, Narrative = 0, round = 12)

  expect_error(
    posterior_dens(b, prior = "wrong"),
    "prior is wrongly specified"
  )
})

test_that("posterior_dens returns K+1 ggplots with correct names (binomial prior)", {
  skip_if_not_installed("ggplot2")

  set.seed(2)

  n <- 12
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y  <- 1 + 0.6*x1 - 0.2*x2 + rnorm(n, sd = 0.4)

  data <- cbind(y = y, x1 = x1, x2 = x2)
  ms <- model_space(data, M = 2, g = "None", HC = FALSE)
  b  <- bma(ms, EMS = 1, dilution = 0, Narrative = 0, round = 12)

  plots <- posterior_dens(b, prior = "binomial")

  expect_type(plots, "list")
  expect_equal(length(plots), 3)               # Const + 2 regressors
  expect_equal(names(plots), c("Const","x1","x2"))

  for (p in plots) {
    expect_s3_class(p, "ggplot")
  }

  # very light sanity: y-axis label is "Density"
  expect_equal(plots[[1]]$labels$y, "Density")
})

test_that("posterior_dens works with beta prior and returns same structure", {
  skip_if_not_installed("ggplot2")

  set.seed(3)

  n <- 12
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y  <- 1 + rnorm(n)

  data <- cbind(y = y, x1 = x1, x2 = x2)
  ms <- model_space(data, M = 2, g = "None", HC = FALSE)
  b  <- bma(ms, EMS = 1, dilution = 0, Narrative = 0, round = 12)

  plots <- posterior_dens(b, prior = "beta")

  expect_type(plots, "list")
  expect_equal(length(plots), 3)
  expect_equal(names(plots), c("Const","x1","x2"))
  expect_s3_class(plots[[2]], "ggplot")
})
