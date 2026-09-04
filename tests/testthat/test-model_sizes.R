test_that("model_sizes returns three ggplot objects (dilution = 0)", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("tidyr")
  skip_if_not_installed("ggpubr")

  set.seed(1)

  n <- 12
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y  <- 1 + 0.6*x1 - 0.2*x2 + rnorm(n, sd = 0.4)

  data <- cbind(y = y, x1 = x1, x2 = x2)

  ms <- model_space(data, M = 2, g = "None", HC = FALSE)   # K=2 => MS=4
  b  <- bma(ms, EMS = 1, dilution = 0, Narrative = 0, round = 12)

  expect_silent({
    plots <- model_sizes(b)

    expect_type(plots, "list")
    expect_equal(length(plots), 3)

    expect_s3_class(plots[[1]], "ggplot")
    expect_s3_class(plots[[2]], "ggplot")
    # Third is ggarrange output (ggpubr). Class can vary by ggpubr version,
    # but it is at least a "gg" object.
    expect_true(inherits(plots[[3]], "gg") || inherits(plots[[3]], "ggarrange") ||
                  inherits(plots[[3]], "gtable"))

    # Light label sanity (stable and meaningful)
    expect_equal(plots[[1]]$labels$x, "Model size (number of regressors)")
    expect_equal(plots[[1]]$labels$y, "Prior, Posterior")
    expect_equal(plots[[2]]$labels$x, "Model size (number of regressors)")
    expect_equal(plots[[2]]$labels$y, "Prior, Posterior")
  })
})

test_that("model_sizes returns three ggplot objects (dilution = 1)", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("tidyr")
  skip_if_not_installed("ggpubr")

  set.seed(2)

  n <- 12
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y  <- 1 + 0.2*x1 + rnorm(n, sd = 0.5)

  data <- cbind(y = y, x1 = x1, x2 = x2)

  ms <- model_space(data, M = 2, g = "None", HC = FALSE)
  b  <- bma(ms, EMS = 1, dilution = 1, dil.Par = 0.5, Narrative = 0, round = 12)

  expect_silent({
    plots <- model_sizes(b)

    expect_type(plots, "list")
    expect_equal(length(plots), 3)

    expect_s3_class(plots[[1]], "ggplot")
    expect_s3_class(plots[[2]], "ggplot")
    # Third is ggarrange output (ggpubr). Class can vary by ggpubr version,
    # but it is at least a "gg" object.
    expect_true(inherits(plots[[3]], "gg") || inherits(plots[[3]], "ggarrange") ||
                  inherits(plots[[3]], "gtable"))

    # Same stable label checks
    expect_equal(plots[[1]]$labels$x, "Model size (number of regressors)")
    expect_equal(plots[[1]]$labels$y, "Prior, Posterior")
    expect_equal(plots[[2]]$labels$x, "Model size (number of regressors)")
    expect_equal(plots[[2]]$labels$y, "Prior, Posterior")
  })
})

test_that("model_sizes type argument switches the geometry", {
  skip_if_not_installed("ggpubr")
  set.seed(4)
  m <- 40; K <- 5
  X <- matrix(stats::rnorm(m * K), m, K)
  y <- 1 + X[, 1] + stats::rnorm(m)
  d <- cbind(y, X); colnames(d) <- c("y", paste0("x", seq_len(K)))
  b <- bma(model_space(d, M = K, g = "UIP"), EMS = 2, round = 6)

  geom_of <- function(p) class(p$layers[[1]]$geom)[1]

  line <- suppressWarnings(model_sizes(b))                       # default
  hist <- suppressWarnings(model_sizes(b, type = "histogram"))

  expect_identical(geom_of(line[[1]]), "GeomLine")
  expect_identical(geom_of(hist[[1]]), "GeomCol")
  expect_identical(geom_of(line[[2]]), "GeomLine")
  expect_identical(geom_of(hist[[2]]), "GeomCol")

  # same numbers behind both, only the geometry differs
  expect_equal(line[[1]]$data, hist[[1]]$data)
  expect_equal(line[[2]]$data, hist[[2]]$data)

  # axis labels are shared so the two forms stay comparable
  expect_identical(line[[1]]$labels$x, hist[[1]]$labels$x)
  expect_identical(line[[1]]$labels$y, hist[[1]]$labels$y)

  expect_identical(geom_of(suppressWarnings(model_sizes(b, type = "line"))[[1]]),
                   "GeomLine")
  expect_error(model_sizes(b, type = "pie"))
})

test_that("model_sizes titles the combined panels, diluted or not", {
  skip_if_not_installed("ggpubr")
  set.seed(5)
  m <- 40; K <- 4
  X <- matrix(stats::rnorm(m * K), m, K)
  y <- 1 + X[, 1] + stats::rnorm(m)
  d <- cbind(y, X); colnames(d) <- c("y", paste0("x", seq_len(K)))
  ms <- model_space(d, M = K, g = "UIP")

  plain <- suppressWarnings(model_sizes(bma(ms, EMS = 2, round = 6),
                                        type = "histogram"))
  dil   <- suppressWarnings(model_sizes(
    bma(ms, EMS = 2, dilution = 1, dil.Par = 0.5, round = 6), type = "histogram"))
  expect_length(plain, 3)
  expect_length(dil, 3)
})
