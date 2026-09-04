# Tests for MC^3 model-space sampling.

make_data <- function(seed = 11, m = 50, K = 6) {
  set.seed(seed)
  X <- matrix(stats::rnorm(m * K), m, K)
  y <- 2 + 1.5 * X[, 1] - X[, 2] + stats::rnorm(m, sd = 1)
  d <- cbind(y, X)
  colnames(d) <- c("y", paste0("x", seq_len(K)))
  d
}

model_keys <- function(M, K) apply(M[, 1:K, drop = FALSE], 1,
                                   function(r) paste0(as.integer(r), collapse = ""))

test_that("draws and burn never switch the method on by themselves", {
  d <- make_data()
  expect_error(model_space(d, M = 3, draws = 100), "only when mc3 = TRUE")
  expect_error(model_space(d, burn = 100),         "only when mc3 = TRUE")
})

test_that("MC3 refuses a constrained model space rather than overriding M", {
  d <- make_data()
  expect_error(model_space(d, M = 3, mc3 = TRUE, draws = 100, burn = 50),
               "full model space")
  # M = K is fine, and so is omitting M
  expect_silent(suppressMessages(
    model_space(d, M = 6, mc3 = TRUE, draws = 200, burn = 100)))
  expect_message(model_space(d, mc3 = TRUE, draws = 200, burn = 100),
                 "full model space")
})

test_that("MC3 returns the documented structure", {
  d <- make_data()
  ms <- suppressMessages(model_space(d, mc3 = TRUE, draws = 500, burn = 200))
  expect_length(ms, 6)
  expect_identical(ms[[6]]$method, "mc3")
  expect_equal(ms[[4]], 6)                 # M = K
  expect_equal(ms[[5]], 6)                 # K
  expect_equal(ms[[3]], nrow(ms[[2]]))     # MS = distinct models visited
  expect_true(ms[[6]]$acceptance > 0 && ms[[6]]$acceptance <= 1)
  expect_equal(sum(ms[[6]]$visits), 500)   # retained draws
  expect_equal(ncol(ms[[2]]), 3 * 6 + 6)
})

test_that("MC3 and enumeration fit identical rows for the same model", {
  d <- make_data(); K <- 6
  ms_e <- model_space(d, M = K, g = "UIP")
  ms_m <- suppressMessages(model_space(d, mc3 = TRUE, draws = 3000, burn = 1000,
                                       g = "UIP"))
  ke <- model_keys(ms_e[[2]], K); km <- model_keys(ms_m[[2]], K)
  common <- intersect(ke, km)
  expect_gt(length(common), 5)
  A <- unname(ms_e[[2]][match(common, ke), , drop = FALSE])
  B <- unname(ms_m[[2]][match(common, km), , drop = FALSE])
  expect_equal(A, B, tolerance = 0)
  expect_identical(colnames(ms_e[[2]]), colnames(ms_m[[2]]))
})

test_that("MC3 visit frequencies converge to the exact posterior", {
  skip_on_cran()
  d <- make_data(); K <- 6
  b_e <- bma(model_space(d, M = K, g = "UIP"), EMS = K/2, round = 12)
  set.seed(99)
  ms_m <- suppressMessages(model_space(d, mc3 = TRUE, draws = 100000,
                                       burn = 10000, g = "UIP"))
  b_m <- suppressMessages(bma(ms_m, EMS = K/2, round = 12))

  ex <- data.frame(k = model_keys(b_e[[10]], K), pmp = b_e[[10]][, K+1],
                   stringsAsFactors = FALSE)
  mc <- data.frame(k = model_keys(b_m[[10]], K), pmp = b_m[[10]][, K+1],
                   freq = ms_m[[6]]$visits / sum(ms_m[[6]]$visits),
                   stringsAsFactors = FALSE)
  j <- merge(ex, mc, by = "k")

  # With EMS = K/2 the binomial model prior is uniform over models, which is
  # what the chain targets, so visit frequencies estimate the exact PMPs.
  expect_gt(stats::cor(j$pmp.x, j$freq), 0.99)
  expect_lt(0.5 * sum(abs(j$pmp.x - j$freq)), 0.05)   # total variation

  # PIPs are what users actually read
  expect_equal(unname(b_e[[1]][-1, "PIP"]), unname(b_m[[1]][-1, "PIP"]),
               tolerance = 1e-2)
})

test_that("EBA is withheld for an MC3 model space", {
  d <- make_data()
  ms <- suppressMessages(model_space(d, mc3 = TRUE, draws = 1000, burn = 500))
  expect_message(bma(ms, round = 12), "not available for an MC3 model space")
  b <- suppressMessages(bma(ms, round = 12))
  expect_null(b[[3]])
  # everything else still present
  expect_true(is.matrix(b[[1]]))
  expect_false(is.null(b[[10]]))
  expect_identical(b[[16]]$method, "mc3")
})

test_that("dilution under MC3 warns that it is reweighting", {
  d <- make_data()
  ms <- suppressMessages(model_space(d, mc3 = TRUE, draws = 1000, burn = 500))
  expect_warning(suppressMessages(bma(ms, dilution = 1, dil.Par = 0.5, round = 12)),
                 "reweighting the visited")
})

test_that("enumerated model spaces are unaffected by the MC3 additions", {
  d <- make_data()
  ms <- model_space(d, M = 3, g = "UIP")
  expect_length(ms, 5)                      # no 6th element
  b <- bma(ms, round = 12)
  expect_false(is.null(b[[3]]))             # EBA still produced
  expect_null(b[[16]])                      # no chain diagnostics
})
