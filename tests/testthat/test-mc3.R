# Tests for MC^3 model-space sampling.

make_data <- function(seed = 11, m = 50, K = 6) {
  set.seed(seed)
  X <- matrix(stats::rnorm(m * K), m, K)
  y <- 2 + 1.5 * X[, 1] - X[, 2] + stats::rnorm(m, sd = 1)
  d <- cbind(y, X)
  colnames(d) <- c("y", paste0("x", seq_len(K)))
  d
}

# Short chains are used throughout for speed: these tests check structure, not
# convergence, so the convergence warning they legitimately trigger is noise.
quiet_mc3 <- function(expr) suppressWarnings(suppressMessages(expr))

model_keys <- function(M, K) apply(M[, 1:K, drop = FALSE], 1,
                                   function(r) paste0(as.integer(r), collapse = ""))

test_that("draws and burn never switch the method on by themselves", {
  d <- make_data()
  expect_error(model_space(d, M = 3, draws = 100), "only when mc3 = TRUE")
  expect_error(model_space(d, burn = 100),         "only when mc3 = TRUE")
})

test_that("MC3 accepts any admissible M", {
  d <- make_data()
  # Reduced model spaces are supported; the sampler corrects the proposal.
  expect_silent(suppressWarnings(suppressMessages(
    model_space(d, M = 3, mc3 = TRUE, draws = 200, burn = 100))))
  # M = K is fine, and so is omitting M
  expect_silent(suppressWarnings(suppressMessages(
    model_space(d, M = 6, mc3 = TRUE, draws = 200, burn = 100))))
  msgs <- capture_messages(suppressWarnings(
    model_space(d, mc3 = TRUE, draws = 200, burn = 100)))
  expect_match(paste(msgs, collapse = ""), "full model space")
  # M = 0 leaves only the null model and needs no sampling
  expect_error(model_space(d, M = 0, mc3 = TRUE, draws = 100, burn = 50),
               "requires M")
})

test_that("MC3 returns the documented structure", {
  d <- make_data()
  ms <- quiet_mc3(model_space(d, mc3 = TRUE, draws = 500, burn = 200))
  expect_length(ms, 6)
  expect_s3_class(ms, "model_space")
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
  ms_m <- quiet_mc3(model_space(d, mc3 = TRUE, draws = 3000, burn = 1000,
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
  ms_m <- quiet_mc3(model_space(d, mc3 = TRUE, draws = 100000,
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
  ms <- quiet_mc3(model_space(d, mc3 = TRUE, draws = 1000, burn = 500))
  expect_match(paste(capture_messages(bma(ms, round = 12)), collapse = ""),
               "not available for an MC3 model space")
  b <- suppressMessages(bma(ms, round = 12))
  expect_null(b[[3]])
  # everything else still present
  expect_true(is.matrix(b[[1]]))
  expect_false(is.null(b[[10]]))
  expect_identical(b[[16]]$method, "mc3")
})

test_that("dilution under MC3 warns that it is reweighting", {
  d <- make_data()
  ms <- quiet_mc3(model_space(d, mc3 = TRUE, draws = 1000, burn = 500))
  expect_warning(suppressMessages(bma(ms, dilution = 1, dil.Par = 0.5, round = 12)),
                 "reweighting the visited")
})

test_that("enumerated model spaces are unaffected by the MC3 additions", {
  d <- make_data()
  ms <- model_space(d, M = 3, g = "UIP")
  expect_length(ms, 6)                      # info element, method "enumeration"
  expect_identical(ms$info$method, "enumeration")
  b <- bma(ms, round = 12)
  expect_false(is.null(b[[3]]))             # EBA still produced
  expect_null(b[[16]])                      # no chain diagnostics
})

test_that("chain diagnostics are reported, not just stored", {
  d <- make_data()
  msgs <- paste(capture_messages(suppressWarnings(
    model_space(d, mc3 = TRUE, draws = 1000, burn = 500))), collapse = "")
  expect_match(msgs, "cor\\(analytic PMP, visit frequency\\)")
  expect_match(msgs, "acceptance")
})

test_that("a short chain warns about non-convergence, a long one does not", {
  skip_on_cran()
  # Diffuse signal over ten regressors: too few draws to settle, so the
  # analytic and frequency estimates disagree and the diagnostic says so.
  set.seed(1)
  m <- 45; K <- 10
  X <- matrix(stats::rnorm(m * K), m, K)
  y <- 1 + X[, 1] + stats::rnorm(m, sd = 2)
  d <- cbind(y, X); colnames(d) <- c("y", paste0("x", seq_len(K)))

  set.seed(101)
  expect_warning(suppressMessages(model_space(d, mc3 = TRUE, draws = 300, burn = 30)),
                 "may not have converged")

  # The same chain run long enough must stop warning, otherwise the threshold
  # is firing on everything and tells the user nothing.
  set.seed(101)
  expect_no_warning(suppressMessages(
    model_space(d, mc3 = TRUE, draws = 100000, burn = 10000)))
})

test_that("slot 3 is not labeled an EBA table when it is NULL", {
  d <- make_data()
  ms <- quiet_mc3(model_space(d, mc3 = TRUE, draws = 1000, burn = 500))
  b  <- suppressMessages(bma(ms, round = 12))
  expect_null(b[[3]])
  expect_match(names(b)[3], "not available")
})

test_that("constrained MC3 (M < K) only visits admissible models", {
  d <- make_data(); K <- 6; M <- 3
  ms <- suppressWarnings(suppressMessages(
    model_space(d, M = M, mc3 = TRUE, draws = 5000, burn = 1000)))
  sizes <- rowSums(ms[[2]][, 1:K, drop = FALSE])
  expect_true(all(sizes <= M))
  expect_equal(ms[[4]], M)
})

test_that("constrained MC3 reproduces the exact constrained posterior", {
  skip_on_cran()
  # Diffuse signal so posterior mass sits on both sides of the M boundary;
  # this is where the |nbd(g)|/|nbd(g')| correction actually matters.
  set.seed(4)
  m <- 60; K <- 7; M <- 3
  X <- matrix(stats::rnorm(m * K), m, K)
  y <- 0.45*X[,1] + 0.35*X[,2] + 0.25*X[,3] + stats::rnorm(m, sd = 1.6)
  d <- cbind(y, X); colnames(d) <- c("y", paste0("x", seq_len(K)))

  b_e <- bma(model_space(d, M = M, g = "UIP"), EMS = K/2, round = 12)
  set.seed(21)
  ms  <- suppressWarnings(suppressMessages(
    model_space(d, M = M, mc3 = TRUE, draws = 200000, burn = 20000, g = "UIP")))
  b_m <- suppressMessages(bma(ms, EMS = K/2, round = 12))

  ex <- data.frame(k = model_keys(b_e[[10]], K), pmp = b_e[[10]][, K+1],
                   stringsAsFactors = FALSE)
  mc <- data.frame(k = model_keys(b_m[[10]], K),
                   freq = ms[[6]]$visits / sum(ms[[6]]$visits),
                   stringsAsFactors = FALSE)
  j <- merge(ex, mc, by = "k")

  expect_lt(0.5 * sum(abs(j$pmp - j$freq)), 0.03)   # total variation

  # The mass at the boundary is what the correction fixes: without it this
  # ratio collapses towards M/K rather than sitting at 1.
  r <- vapply(j$k, function(s) sum(as.integer(strsplit(s, "")[[1]])), 0)
  expect_equal(sum(j$freq[r == M]) / sum(j$pmp[r == M]), 1, tolerance = 0.05)

  expect_equal(unname(b_e[[1]][-1, "PIP"]), unname(b_m[[1]][-1, "PIP"]),
               tolerance = 2e-2)
})

test_that("a barely-moving chain is flagged", {
  # Strong signal with the mode exactly at the boundary: the chain sticks.
  set.seed(11); m <- 60; K <- 8
  X <- matrix(stats::rnorm(m*K), m, K)
  y <- 2 + 1.5*X[,1] - X[,2] + 0.7*X[,3] + stats::rnorm(m, sd = 1)
  d <- cbind(y, X); colnames(d) <- c("y", paste0("x", seq_len(K)))
  expect_warning(suppressMessages(
    model_space(d, M = 3, mc3 = TRUE, draws = 5000, burn = 1000, g = "UIP")),
    "barely moved")
})

test_that("S3 methods are available on both objects", {
  d <- make_data()
  ms <- model_space(d, M = 3, g = "UIP")
  b  <- bma(ms, EMS = 3, round = 6)

  expect_s3_class(ms, "model_space")
  expect_s3_class(b, "bma")

  # print methods return their argument invisibly
  expect_output(print(ms), "model space")
  expect_output(print(b),  "Bayesian model averaging")
  # capture.output() keeps print()'s own output out of the test log while the
  # visibility of its return value is checked.
  invisible(capture.output(vis <- withVisible(print(ms))$visible))
  expect_false(vis)

  sm <- summary(ms)
  expect_s3_class(sm, "summary.model_space")
  expect_equal(sum(sm$sizes$models), ms$MS)   # every model counted once
  expect_output(print(sm), "models by size")

  sb <- summary(b)
  expect_s3_class(sb, "summary.bma")
  expect_output(print(sb), "posterior inclusion probability")
  # ordered by PIP, intercept aside
  pip <- sb$table[-nrow(sb$table), "PIP"]
  expect_false(is.unsorted(rev(pip[!is.na(pip)])))

  cf <- coef(b)
  expect_type(cf, "double")
  expect_length(cf, ms$K + 1)
  expect_identical(names(cf)[1], "CONST")
  expect_equal(unname(cf), unname(b[[1]][, "PM"]))
  expect_equal(unname(coef(b, conditional = TRUE)), unname(b[[1]][, "PMcon"]))
  expect_false(identical(coef(b, "uniform"), coef(b, "random")))
  expect_error(coef(b, prior = "nonsense"))
})

test_that("positional access into both objects still works", {
  d <- make_data()
  ms <- model_space(d, M = 3, g = "UIP")
  expect_true(is.matrix(ms[[2]]))
  expect_equal(ms[[5]], 6)
  b <- bma(ms, EMS = 3, round = 6)
  expect_true(is.matrix(b[[1]]))
  expect_equal(b[[6]], 6)
})
