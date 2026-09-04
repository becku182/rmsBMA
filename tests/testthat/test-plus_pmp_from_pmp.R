# tests/testthat/test-plus_pmp_from_pmp.R

test_that("plus_pmp_from_pmp returns correct structure", {
  MS <- 3
  K  <- 2

  pmp_u <- c(0.2, 0.3, 0.5)
  pmp_r <- c(0.1, 0.1, 0.8)

  betas <- matrix(c(
    1,  1,  0,
    2,  0, -1,
    3,  2,  2
  ), nrow = MS, byrow = TRUE) # (MS x (K+1)) = 3x3

  VAR <- matrix(1, nrow = MS, ncol = K + 1)
  DF <- c(10, 10, 10)

  Reg_ID <- matrix(c(
    1, 0,
    0, 1,
    1, 1
  ), nrow = MS, byrow = TRUE)

  out <- plus_pmp_from_pmp(pmp_u, pmp_r, betas, VAR, DF, Reg_ID)

  expect_type(out, "list")
  expect_true(all(c("Plus_PMP_uniform", "Plus_PMP_random") %in% names(out)))

  expect_true(is.matrix(out$Plus_PMP_uniform))
  expect_true(is.matrix(out$Plus_PMP_random))

  expect_equal(dim(out$Plus_PMP_uniform), c(K + 1, 1))
  expect_equal(dim(out$Plus_PMP_random),  c(K + 1, 1))
})

test_that("plus_pmp_from_pmp matches manual calculation on a small example", {
  # MS=3 models, K=2 slopes => betas is 3x3: (Intercept, x1, x2)
  MS <- 3
  K  <- 2

  pmp_u <- c(0.2, 0.3, 0.5)
  pmp_r <- c(0.1, 0.1, 0.8)

  DF <- c(10, 5, 20)

  betas <- matrix(c(
    1,  1,  0,   # model 1: x1 included
    2,  0, -1,   # model 2: x2 included
    3,  2,  2    # model 3: x1 and x2 included
  ), nrow = MS, byrow = TRUE)

  VAR <- matrix(c(
    1,  1,  1,
    4,  1,  1,
    9,  1,  4
  ), nrow = MS, byrow = TRUE)

  Reg_ID <- matrix(c(
    1L, 0L,
    0L, 1L,
    1L, 1L
  ), nrow = MS, byrow = TRUE)

  out <- plus_pmp_from_pmp(pmp_u, pmp_r, betas, VAR, DF, Reg_ID)

  # ---- manual reference ----
  # intercept: all models contribute
  t0 <- betas[, 1] / sqrt(VAR[, 1])
  p0 <- stats::pt(t0, df = DF)
  ref_u_0 <- sum(pmp_u * p0)
  ref_r_0 <- sum(pmp_r * p0)

  # x1: models 1 and 3 only
  rows_x1 <- c(1, 3)
  t1 <- betas[rows_x1, 2] / sqrt(VAR[rows_x1, 2])
  p1 <- stats::pt(t1, df = DF[rows_x1])
  ref_u_1 <- sum(pmp_u[rows_x1] * p1)
  ref_r_1 <- sum(pmp_r[rows_x1] * p1)

  # x2: models 2 and 3 only
  rows_x2 <- c(2, 3)
  t2 <- betas[rows_x2, 3] / sqrt(VAR[rows_x2, 3])
  p2 <- stats::pt(t2, df = DF[rows_x2])
  ref_u_2 <- sum(pmp_u[rows_x2] * p2)
  ref_r_2 <- sum(pmp_r[rows_x2] * p2)

  expect_equal(unname(out$Plus_PMP_uniform[1, 1]), ref_u_0, tolerance = 1e-12)
  expect_equal(unname(out$Plus_PMP_random[1, 1]),  ref_r_0, tolerance = 1e-12)

  expect_equal(unname(out$Plus_PMP_uniform[2, 1]), ref_u_1, tolerance = 1e-12)
  expect_equal(unname(out$Plus_PMP_random[2, 1]),  ref_r_1, tolerance = 1e-12)

  expect_equal(unname(out$Plus_PMP_uniform[3, 1]), ref_u_2, tolerance = 1e-12)
  expect_equal(unname(out$Plus_PMP_random[3, 1]),  ref_r_2, tolerance = 1e-12)
})

test_that("plus_pmp_from_pmp ignores excluded regressors and nonpositive variances for slopes", {
  MS <- 2
  K  <- 2

  pmp_u <- c(0.4, 0.6)
  pmp_r <- c(0.5, 0.5)
  DF <- c(10, 10)

  betas <- matrix(c(
    0,  1,  100,
    0, -1, -100
  ), nrow = MS, byrow = TRUE)

  # Make VAR for x2 invalid (0) so ok = FALSE for x2 even if included
  VAR <- matrix(c(
    1, 1, 0,
    1, 1, 0
  ), nrow = MS, byrow = TRUE)

  # Include x2 everywhere, but it should contribute 0 due to VAR=0
  Reg_ID <- matrix(c(
    1L, 1L,
    1L, 1L
  ), nrow = MS, byrow = TRUE)

  out <- plus_pmp_from_pmp(pmp_u, pmp_r, betas, VAR, DF, Reg_ID)

  # x2 should be exactly 0 contribution because V<=0 => ok FALSE
  expect_equal(unname(out$Plus_PMP_uniform[3, 1]), 0)
  expect_equal(unname(out$Plus_PMP_random[3, 1]),  0)

  # x1 is valid and included: should be nonzero
  expect_true(unname(out$Plus_PMP_uniform[2, 1]) != 0)
  expect_true(unname(out$Plus_PMP_random[2, 1])  != 0)
})

test_that("plus_pmp_from_pmp errors on dimension/length mismatches (stopifnot)", {
  MS <- 3
  K  <- 2

  pmp_u <- rep(1 / MS, MS)
  pmp_r <- rep(1 / MS, MS)
  DF <- rep(10, MS)

  betas <- matrix(rnorm(MS * (K + 1)), nrow = MS)
  VAR   <- matrix(abs(rnorm(MS * (K + 1))) + 0.1, nrow = MS)
  Reg_ID <- matrix(sample(0:1, MS * K, replace = TRUE), nrow = MS)

  # wrong length pmp
  expect_error(plus_pmp_from_pmp(pmp_u[-1], pmp_r, betas, VAR, DF, Reg_ID))

  # wrong DF length
  expect_error(plus_pmp_from_pmp(pmp_u, pmp_r, betas, VAR, DF[-1], Reg_ID))

  # wrong betas dimensions
  expect_error(plus_pmp_from_pmp(pmp_u, pmp_r, betas[, -1, drop = FALSE], VAR, DF, Reg_ID))

  # wrong VAR dimensions
  expect_error(plus_pmp_from_pmp(pmp_u, pmp_r, betas, VAR[, -1, drop = FALSE], DF, Reg_ID))
})
