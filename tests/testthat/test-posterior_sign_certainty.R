# tests/testthat/test-posterior_sign_certainty.R

test_that("posterior_sign_certainty returns correct structure", {
  MS <- 3
  K  <- 2

  pmp_u <- c(0.2, 0.3, 0.5)
  pmp_r <- c(0.1, 0.1, 0.8)

  betas <- matrix(c(
    1,  1,  0,
    2,  0, -1,
    3,  2,  2
  ), nrow = MS, byrow = TRUE)

  VAR <- matrix(1, nrow = MS, ncol = K + 1)
  DF <- c(10, 10, 10)

  Reg_ID <- matrix(c(
    1, 0,
    0, 1,
    1, 1
  ), nrow = MS, byrow = TRUE)

  out <- posterior_sign_certainty(pmp_u, pmp_r, betas, VAR, DF, Reg_ID)

  expect_type(out, "list")
  expect_true(all(c("PSC_uniform","PSC_random","PostMean_uniform","PostMean_random") %in% names(out)))

  expect_equal(dim(out$PSC_uniform), c(K + 1, 1))
  expect_equal(dim(out$PSC_random),  c(K + 1, 1))
  expect_equal(dim(out$PostMean_uniform), c(K + 1, 1))
  expect_equal(dim(out$PostMean_random),  c(K + 1, 1))
})

test_that("posterior_sign_certainty matches manual calculations on a small example", {
  # MS=3 models, K=2 slopes
  MS <- 3
  K  <- 2

  w_u <- c(0.2, 0.3, 0.5)
  w_r <- c(0.1, 0.1, 0.8)
  DF  <- c(10, 5, 20)

  # columns: (Intercept, x1, x2)
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

  out <- posterior_sign_certainty(w_u, w_r, betas, VAR, DF, Reg_ID)

  # ----- manual posterior means -----
  mu_u0 <- sum(w_u * betas[, 1])
  mu_r0 <- sum(w_r * betas[, 1])

  # slopes: sum w * beta * inclusion
  B <- betas[, 2:3, drop = FALSE]
  mu_u <- c(mu_u0, colSums((B * Reg_ID) * w_u))
  mu_r <- c(mu_r0, colSums((B * Reg_ID) * w_r))

  expect_equal(unname(out$PostMean_uniform[1,1]), mu_u[1], tolerance = 1e-12)
  expect_equal(unname(out$PostMean_random[1,1]),  mu_r[1], tolerance = 1e-12)
  expect_equal(as.numeric(out$PostMean_uniform[,1]), mu_u, tolerance = 1e-12)
  expect_equal(as.numeric(out$PostMean_random[,1]),  mu_r, tolerance = 1e-12)

  # ----- manual S terms -----
  # intercept uses pt(beta/sqrt(var), df)
  t0 <- betas[, 1] / sqrt(VAR[, 1])
  c0 <- stats::pt(t0, df = DF)
  S_u0 <- sum(w_u * c0)
  S_r0 <- sum(w_r * c0)

  # slopes: excluded contribute 0.5
  C <- matrix(0.5, nrow = MS, ncol = K)

  # included entries use pt(beta/sqrt(var), df_i)
  for (i in 1:MS) {
    for (k in 1:K) {
      if (Reg_ID[i, k] == 1L) {
        t1 <- B[i, k] / sqrt(VAR[i, k + 1])
        C[i, k] <- stats::pt(t1, df = DF[i])
      }
    }
  }

  S_u <- c(S_u0, colSums(C * w_u))
  S_r <- c(S_r0, colSums(C * w_r))

  # PSC: if mu>0 => S ; if mu<0 => 1-S ; if mu==0 => 0.5
  PSC_u <- ifelse(mu_u > 0, S_u, ifelse(mu_u < 0, 1 - S_u, 0.5))
  PSC_r <- ifelse(mu_r > 0, S_r, ifelse(mu_r < 0, 1 - S_r, 0.5))

  expect_equal(as.numeric(out$PSC_uniform[,1]), PSC_u, tolerance = 1e-12)
  expect_equal(as.numeric(out$PSC_random[,1]),  PSC_r, tolerance = 1e-12)
})

test_that("posterior_sign_certainty returns 0.5 when posterior mean is exactly zero", {
  MS <- 2
  K  <- 1

  w_u <- c(0.5, 0.5)
  w_r <- c(0.5, 0.5)
  DF  <- c(10, 10)

  # slope included in both models; choose betas that cancel exactly in mean
  betas <- matrix(c(
    0,  1,
    0, -1
  ), nrow = MS, byrow = TRUE)

  VAR <- matrix(1, nrow = MS, ncol = K + 1)
  Reg_ID <- matrix(c(1L, 1L), ncol = 1)

  out <- posterior_sign_certainty(w_u, w_r, betas, VAR, DF, Reg_ID)

  # intercept mean is 0 too => PSC=0.5
  expect_equal(unname(out$PSC_uniform[1,1]), 0.5)
  expect_equal(unname(out$PSC_uniform[2,1]), 0.5)
})

test_that("posterior_sign_certainty errors when intercept variance is nonpositive", {
  MS <- 2
  K  <- 1

  w <- c(0.5, 0.5)
  DF <- c(10, 10)

  betas <- matrix(c(
    1, 0.2,
    2, 0.1
  ), nrow = MS, byrow = TRUE)

  VAR <- matrix(c(
    0, 1,   # intercept variance = 0 (invalid)
    1, 1
  ), nrow = MS, byrow = TRUE)

  Reg_ID <- matrix(c(1L, 1L), ncol = 1)

  expect_error(
    posterior_sign_certainty(w, w, betas, VAR, DF, Reg_ID),
    "Intercept variance"
  )
})

test_that("posterior_sign_certainty errors on dimension/length mismatch (stopifnot)", {
  MS <- 3
  K  <- 2

  w_u <- rep(1 / MS, MS)
  w_r <- rep(1 / MS, MS)
  DF  <- rep(10, MS)

  betas <- matrix(rnorm(MS * (K + 1)), nrow = MS)
  VAR   <- matrix(abs(rnorm(MS * (K + 1))) + 0.1, nrow = MS)
  Reg_ID <- matrix(sample(0:1, MS * K, replace = TRUE), nrow = MS)

  expect_error(posterior_sign_certainty(w_u[-1], w_r, betas, VAR, DF, Reg_ID))
  expect_error(posterior_sign_certainty(w_u, w_r, betas, VAR, DF[-1], Reg_ID))
  expect_error(posterior_sign_certainty(w_u, w_r, betas[, -1, drop = FALSE], VAR, DF, Reg_ID))
  expect_error(posterior_sign_certainty(w_u, w_r, betas, VAR[, -1, drop = FALSE], DF, Reg_ID))
})
