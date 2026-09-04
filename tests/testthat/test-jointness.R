test_that("jointness returns a square matrix with names and NA diagonal", {
  set.seed(101)

  n <- 25
  x1 <- rnorm(n); x2 <- rnorm(n); x3 <- rnorm(n); x4 <- rnorm(n)
  y  <- 1 + 0.8*x1 - 0.5*x3 + rnorm(n, sd = 0.4)

  data <- cbind(y = y, x1 = x1, x2 = x2, x3 = x3, x4 = x4)
  ms <- model_space(data, M = 4, g = "None", HC = FALSE)
  b  <- bma(ms, EMS = 2, dilution = 0, Narrative = 0, round = 12)

  out <- jointness(b, measure = "PPI", round = 6)

  expect_true(is.matrix(out))
  expect_equal(dim(out), c(4, 4))
  expect_equal(rownames(out), c("x1","x2","x3","x4"))
  expect_equal(colnames(out), c("x1","x2","x3","x4"))

  expect_true(all(is.na(diag(out))))
})

test_that("jointness errors on unknown measure", {
  # minimal fake bma_list to hit the measure branch
  fake <- list(NULL, NULL, NULL, NULL,
               c("x1","x2"),  # [[5]]
               2,             # [[6]] R
               4,             # [[7]] M
               NULL, NULL,
               cbind(matrix(c(0,1,1,1, 1,0,1,0), ncol = 2, byrow = TRUE),  # Z
                     c(0.5,0.5,0.5,0.5),                                 # wU
                     c(0.5,0.5,0.5,0.5)))                                # wR
  expect_error(jointness(fake, measure = "NOPE"), "Unknown measure")
})

test_that("jointness PPI matches joint inclusion probabilities from for_jointness (upper=uniform, lower=random)", {
  set.seed(202)

  n <- 22
  x1 <- rnorm(n); x2 <- rnorm(n); x3 <- rnorm(n); x4 <- rnorm(n)
  y  <- 1 + 0.7*x1 + 0.3*x2 - 0.4*x4 + rnorm(n, sd = 0.5)

  data <- cbind(y = y, x1 = x1, x2 = x2, x3 = x3, x4 = x4)
  ms <- model_space(data, M = 4, g = "None", HC = FALSE)
  b  <- bma(ms, EMS = 2, dilution = 0, Narrative = 0, round = 12)

  # what jointness() uses internally:
  forJ <- b[[10]]
  R <- b[[6]]
  Z  <- as.matrix(forJ[, 1:R])
  wU <- as.numeric(forJ[, R+1])
  wR <- as.numeric(forJ[, R+2])

  # implied PPI matrices:
  Pab_U <- crossprod(Z * wU, Z)  # R x R
  Pab_R <- crossprod(Z * wR, Z)

  out <- jointness(b, measure = "PPI", round = 12)

  # Check a couple of pairs explicitly:
  # Above diagonal => uniform; below diagonal => random
  # (1,2) is above diag, should equal Pab_U[1,2]
  expect_equal(unname(out[1,2]), unname(Pab_U[1,2]), tolerance = 1e-10)
  # (2,1) is below diag, should equal Pab_R[2,1]
  expect_equal(unname(out[2,1]), unname(Pab_R[2,1]), tolerance = 1e-10)

  # Another pair
  expect_equal(unname(out[1,4]), unname(Pab_U[1,4]), tolerance = 1e-10)
  expect_equal(unname(out[4,1]), unname(Pab_R[4,1]), tolerance = 1e-10)

  # Diagonal is forced NA
  expect_true(is.na(out[1,1]))
})

test_that("jointness returns finite values for all supported measures (where defined)", {
  set.seed(303)

  n <- 30
  x1 <- rnorm(n); x2 <- rnorm(n); x3 <- rnorm(n); x4 <- rnorm(n)
  y  <- 1 + 0.6*x1 - 0.2*x2 + 0.4*x3 + rnorm(n, sd = 0.6)

  data <- cbind(y = y, x1 = x1, x2 = x2, x3 = x3, x4 = x4)
  ms <- model_space(data, M = 4, g = "None", HC = FALSE)
  b  <- bma(ms, EMS = 2, dilution = 0, Narrative = 0, round = 12)

  for (measure in c("HCGHM","LS","DW","PPI")) {
    out <- jointness(b, measure = measure, rho = 0.5, round = 6)

    # ignore diagonal NA; require finite elsewhere (may still get Inf/NaN if probabilities hit 0)
    vals <- out[!is.na(out)]
    expect_true(all(is.finite(vals)))
  }
})

test_that("jointness respects rounding", {
  set.seed(404)

  n <- 25
  x1 <- rnorm(n); x2 <- rnorm(n); x3 <- rnorm(n); x4 <- rnorm(n)
  y  <- 1 + 0.9*x1 + rnorm(n, sd = 0.7)

  data <- cbind(y = y, x1 = x1, x2 = x2, x3 = x3, x4 = x4)
  ms <- model_space(data, M = 4, g = "None", HC = FALSE)
  b  <- bma(ms, EMS = 2, dilution = 0, Narrative = 0, round = 12)

  out3 <- jointness(b, measure = "PPI", round = 3)
  out6 <- jointness(b, measure = "PPI", round = 6)

  # Values rounded to 3 decimals should match rounding of the 6-decimal version
  expect_equal(out3, round(out6, 3))
})

test_that("jointness HCGHM changes with rho (sensitivity check)", {
  set.seed(505)

  n <- 25
  x1 <- rnorm(n); x2 <- rnorm(n); x3 <- rnorm(n); x4 <- rnorm(n)
  y  <- 1 + 0.5*x1 - 0.5*x2 + rnorm(n, sd = 0.5)

  data <- cbind(y = y, x1 = x1, x2 = x2, x3 = x3, x4 = x4)
  ms <- model_space(data, M = 4, g = "None", HC = FALSE)
  b  <- bma(ms, EMS = 2, dilution = 0, Narrative = 0, round = 12)

  out_a <- jointness(b, measure = "HCGHM", rho = 0.1, round = 6)
  out_b <- jointness(b, measure = "HCGHM", rho = 0.9, round = 6)

  # should differ somewhere off-diagonal
  expect_true(any(abs(out_a - out_b) > 1e-8, na.rm = TRUE))
})
