test_that("eba returns correct bounds on a small hand-checkable example", {
  betas <- rbind(
    c( 1.0,  0.0,  0.0),  # model 1: none
    c( 2.0, -1.0,  0.0),  # model 2: x1
    c( 3.0,  0.0,  2.0),  # model 3: x2
    c( 4.0, -2.0,  3.0)   # model 4: x1 + x2
  )

  VAR <- rbind(
    c(1.00, 0.00, 0.00),
    c(4.00, 1.00, 0.00),
    c(9.00, 0.00, 4.00),
    c(16.0, 9.00, 1.00)
  )

  Reg_ID <- rbind(
    c(0L, 0L),
    c(1L, 0L),
    c(0L, 1L),
    c(1L, 1L)
  )

  out <- eba(betas, VAR, Reg_ID)

  expect_true(is.matrix(out))
  expect_equal(dim(out), c(3, 5))
  expect_equal(colnames(out), c("Lower_bound", "Low", "Mean_coef", "High", "Upper_bound"))

  # ---- intercept ----
  expect_equal(unname(out[1, "Low"]), 1)
  expect_equal(unname(out[1, "High"]), 4)
  expect_equal(unname(out[1, "Mean_coef"]), mean(c(1, 2, 3, 4)))
  expect_equal(unname(out[1, "Lower_bound"]), -1)   # 1 - 2*1
  expect_equal(unname(out[1, "Upper_bound"]), 12)   # 4 + 2*4

  # ---- x1 (models 2 and 4 only) ----
  expect_equal(unname(out[2, "Low"]), -2)
  expect_equal(unname(out[2, "High"]), -1)
  expect_equal(unname(out[2, "Mean_coef"]), mean(c(-1, -2)))
  expect_equal(unname(out[2, "Lower_bound"]), -8)   # -2 - 2*3
  expect_equal(unname(out[2, "Upper_bound"]), 1)    # -1 + 2*1

  # ---- x2 (models 3 and 4 only) ----
  expect_equal(unname(out[3, "Low"]), 2)
  expect_equal(unname(out[3, "High"]), 3)
  expect_equal(unname(out[3, "Mean_coef"]), mean(c(2, 3)))
  expect_equal(unname(out[3, "Lower_bound"]), -2)   # 2 - 2*2
  expect_equal(unname(out[3, "Upper_bound"]), 5)    # 3 + 2*1
})

test_that("eba leaves regressor rows as NA if never included", {
  set.seed(1)
  MS <- 5
  K <- 2

  betas <- matrix(rnorm(MS * (K + 1)), nrow = MS)
  VAR   <- matrix(runif(MS * (K + 1), min = 0.1, max = 2), nrow = MS)

  Reg_ID <- cbind(
    rep(0L, MS),                 # x1 never included
    c(0L, 1L, 0L, 1L, 0L)         # x2 included sometimes
  )

  out <- eba(betas, VAR, Reg_ID)

  expect_true(all(is.na(out[2, ])))       # x1 row
  expect_false(all(is.na(out[3, ])))      # x2 row
})

test_that("eba respects var_tol and drops intercept entries with VAR <= var_tol", {
  betas <- rbind(
    c(1,  0),
    c(2,  5),
    c(3, -1)
  )
  VAR <- rbind(
    c(0,   1),
    c(4,   1),
    c(0,   1)
  )
  Reg_ID <- matrix(c(1L, 1L, 1L), ncol = 1)

  out <- eba(betas, VAR, Reg_ID, var_tol = 0)

  # intercept keeps only row 2 (VAR=4)
  expect_equal(unname(out[1, "Low"]), 2)
  expect_equal(unname(out[1, "High"]), 2)
  expect_equal(unname(out[1, "Mean_coef"]), 2)
  expect_equal(unname(out[1, "Lower_bound"]), -2)  # 2 - 2*2
  expect_equal(unname(out[1, "Upper_bound"]), 6)   # 2 + 2*2
})

test_that("eba errors on invalid inputs", {
  MS <- 4
  K <- 2

  betas <- matrix(rnorm(MS * (K + 1)), nrow = MS)
  VAR   <- matrix(abs(rnorm(MS * (K + 1))), nrow = MS)
  Reg_ID <- matrix(sample(0:1, MS * K, replace = TRUE), nrow = MS)

  expect_error(eba(betas, VAR, Reg_ID, var_tol = -1), "nonnegative")
  expect_error(eba(betas, VAR, Reg_ID, var_tol = NA), "nonnegative")
  expect_error(eba(betas, VAR, Reg_ID, var_tol = c(0, 1)), "single")

  expect_error(eba(betas[-1, , drop = FALSE], VAR, Reg_ID), "nrow\\(betas\\)")
  expect_error(eba(betas, VAR[, -1, drop = FALSE], Reg_ID), "same dimensions")
  expect_error(eba(betas[, -1, drop = FALSE], VAR[, -1, drop = FALSE], Reg_ID), "ncol\\(betas\\)")

  VAR_bad <- VAR
  VAR_bad[1, 1] <- -0.1
  expect_error(eba(betas, VAR_bad, Reg_ID), "negative values")
})

test_that("eba errors if no valid intercept entries exist after var_tol filtering", {
  betas <- matrix(c(1, 0, 2, 0, 3, 0), ncol = 2, byrow = TRUE) # MS=3, K=1
  VAR   <- matrix(c(0, 1, 0, 1, 0, 1), ncol = 2, byrow = TRUE)
  Reg_ID <- matrix(c(0L, 1L, 1L), ncol = 1)

  expect_error(eba(betas, VAR, Reg_ID, var_tol = 0), "No valid intercept entries")
})
