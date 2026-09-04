testthat::test_that("best_models returns expected structure and ranking (criterion = 1)", {

  R <- 2L
  K <- R + 1L
  MS <- 4L
  x_names <- c("x1", "x2")

  Reg_ID <- rbind(
    c(0, 0),
    c(1, 0),
    c(0, 1),
    c(1, 1)
  )

  betas <- rbind(
    c(1.0,  0.0,  0.0),
    c(1.2,  0.5,  0.0),
    c(0.8,  0.0, -0.3),
    c(1.1,  0.4, -0.2)
  )

  ses <- rbind(
    c(0.10, 0.10, 0.10),
    c(0.10, 0.20, 0.10),
    c(0.10, 0.10, 0.20),
    c(0.10, 0.20, 0.20)
  )

  DF <- c(10, 10, 10, 10)
  R2 <- c(0.10, 0.40, 0.30, 0.50)

  pmp_uniform <- c(0.05, 0.60, 0.20, 0.15)
  pmp_random  <- c(0.10, 0.10, 0.70, 0.10)

  # Matches bma(): cbind(Reg_ID, betas, SE, DF, R2, pmp_uniform, pmp_random)
  for_best_models <- cbind(Reg_ID, betas, ses, DF, R2, pmp_uniform, pmp_random)

  bma_list <- vector("list", 11)
  bma_list[[5]]  <- x_names
  bma_list[[6]]  <- R
  bma_list[[7]]  <- MS
  bma_list[[11]] <- for_best_models

  out <- best_models(bma_list, criterion = 1, best = 3, round = 3, estimate = TRUE)

  testthat::expect_type(out, "list")
  testthat::expect_length(out, 6)

  inclusion <- out[[1]]
  est_tbl   <- out[[2]]

  # Row names may be stored as a (K+2) x 1 matrix, so coerce safely:
  rn_incl <- rownames(inclusion)
  rn_incl <- if (is.matrix(rn_incl)) as.character(rn_incl[, 1]) else as.character(rn_incl)

  rn_est <- rownames(est_tbl)
  rn_est <- if (is.matrix(rn_est)) as.character(rn_est[, 1]) else as.character(rn_est)

  testthat::expect_equal(rn_incl, c("Const", "x1", "x2", "PMP", "R^2"))
  testthat::expect_equal(rn_est,  c("Const", "x1", "x2", "PMP", "R^2"))

  testthat::expect_equal(dim(inclusion), c(K + 2L, 3L))
  testthat::expect_equal(dim(est_tbl),   c(K + 2L, 3L))

  # Ranking check (criterion = 1 uses PMP_uniform)
  # Highest is model 2 -> includes x1=1, x2=0
  testthat::expect_equal(as.numeric(inclusion["x1", 1]), 1)
  testthat::expect_equal(as.numeric(inclusion["x2", 1]), 0)

  # PMP row should be decreasing
  pmp_row <- as.numeric(inclusion["PMP", ])
  testthat::expect_true(all(diff(pmp_row) <= 0))
})

testthat::test_that("best_models uses criterion=2 (PMP_random) for ranking", {

  R <- 2L
  K <- R + 1L
  MS <- 4L
  x_names <- c("x1", "x2")

  Reg_ID <- rbind(
    c(0, 0),
    c(1, 0),
    c(0, 1),
    c(1, 1)
  )

  # Keep simple, but correct dimensions:
  betas <- rbind(
    c(1, 0, 0),
    c(1, 1, 0),
    c(1, 0, 1),
    c(1, 1, 1)
  )
  ses <- matrix(0.1, nrow = MS, ncol = K)

  DF <- rep(10, MS)
  R2 <- c(0.1, 0.2, 0.3, 0.4)

  pmp_uniform <- c(0.25, 0.25, 0.25, 0.25)
  pmp_random  <- c(0.05, 0.10, 0.80, 0.05)  # model 3 best

  for_best_models <- cbind(Reg_ID, betas, ses, DF, R2, pmp_uniform, pmp_random)

  bma_list <- vector("list", 11)
  bma_list[[5]]  <- x_names
  bma_list[[6]]  <- R
  bma_list[[7]]  <- MS
  bma_list[[11]] <- for_best_models

  # Use best >= 2 to avoid stricter row.names<- edge-case with matrix RHS
  out <- best_models(bma_list, criterion = 2, best = 2, round = 3, estimate = TRUE)
  inclusion <- out[[1]]

  # Under PMP_random the best model is model 3: includes x2 only
  testthat::expect_equal(as.numeric(inclusion["x1", 1]), 0)
  testthat::expect_equal(as.numeric(inclusion["x2", 1]), 1)
})

testthat::test_that("best_models clamps best > M to best = M and still returns tables", {

  R <- 2L
  K <- R + 1L
  MS <- 4L
  x_names <- c("x1", "x2")

  Reg_ID <- rbind(c(0,0), c(1,0), c(0,1), c(1,1))
  betas  <- matrix(0, nrow = MS, ncol = K); betas[,1] <- 1
  ses    <- matrix(0.1, nrow = MS, ncol = K)
  DF     <- rep(10, MS)
  R2     <- rep(0.1, MS)
  pmp_uniform <- c(0.4, 0.3, 0.2, 0.1)
  pmp_random  <- c(0.1, 0.2, 0.3, 0.4)

  for_best_models <- cbind(Reg_ID, betas, ses, DF, R2, pmp_uniform, pmp_random)

  bma_list <- vector("list", 11)
  bma_list[[5]]  <- x_names
  bma_list[[6]]  <- R
  bma_list[[7]]  <- MS
  bma_list[[11]] <- for_best_models

  testthat::expect_message(
    out <- best_models(bma_list, criterion = 1, best = 999, round = 3, estimate = TRUE),
    "best > M"
  )

  inclusion <- out[[1]]
  testthat::expect_equal(ncol(inclusion), MS)
})
