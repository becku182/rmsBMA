# tests/testthat/test-bma.R

test_that("bma returns a well-formed object on a tiny model space (K=2)", {
  set.seed(123)

  # --- tiny data: K=2 regressors, M=2 => MS=2^2=4 models
  n <- 12
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y  <- 1 + 0.5 * x1 - 0.25 * x2 + rnorm(n, sd = 0.2)

  data <- cbind(y = y, x1 = x1, x2 = x2)

  ms <- model_space(data, M = 2, g = "None", HC = FALSE)

  out <- bma(ms, EMS = 1, dilution = 0, Narrative = 0, round = 10)

  # --- list shape
  expect_type(out, "list")
  expect_length(out, 15)

  # --- extracted scalars
  expect_equal(out[[6]], 2)          # K
  expect_equal(out[[7]], 4)          # MS
  expect_equal(out[[8]], 1)          # EMS we passed
  expect_equal(out[[9]], 0)          # dilution flag

  # --- names of regressors
  expect_equal(out[[5]], c("x1", "x2"))

  # --- main tables dimensions: (K+1) rows, 7 columns
  uniform_table <- out[[1]]
  random_table  <- out[[2]]

  expect_true(is.matrix(uniform_table))
  expect_true(is.matrix(random_table))
  expect_equal(dim(uniform_table), c(3, 7))
  expect_equal(dim(random_table),  c(3, 7))

  # row/col names present and correct
  expect_equal(rownames(uniform_table), c("CONST", "x1", "x2"))
  expect_equal(colnames(uniform_table), c("PIP","PM","PSD","PMcon","PSDcon","P(+)","PSC"))

  # intercept PIP is NA by construction
  expect_true(is.na(uniform_table["CONST","PIP"]))
  expect_true(is.na(random_table["CONST","PIP"]))

  # PIPs for slopes in [0,1]
  expect_true(all(uniform_table[c("x1","x2"),"PIP"] >= 0))
  expect_true(all(uniform_table[c("x1","x2"),"PIP"] <= 1))
  expect_true(all(random_table[c("x1","x2"),"PIP"] >= 0))
  expect_true(all(random_table[c("x1","x2"),"PIP"] <= 1))

  # PSC and P(+) are probabilities in [0,1]
  expect_true(all(uniform_table[,"PSC"] >= 0 & uniform_table[,"PSC"] <= 1))
  expect_true(all(random_table[,"PSC"]  >= 0 & random_table[,"PSC"]  <= 1))
  expect_true(all(uniform_table[,"P(+)"] >= 0 & uniform_table[,"P(+)"] <= 1))
  expect_true(all(random_table[,"P(+)"]  >= 0 & random_table[,"P(+)"]  <= 1))

  # --- EBA object has K+1 rows and required columns
  eba_object <- out[[3]]
  expect_true(is.data.frame(eba_object))
  expect_equal(nrow(eba_object), 3)
  expect_true(all(c("Variable","Lower bound","Minimum","Mean","Maximum","Upper bound","EBA test","%(+)") %in% names(eba_object)))

  # --- pms_table has 2 rows (Binomial, Binomial-beta) and 2 cols
  pms_table <- out[[4]]
  expect_true(is.matrix(pms_table))
  expect_equal(dim(pms_table), c(2, 2))
  expect_equal(rownames(pms_table), c("Binomial","Binomial-beta"))
  expect_equal(colnames(pms_table), c("Prior model size","Posterior model size"))

  # --- for_model_pmp must have 4 columns: prior_u, prior_r, pmp_u, pmp_r
  for_model_pmp <- out[[12]]
  expect_equal(dim(for_model_pmp), c(4, 4))
  # posterior model probabilities sum to 1
  expect_equal(sum(for_model_pmp[, 3]), 1, tolerance = 1e-10)
  expect_equal(sum(for_model_pmp[, 4]), 1, tolerance = 1e-10)

  # --- for_model_sizes must have M+1 rows (=3) and 4 cols
  for_model_sizes <- out[[13]]
  expect_equal(dim(for_model_sizes), c(3, 4))
  # prior size probs sum to 1
  expect_equal(sum(for_model_sizes[, 1]), 1, tolerance = 1e-10)
  expect_equal(sum(for_model_sizes[, 2]), 1, tolerance = 1e-10)
  # posterior size probs sum to 1
  expect_equal(sum(for_model_sizes[, 3]), 1, tolerance = 1e-10)
  expect_equal(sum(for_model_sizes[, 4]), 1, tolerance = 1e-10)

  # --- for_jointness has MS rows and (K + 2) cols: Reg_ID + pmp_u + pmp_r
  for_jointness <- out[[10]]
  expect_equal(dim(for_jointness), c(4, 4))

  # --- for_best_models has MS rows and (K + (K+1) + (K+1) + 2 + 2) cols
  # = K Reg_ID + betas(K+1) + SE(K+1) + DF + R2 + pmp_u + pmp_r
  # = K + (K+1) + (K+1) + 4 = 2 + 3 + 3 + 4 = 12
  for_best_models <- out[[11]]
  expect_equal(nrow(for_best_models), 4)
  expect_equal(ncol(for_best_models), 12)

  # --- alphas length is MS (stored as MS x 1 matrix)
  alphas <- out[[14]]
  expect_true(is.matrix(alphas))
  expect_equal(dim(alphas), c(4, 1))
})

test_that("bma Narrative prior downweights models with multiple substitutes (tiny K=2)", {
  set.seed(1)

  n <- 10
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y  <- 1 + x1 + rnorm(n, sd = 0.5)
  data <- cbind(y = y, x1 = x1, x2 = x2)

  ms <- model_space(data, M = 2, g = "None", HC = FALSE)

  # Nar_vec groups both regressors together => group 1 contains x1 and x2
  Nar_vec <- c(1L, 1L)

  # Without narrative (baseline)
  out0 <- bma(ms, EMS = 1, Narrative = 0, dilution = 0, round = 12)

  # With narrative, scalar p < 1 (penalize having both x1 and x2 together)
  out1 <- bma(ms, EMS = 1, Narrative = 1, p = 0.5, Nar_vec = Nar_vec, dilution = 0, round = 12)

  pmp0 <- out0[[12]][, 3]  # pmp_uniform
  pmp1 <- out1[[12]][, 3]

  # Model with both regressors included is the one with rowSums(Reg_ID)==2.
  Reg_ID <- ms[[2]][, 1:2, drop = FALSE]
  idx_both <- which(rowSums(Reg_ID) == 2)

  # Narrative prior should not increase probability of the "both included" model
  # (it should be penalized relative to other models)
  expect_true(pmp1[idx_both] <= pmp0[idx_both] + 1e-12)

  # Posterior probabilities still sum to 1
  expect_equal(sum(pmp1), 1, tolerance = 1e-10)
})

test_that("bma errors when both dilution and Narrative are enabled", {
  set.seed(2)

  n <- 8
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y  <- 1 + rnorm(n)
  data <- cbind(y = y, x1 = x1, x2 = x2)

  ms <- model_space(data, M = 2, g = "None", HC = FALSE)

  expect_error(
    bma(ms, dilution = 1, Narrative = 1, Nar_vec = c(1L,1L), p = 0.7),
    "CANNOT CHOOSE BOTH"
  )
})

test_that("bma errors when Narrative=1 but Nar_vec missing or wrong length", {
  set.seed(3)

  n <- 8
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y  <- 1 + rnorm(n)
  data <- cbind(y = y, x1 = x1, x2 = x2)

  ms <- model_space(data, M = 2, g = "None", HC = FALSE)

  expect_error(
    bma(ms, Narrative = 1, p = 0.7, Nar_vec = NULL),
    "Nar_vec"
  )

  expect_error(
    bma(ms, Narrative = 1, p = 0.7, Nar_vec = c(1L)),
    "Nar_vec should have K elements"
  )
})

test_that("bma errors when Narrative=1 and p is a vector with wrong length vs groups", {
  set.seed(4)

  n <- 10
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y  <- 1 + rnorm(n)
  data <- cbind(y = y, x1 = x1, x2 = x2)

  ms <- model_space(data, M = 2, g = "None", HC = FALSE)

  # Two groups implied (1 and 2), but p length mismatch
  Nar_vec <- c(1L, 2L)

  expect_error(
    bma(ms, Narrative = 1, Nar_vec = Nar_vec, p = c(0.7, 0.5, 0.2)),
    "p is misspecified"
  )
})


test_that("bma works and is internally consistent on MS=16 (K=4, full model space)", {
  set.seed(202)

  n <- 25
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x3 <- rnorm(n)
  x4 <- rnorm(n)

  # build y with some signal
  y <- 1 + 0.7*x1 - 0.5*x3 + rnorm(n, sd = 0.4)

  data <- cbind(y = y, x1 = x1, x2 = x2, x3 = x3, x4 = x4)

  # full model space: K=4, M=4 => MS=16
  ms <- model_space(data, M = 4, g = "None", HC = FALSE)

  out <- bma(ms, EMS = 2, dilution = 0, Narrative = 0, round = 12)

  # --- basic extracted scalars
  expect_equal(out[[6]], 4)    # K
  expect_equal(out[[7]], 16)   # MS
  expect_equal(out[[8]], 2)    # EMS

  # --- Tables are (K+1) x 7
  uniform_table <- out[[1]]
  random_table  <- out[[2]]

  expect_equal(dim(uniform_table), c(5, 7))
  expect_equal(dim(random_table),  c(5, 7))
  expect_equal(rownames(uniform_table), c("CONST","x1","x2","x3","x4"))
  expect_equal(colnames(uniform_table), c("PIP","PM","PSD","PMcon","PSDcon","P(+)","PSC"))

  # --- Posterior model probabilities sum to 1
  for_model_pmp <- out[[12]]
  expect_equal(dim(for_model_pmp), c(16, 4))
  expect_equal(sum(for_model_pmp[, 3]), 1, tolerance = 1e-10) # pmp_uniform
  expect_equal(sum(for_model_pmp[, 4]), 1, tolerance = 1e-10) # pmp_random

  # --- Posterior model size probabilities sum to 1 and have length M+1 = 5
  for_model_sizes <- out[[13]]
  expect_equal(dim(for_model_sizes), c(5, 4))
  expect_equal(sum(for_model_sizes[, 1]), 1, tolerance = 1e-10) # prior uniform sizes
  expect_equal(sum(for_model_sizes[, 2]), 1, tolerance = 1e-10) # prior random sizes
  expect_equal(sum(for_model_sizes[, 3]), 1, tolerance = 1e-10) # post uniform sizes
  expect_equal(sum(for_model_sizes[, 4]), 1, tolerance = 1e-10) # post random sizes

  # --- PIP range checks (slopes only)
  expect_true(all(uniform_table[c("x1","x2","x3","x4"), "PIP"] >= 0))
  expect_true(all(uniform_table[c("x1","x2","x3","x4"), "PIP"] <= 1))
  expect_true(all(random_table[c("x1","x2","x3","x4"), "PIP"] >= 0))
  expect_true(all(random_table[c("x1","x2","x3","x4"), "PIP"] <= 1))

  # --- PSC and P(+) are probabilities
  expect_true(all(uniform_table[, "PSC"] >= 0 & uniform_table[, "PSC"] <= 1))
  expect_true(all(random_table[, "PSC"]  >= 0 & random_table[, "PSC"]  <= 1))
  expect_true(all(uniform_table[, "P(+)"] >= 0 & uniform_table[, "P(+)"] <= 1))
  expect_true(all(random_table[, "P(+)"]  >= 0 & random_table[, "P(+)"]  <= 1))

  # --- Jointness table dimensions: MS rows, K + 2 cols
  for_jointness <- out[[10]]
  expect_equal(dim(for_jointness), c(16, 6))

  # --- for_best_models dimensions: MS rows, K + (K+1) + (K+1) + 4 cols = 4 + 5 + 5 + 4 = 18
  for_best_models <- out[[11]]
  expect_equal(nrow(for_best_models), 16)
  expect_equal(ncol(for_best_models), 18)

  # --- EBA object has K+1 rows and required columns
  eba_object <- out[[3]]
  expect_true(is.data.frame(eba_object))
  expect_equal(nrow(eba_object), 5)
  expect_true(all(c("Variable","Lower bound","Minimum","Mean","Maximum","Upper bound","EBA test","%(+)") %in% names(eba_object)))

  # --- Consistency identity check (slopes only):
  # PM ≈ PIP * PMcon  (in your code: con_PM = PM / PIP)
  # Use tolerance because of rounding (round=12 reduces this problem).
  for (v in c("x1","x2","x3","x4")) {
    pip <- uniform_table[v, "PIP"]
    pm  <- uniform_table[v, "PM"]
    pmc <- uniform_table[v, "PMcon"]

    # skip if pip is extremely small (avoid 0 * Inf style numerical artefacts)
    if (is.finite(pip) && pip > 1e-12) {
      expect_equal(pm, pip * pmc, tolerance = 1e-7)
    }
  }
})

test_that("bma Narrative prior penalizes multi-substitute models in MS=16 space", {
  set.seed(303)

  n <- 25
  X <- replicate(4, rnorm(n))
  colnames(X) <- paste0("x", 1:4)
  y <- 1 + 0.8*X[,1] + rnorm(n, sd = 0.5)
  data <- cbind(y = y, X)

  ms <- model_space(data, M = 4, g = "None", HC = FALSE)

  # Two groups: {x1,x2} are substitutes (group 1), {x3,x4} substitutes (group 2)
  Nar_vec <- c(1L, 1L, 2L, 2L)

  out0 <- bma(ms, EMS = 2, Narrative = 0, dilution = 0, round = 12)
  out1 <- bma(ms, EMS = 2, Narrative = 1, p = c(0.6, 0.5), Nar_vec = Nar_vec, dilution = 0, round = 12)

  pmp0 <- out0[[12]][, 3]  # pmp_uniform
  pmp1 <- out1[[12]][, 3]

  Reg_ID <- ms[[2]][, 1:4, drop = FALSE]

  # Models that include BOTH x1 and x2 (two substitutes in group 1) should be penalized
  idx_two_in_g1 <- which(Reg_ID[,1] == 1L & Reg_ID[,2] == 1L)

  # At least one of these models should have (weakly) lower posterior after narrative penalty
  expect_true(any(pmp1[idx_two_in_g1] <= pmp0[idx_two_in_g1] + 1e-12))

  # Still a valid posterior distribution
  expect_equal(sum(pmp1), 1, tolerance = 1e-10)
})
