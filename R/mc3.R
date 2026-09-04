# Internal machinery for MC^3 (Markov Chain Monte Carlo Model Composition)
# sampling over the model space. See model_space(mc3 = TRUE).
#
# Nothing in this file is exported.

#' Key identifying a model by its inclusion vector
#' @param incl integer inclusion vector of length K
#' @return a single character string
#' @keywords internal
#' @noRd
mc3_key <- function(incl) paste0(as.integer(incl), collapse = "")

#' Estimate one model and return it as a model_space row
#'
#' Mirrors exactly what model_space() writes for a single model, so that a row
#' produced here is indistinguishable from the corresponding enumerated row.
#' The layout of the returned vector (length 3K + 6) is:
#'
#' ```
#'   1:K              inclusion indicators
#'   (K+1):(2K+1)     coefficients (constant first)
#'   (2K+2):(3K+2)    standard errors (constant first)
#'   3K+3             log marginal likelihood
#'   3K+4             R^2
#'   3K+5             degrees of freedom
#'   3K+6             dilution-prior term
#' ```
#'
#' @param y response, as a one-column matrix
#' @param x matrix of all K candidate regressors
#' @param incl inclusion vector of length K
#' @param K total number of regressors
#' @param g_none TRUE when g = "None" (plain OLS / Leamer marginal likelihood)
#' @param g_val resolved numeric value of g, ignored when g_none is TRUE
#' @param HC TRUE for heteroscedasticity-consistent standard errors
#' @return numeric vector of length 3K + 6
#' @keywords internal
#' @noRd
mc3_fit_row <- function(y, x, incl, K, g_none, g_val, HC) {

  row <- numeric(3 * K + 6)
  row[1:K] <- as.integer(incl)

  if (sum(incl) == 0L) {
    # Model with no regressors. Note that model_space() uses the g-prior
    # constant-only fit whenever g != "None", regardless of HC.
    fit <- if (g_none) fast_ols_const(y) else g_regression_fast_const(y, g_val)
    row[K + 1]     <- as.numeric(fit[1])
    row[2 * K + 2] <- as.numeric(fit[2])
  } else {
    x_ms <- subset_design(x, incl)
    fit <-
      if (g_none) {
        if (HC) fast_ols_HC(y, x_ms) else fast_ols(y, x_ms)
      } else {
        model_data <- cbind(y, x_ms)
        if (HC) g_regression_fast_HC(model_data, g = g_val)
        else    g_regression_fast(model_data, g = g_val)
      }
    row[(K + 1):(2 * K + 1)]     <- coef_to_full(fit[[1]], incl)
    row[(2 * K + 2):(3 * K + 2)] <- coef_to_full(fit[[2]], incl)
  }

  row[3 * K + 3] <- as.numeric(fit[3])   # log marginal likelihood
  row[3 * K + 4] <- as.numeric(fit[4])   # R^2
  row[3 * K + 5] <- as.numeric(fit[5])   # degrees of freedom
  row[3 * K + 6] <- as.numeric(fit[6])   # dilution term
  row
}

#' MC^3 sampler over the full model space (M = K)
#'
#' Metropolis-Hastings random walk over models, in the manner of
#' Madigan and York (1995). From the current model a neighbour is proposed by
#' flipping one of the K inclusion indicators, chosen uniformly at random.
#'
#' The proposal is symmetric here BECAUSE M = K: every model has exactly K
#' single-flip neighbours, so q(g' | g) = q(g | g') = 1/K and the Hastings
#' ratio is 1. This is precisely the property that fails once M < K, where
#' models of size M have only M neighbours (drops only) while smaller models
#' have K, and the acceptance ratio then requires the correction
#' |nbd(g)| / |nbd(g')|. Constrained MC^3 must not reuse this function without
#' adding that term.
#'
#' The chain targets the posterior under a UNIFORM prior over models, i.e.
#' proportional to the marginal likelihood alone. Model priors are applied
#' afterwards in bma(), which renormalises over the visited models. Note that
#' the binomial model prior with EMS = K/2 is exactly uniform over models, so
#' for the default EMS the reweighting is exact.
#'
#' @param y response, as a one-column matrix
#' @param x matrix of all K candidate regressors
#' @param K total number of regressors
#' @param draws number of retained post-burn-in draws
#' @param burn number of burn-in draws, discarded
#' @param g_none TRUE when g = "None"
#' @param g_val resolved numeric value of g
#' @param HC TRUE for heteroscedasticity-consistent standard errors
#' @return list with the visited-model matrix and chain diagnostics
#' @keywords internal
#' @noRd
mc3_sample <- function(y, x, K, draws, burn, g_none, g_val, HC) {

  total <- burn + draws
  cache  <- new.env(hash = TRUE, parent = emptyenv())  # key -> fitted row
  visits <- new.env(hash = TRUE, parent = emptyenv())  # key -> post-burn-in count

  fitted_row <- function(key, incl) {
    if (exists(key, envir = cache, inherits = FALSE)) return(get(key, envir = cache))
    r <- mc3_fit_row(y, x, incl, K, g_none, g_val, HC)
    assign(key, r, envir = cache)
    r
  }

  curr      <- integer(K)                       # start from the null model
  curr_key  <- mc3_key(curr)
  curr_ll   <- fitted_row(curr_key, curr)[3 * K + 3]
  if (!is.finite(curr_ll))
    stop("MC3 could not start: the null model has a non-finite marginal likelihood.")

  n_accept  <- 0L
  size_path <- integer(total)

  for (it in seq_len(total)) {

    j       <- sample.int(K, 1L)
    prop    <- curr
    prop[j] <- 1L - prop[j]
    prop_key <- mc3_key(prop)
    prop_ll  <- fitted_row(prop_key, prop)[3 * K + 3]

    # Symmetric proposal and uniform model prior, so the acceptance ratio is
    # just the marginal likelihood ratio. Non-finite proposals are rejected
    # (e.g. a singular design).
    if (is.finite(prop_ll) && log(stats::runif(1)) < (prop_ll - curr_ll)) {
      curr     <- prop
      curr_key <- prop_key
      curr_ll  <- prop_ll
      n_accept <- n_accept + 1L
    }

    size_path[it] <- sum(curr)

    if (it > burn) {
      cnt <- if (exists(curr_key, envir = visits, inherits = FALSE))
        get(curr_key, envir = visits) else 0L
      assign(curr_key, cnt + 1L, envir = visits)
    }
  }

  keys   <- ls(visits)
  counts <- vapply(keys, function(k) get(k, envir = visits), numeric(1))
  rows   <- lapply(keys, function(k) get(k, envir = cache))
  ols_results <- do.call(rbind, rows)

  # Deterministic ordering: by model size, then by key. Puts the null model
  # first when it was visited, mirroring the enumerated layout.
  ord <- order(rowSums(ols_results[, 1:K, drop = FALSE]), keys)
  ols_results <- ols_results[ord, , drop = FALSE]
  counts <- counts[ord]

  # Convergence diagnostic: visit frequencies against analytic posterior mass
  # over the same visited set. These should agree closely once the chain has
  # converged; a low correlation means it has not.
  ll  <- ols_results[, 3 * K + 3]
  aw  <- exp(ll - max(ll)); aw <- aw / sum(aw)
  vf  <- counts / sum(counts)
  cor_pmp <- if (length(ll) > 2L) suppressWarnings(stats::cor(aw, vf)) else NA_real_

  list(ols_results   = ols_results,
       visits        = as.numeric(counts),
       n_models      = nrow(ols_results),
       acceptance    = n_accept / total,
       cor_pmp       = cor_pmp,
       draws         = draws,
       burn          = burn,
       mean_size     = mean(size_path[(burn + 1L):total]))
}

#' Column names for a model_space results matrix
#'
#' Kept here so that the enumerated and MC^3 paths cannot drift apart.
#'
#' @param x_names character vector of regressor names
#' @param K total number of regressors
#' @return character vector of length 3K + 6
#' @keywords internal
#' @noRd
model_space_colnames <- function(x_names, K) {
  c(x_names,
    "Coef_Const", paste0("Coef_", x_names),
    "SE_Const",   paste0("SE_",   x_names),
    "log_like", "R^2", "DF", "Dilut")
}
