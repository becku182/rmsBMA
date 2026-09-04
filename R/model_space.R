#' Calculation of the model space
#'
#' This function calculates all possible models with M regressors that can be constructed out of K regressors.
#'
#' @param data Data set to work with. The first column is the data for the dependent variable, and the other columns is the data for the regressors.
#' @param M Maximum number of regressor in the estimated models (default is K - total number of regressors).
#' @param g Value for g in the g prior. Either a number above zero specified by the user or: \cr
#' a) "UIP" for Unit Information Prior (Kass and Wasserman, 1995) \cr
#' b) "RIC" for Risk Inflation Criterion (Foster and George, 1994) \cr
#' c) "Benchmark" for benchmark prior of Fernandez, Ley and Steel (2001) \cr
#' d) "HQ" for prior mimicking Hannan-Quinn information criterion \cr
#' e) "rootUIP" for prior given by the square root of Unit Information Prior \cr
#' f) "None" for the case with no g prior and simple ols regression.
#' In this case the marginal likelihood is calculated according to formula proposed by Leamer (1978).
#' @param HC Logical indicator (default = FALSE) specifying whether a
#' heteroscedasticity-consistent covariance matrix should be used
#' for the estimation of standard errors (MacKinnon & White 1985).
#' @param mc3 Logical (default = FALSE). If TRUE the model space is explored by
#' MC^3 sampling (Madigan and York, 1995) instead of exhaustive enumeration,
#' which makes large K feasible. Reduced model spaces (M < K) are supported:
#' the sampler applies the proposal correction that the size constraint
#' requires. Note that Extreme Bounds Analysis is not available for an MC^3
#' model space.
#' @param draws Number of retained post-burn-in draws (default 10000). Total
#' iterations are draws + burn. Used only when mc3 = TRUE.
#' @param burn Number of initial draws discarded as burn-in (default: equal to
#' draws). Used only when mc3 = TRUE.
#'
#' @return A list with model_space objects: \cr
#' 1. x_names - vector with names of the regressors \cr
#' 2. ols_results - table with the model space - contains ols objects for all the estimated models\cr
#' 3. MS - size of the model space; under mc3 = TRUE this is instead the number of DISTINCT MODELS VISITED \cr
#' 4. M - maximum number of regressors in a model \cr
#' 5. K- total number of regressors \cr
#' 6. mc3 - chain diagnostics (acceptance rate, visit counts, correlation
#' between visit frequencies and analytic posterior mass, mean model size,
#' size of the full space). Present only when mc3 = TRUE.
#'
#' @references
#' Madigan, D. and York, J. (1995). Bayesian graphical models for discrete data.
#' \emph{International Statistical Review}, 63(2), 215-232.
#'
#' @export
#'
#' @examples
#' x1 <- rnorm(20, mean = 0, sd = 1)
#' x2 <- rnorm(20, mean = 0, sd = 2)
#' x3 <- rnorm(20, mean = 0, sd = 3)
#' x4 <- rnorm(20, mean = 0, sd = 1)
#' x5 <- rnorm(20, mean = 0, sd = 2)
#' x6 <- rnorm(20, mean = 0, sd = 4)
#' e <- rnorm(20, mean = 0, sd = 0.5)
#' y <- 2 + x1 + 2*x2 + e
#' data <- cbind(y,x1,x2,x3,x4,x5,x6)
#' modelSpace <- model_space(data, M = 3)
#'
#' # MC^3 sampling of the full model space
#' sampled <- model_space(data, mc3 = TRUE, draws = 2000, burn = 1000)
#' sampled[[6]]$acceptance
#'

model_space=function(data, M = NULL, g = "UIP", HC = FALSE,
                     mc3 = FALSE, draws = NULL, burn = NULL){

  # collecting data characteristics
  m <- nrow(data) # number of rows in the data
  n <- ncol(data) # number of columns in the data
  K <- n - 1 # number of regressors

  M_supplied <- !is.null(M)
  # What to do if M is not set by the user
  if (is.null(M)){M <- K}

  # collecting names of the used variables
  Var_names <- colnames(data) # names of the variables
  x_names <- Var_names[2:n] # names of the regressors

  if (M > K){
    warning("M > K: setting M = K (total number of regressors).")
    M <- K
  }

  y <- as.matrix(data[,1]) # data on the regressant (dependend variable)
  x <- as.matrix(data[,2:n]) # data on the regressors

  # Total number of models considered
  MS <- sum(choose(K, 0:M))

  ## ---- validation, g resolution, and the MC^3 path ------------------------

  if (!is.logical(HC) || length(HC) != 1 || is.na(HC)) {
    stop("Argument 'HC' must be a single logical value (TRUE or FALSE).")
  }
  if (!is.logical(mc3) || length(mc3) != 1 || is.na(mc3)) {
    stop("Argument 'mc3' must be a single logical value (TRUE or FALSE).")
  }

  # 'draws' and 'burn' tune the sampler; they never switch it on by themselves.
  # Silently moving a user from exact enumeration to sampling would change the
  # inference with no signal at all, so a stray 'draws' is an error, not a hint.
  if (!mc3 && (!is.null(draws) || !is.null(burn))) {
    stop("'draws' and 'burn' apply only when mc3 = TRUE. ",
         "Set mc3 = TRUE to sample the model space.")
  }

  g_label <- if (is.character(g)) g else "user-specified"
  g_none <- identical(g, "None")
  if (!g_none) {
    if (is.null(g) || identical(g, "UIP")) {
      g <- 1 / m
    } else if (identical(g, "RIC")) {
      g <- 1 / (K^2)
    } else if (identical(g, "Benchmark")) {
      g <- 1 / max(m, (K^2))
    } else if (identical(g, "HQ")) {
      g <- 1 / (log(m)^3)
    } else if (identical(g, "rootUIP")) {
      g <- sqrt(1 / m)
    } else if (is.numeric(g)) {
      if (length(g) != 1 || !is.finite(g) || g <= 0)
        stop("g must be strictly positive")
    } else {
      g <- 1 / m
    }
  }

  if (mc3) {
    # Reduced model spaces are supported: mc3_sample() applies the
    # |nbd(g)|/|nbd(g')| proposal correction that M < K requires.
    if (M < 1) {
      stop("MC3 requires M >= 1. A model space with M = 0 contains only the ",
           "null model and needs no sampling.")
    }
    if (!M_supplied) {
      message("mc3 = TRUE: sampling the full model space (M = K = ", K, ").")
    }

    if (is.null(draws)) draws <- 10000L
    if (is.null(burn))  burn  <- draws          # 50% burn-in by default
    draws <- suppressWarnings(as.integer(draws))
    burn  <- suppressWarnings(as.integer(burn))
    if (is.na(draws) || draws < 1L)  stop("'draws' must be a positive integer.")
    if (is.na(burn)  || burn  < 0L)  stop("'burn' must be a non-negative integer.")

    fit <- mc3_sample(y, x, K, M, draws = draws, burn = burn,
                      g_none = g_none, g_val = g, HC = HC)

    ols_results <- fit$ols_results
    colnames(ols_results) <- model_space_colnames(x_names, K)

    mc3_info <- list(method     = "mc3",
                     g          = if (g_none) "None" else g,
                     g_label    = g_label,
                     HC         = HC,
                     draws      = draws,
                     burn       = burn,
                     acceptance = fit$acceptance,
                     cor_pmp    = fit$cor_pmp,
                     visits     = fit$visits,
                     mean_size  = fit$mean_size,
                     space_size = MS)

    # Report the chain diagnostics. These decide whether the run is usable at
    # all, so they are shown by default rather than left for the user to dig
    # out of the returned object.
    message(sprintf(
      paste0("MC3: M = %d of K = %d | %d draws after %d burn-in | ",
             "acceptance %.3f | %d distinct models visited | ",
             "mean model size %.2f\nMC3: cor(analytic PMP, visit frequency) = %s"),
      M, K, draws, burn, fit$acceptance, fit$n_models, fit$mean_size,
      if (is.na(fit$cor_pmp)) "NA (too few models)" else sprintf("%.4f", fit$cor_pmp)))

    # The correlation between the analytic posterior mass and the visit
    # frequencies is the convergence check: the two estimate the same thing
    # and agree only once the chain has settled. A low value means the run is
    # too short, not that the model is wrong.
    if (!is.na(fit$cor_pmp) && fit$n_models >= 10L && fit$cor_pmp < 0.99) {
      warning(sprintf(
        paste0("MC3 may not have converged: cor(analytic PMP, visit frequency) ",
               "= %.4f, below 0.99. Increase 'draws'."), fit$cor_pmp),
        call. = FALSE)
    }

    # A chain that almost never moves has explored almost nothing, however
    # well the visited models agree among themselves. This happens when the
    # posterior is very concentrated, and especially at the M boundary, where
    # only deletions are proposed: reaching a different model of size M
    # requires first accepting a worse smaller model.
    if (fit$acceptance < 0.02) {
      warning(sprintf(
        paste0("MC3 accepted only %.3f of proposals and visited %d distinct ",
               "models. The chain barely moved, so the visited set may be a ",
               "poor picture of the model space even though its diagnostics ",
               "look good%s."),
        fit$acceptance, fit$n_models,
        if (M < K) ", and at the M boundary only deletions are proposed"
        else ""),
        call. = FALSE)
    }

    # Element 3 is now the number of DISTINCT MODELS VISITED, not the size of
    # the model space. The full space size is kept in mc3_info$space_size.
    out <- list(x_names = x_names, ols_results = ols_results,
                MS = fit$n_models, M = M, K = K, info = mc3_info)
    return(structure(out, class = "model_space"))
  }
  ## ---- end MC^3 path ------------------------------------------------------


  id_matrix <- model_matrix(K,M)
  for_results <- matrix(0, nrow = MS, ncol = 2*K+6)
  ols_results <- cbind(id_matrix,for_results)
  like_matrix <- matrix(0, nrow = MS, ncol = 1)

  if (!is.logical(HC) || length(HC) != 1 || is.na(HC)) {
    stop("Argument 'HC' must be a single logical value (TRUE or FALSE).")
  }

  if (identical(g, "None")){
    if (HC){
      # Model with no regressors
      ols1_model <- fast_ols_const(y)
      ols_results[1,K+1] = as.numeric(ols1_model[1]) # extraction of the coefficients
      ols_results[1,2*K+2] = as.numeric(ols1_model[2]) # extraction of the standard errors
      ols_results[1,3*K+3] = as.numeric(ols1_model[3]) # here we extract value of the Likelihood function
      ols_results[1,3*K+4] = as.numeric(ols1_model[4]) # here we extract R2
      ols_results[1,3*K+5] = as.numeric(ols1_model[5]) # here we extract the number of degrees of freedom
      ols_results[1,3*K+6] = as.numeric(ols1_model[6]) # here we extract information for dilution prior

      # All other models
      for (ms in 2:MS){
        x_ms <- subset_design(x, ols_results[ms,1:K])
        model_ms <- fast_ols_HC(y, x_ms) #estimation of the model
        ols_results[ms,(K+1):(2*K+1)] <- coef_to_full(model_ms[[1]], ols_results[ms, 1:K]) # extraction of the coefficients
        ols_results[ms,(2*K+2):(3*K+2)] <- coef_to_full(model_ms[[2]], ols_results[ms, 1:K]) # extraction of the standard errors
        ols_results[ms,3*K+3] = as.numeric(model_ms[3]) #here we extract value of the Likelihood function
        ols_results[ms,3*K+4] = as.numeric(model_ms[4]) #here we extract R2
        ols_results[ms,3*K+5] = as.numeric(model_ms[5]) #here we extract the number of degrees of freedom
        ols_results[ms,3*K+6] = as.numeric(model_ms[6]) # here we extract information for dilution prior
      }
    }else{
      # Model with no regressors
      ols1_model <- fast_ols_const(y)
      ols_results[1,K+1] = as.numeric(ols1_model[1]) # extraction of the coefficients
      ols_results[1,2*K+2] = as.numeric(ols1_model[2]) # extraction of the standard errors
      ols_results[1,3*K+3] = as.numeric(ols1_model[3]) # here we extract value of the Likelihood function
      ols_results[1,3*K+4] = as.numeric(ols1_model[4]) # here we extract R2
      ols_results[1,3*K+5] = as.numeric(ols1_model[5]) # here we extract the number of degrees of freedom
      ols_results[1,3*K+6] = as.numeric(ols1_model[6]) # here we extract information for dilution prior

      # All other models
      for (ms in 2:MS){
        x_ms <- subset_design(x, ols_results[ms,1:K])
        model_ms <- fast_ols(y, x_ms) #estimation of the model
        ols_results[ms,(K+1):(2*K+1)] <- coef_to_full(model_ms[[1]], ols_results[ms, 1:K]) # extraction of the coefficients
        ols_results[ms,(2*K+2):(3*K+2)] <- coef_to_full(model_ms[[2]], ols_results[ms, 1:K]) # extraction of the standard errors
        ols_results[ms,3*K+3] = as.numeric(model_ms[3]) #here we extract value of the Likelihood function
        ols_results[ms,3*K+4] = as.numeric(model_ms[4]) #here we extract R2
        ols_results[ms,3*K+5] = as.numeric(model_ms[5]) #here we extract the number of degrees of freedom
        ols_results[ms,3*K+6] = as.numeric(model_ms[6]) # here we extract information for dilution prior
      }
    }
    }else{
    # --- g prior choices ---
    if (is.null(g) || identical(g, "UIP")) {
      g <- 1 / m
    } else if (identical(g, "RIC")) {
      g <- 1 / (K^2)
    } else if (identical(g, "Benchmark")) {
      g <- 1 / max(m, (K^2))
    } else if (identical(g, "HQ")) {
      g <- 1 / (log(m)^3)
    } else if (identical(g, "rootUIP")) {
      g <- sqrt(1 / m)
    } else if (is.numeric(g)) {
      if (length(g) != 1 || !is.finite(g) || g <= 0)
        stop("g must be strictly positive")
    } else {
      g <- 1 / m
    }
    if (HC){
      # Model with no regressors
      ols1_model <- g_regression_fast_const(y,g)
      ols_results[1,K+1] = as.numeric(ols1_model[1]) # extraction of the coefficients
      ols_results[1,2*K+2] = as.numeric(ols1_model[2]) # extraction of the standard errors
      ols_results[1,3*K+3] = as.numeric(ols1_model[3]) # here we extract value of the Likelihood function
      ols_results[1,3*K+4] = as.numeric(ols1_model[4]) # here we extract R2
      ols_results[1,3*K+5] = as.numeric(ols1_model[5]) # here we extract the number of degrees of freedom
      ols_results[1,3*K+6] = as.numeric(ols1_model[6]) # here we extract information for dilution prior

      # All other models
      for (ms in 2:MS){
        x_ms <- subset_design(x, ols_results[ms,1:K])
        model_data <- cbind(y, x_ms)
        model_ms <- g_regression_fast_HC(model_data,g=g) #estimation of the model
        ols_results[ms,(K+1):(2*K+1)] <- coef_to_full(model_ms[[1]], ols_results[ms, 1:K]) # extraction of the coefficients
        ols_results[ms,(2*K+2):(3*K+2)] <- coef_to_full(model_ms[[2]], ols_results[ms, 1:K]) # extraction of the standard errors
        ols_results[ms,3*K+3] = as.numeric(model_ms[3]) #here we extract value of the Likelihood function
        ols_results[ms,3*K+4] = as.numeric(model_ms[4]) #here we extract R2
        ols_results[ms,3*K+5] = as.numeric(model_ms[5]) #here we extract the number of degrees of freedom
        ols_results[ms,3*K+6] = as.numeric(model_ms[6]) # here we extract information for dilution prior
      }
    }else{
      # Model with no regressors
      ols1_model <- g_regression_fast_const(y,g)
      ols_results[1,K+1] = as.numeric(ols1_model[1]) # extraction of the coefficients
      ols_results[1,2*K+2] = as.numeric(ols1_model[2]) # extraction of the standard errors
      ols_results[1,3*K+3] = as.numeric(ols1_model[3]) # here we extract value of the Likelihood function
      ols_results[1,3*K+4] = as.numeric(ols1_model[4]) # here we extract R2
      ols_results[1,3*K+5] = as.numeric(ols1_model[5]) # here we extract the number of degrees of freedom
      ols_results[1,3*K+6] = as.numeric(ols1_model[6]) # here we extract information for dilution prior

      # All other models
      for (ms in 2:MS){
        x_ms <- subset_design(x, ols_results[ms,1:K])
        model_data <- cbind(y, x_ms)
        model_ms <- g_regression_fast(model_data, g=g) #estimation of the model
        ols_results[ms,(K+1):(2*K+1)] <- coef_to_full(model_ms[[1]], ols_results[ms, 1:K]) # extraction of the coefficients
        ols_results[ms,(2*K+2):(3*K+2)] <- coef_to_full(model_ms[[2]], ols_results[ms, 1:K]) # extraction of the standard errors
        ols_results[ms,3*K+3] = as.numeric(model_ms[3]) #here we extract value of the Likelihood function
        ols_results[ms,3*K+4] = as.numeric(model_ms[4]) #here we extract R2
        ols_results[ms,3*K+5] = as.numeric(model_ms[5]) #here we extract the number of degrees of freedom
        ols_results[ms,3*K+6] = as.numeric(model_ms[6]) # here we extract information for dilution prior
      }
    }

  }


  # NAMES of objects in the ols_results TABLE.
  # Shared with the MC^3 path so the two cannot drift apart.
  colnames(ols_results) <- model_space_colnames(x_names, K)

  # Element 6 is always present and records how the space was built, so that
  # summary(), print() and bma() can describe it without guessing. It is the
  # method field, not the length of the list, that tells MC^3 from enumeration.
  info <- list(method = "enumeration",
               g      = if (g_none) "None" else g,
               g_label = g_label,
               HC     = HC)

  out <- list(x_names = x_names, ols_results = ols_results,
              MS = MS, M = M, K = K, info = info)

  return(structure(out, class = "model_space"))
}
