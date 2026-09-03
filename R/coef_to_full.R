#' Extracting coefficients from g_regression and fast_ols function
#'
#' The function extracts coefficients or standard errors from the results of
#' g_regression and fast_ols function
#'
#' @param model_coefs a vector with estimated coefficients or standard errors
#' @param model_row inclusion vector - row of a model space matrix
#'
#' @return A vector of coefficients coefficients or standard errors from the results of
#' g_regression and fast_ols function
#'
#'
#' @keywords internal

coef_to_full <- function(model_coefs, model_row) {
  model_coefs <- as.numeric(model_coefs)
  model_row   <- as.integer(model_row)

  K <- length(model_row)
  zz <- numeric(K + 1)

  if (length(model_coefs) < 1L) {
    stop("model_coefs must contain at least an intercept in position 1.")
  }

  # intercept
  zz[1] <- model_coefs[1]

  # which regressors are in the model (columns of x)
  idx <- which(model_row == 1L)
  k <- length(idx)

  # slopes from the fitted model
  slopes <- model_coefs[-1L]

  if (length(slopes) != k) {
    stop(sprintf(
      "Mismatch: model_row selects %d regressors, but model_coefs has %d slope coefficients.",
      k, length(slopes)
    ))
  }

  if (k > 0L) {
    zz[idx + 1L] <- slopes
  }

  zz
}
