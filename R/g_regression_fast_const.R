#' Fast regression with g prior for a model with just a constant
#'
#' The function implements Bayesian regression with g prior (Zellner, 1986) for a model with just a constant
#'
#' @param y A vector with data - only the dependent variable.
#' @param g Value for g in the g prior. Default value: g = 0.5.
#'
#' @return A list with g_regression objects: \cr
#' 1) Expected values of coefficients \cr
#' 2) Posterior standard errors \cr
#' 3) Natural logarithm of marginal likelihood \cr
#' 4) R^2 form ols model \cr
#' 5) Degrees of freedom \cr
#' 6) Determinant of the regressors' matrix
#'
#' @export
#'
#' @examples
#' x1 <- rnorm(100, mean = 0, sd = 1)
#' x2 <- rnorm(100, mean = 0, sd = 2)
#' e <- rnorm(100, mean = 0, sd = 5)
#' y <- 2 + x1 + 2*x2 + e
#' g_result <- g_regression_fast_const(y, g = 0.99)
#' g_result[[1]]
#' g_result[[2]]
#'
#' x1 <- rnorm(50, mean = 0, sd = 1)
#' x2 <- rnorm(50, mean = 0, sd = 2)
#' e <- rnorm(50, mean = 0, sd = 0.5)
#' y <- 2 + x1 + 2*x2 + e
#' g_result <- g_regression_fast_const(y, g = 1.1)
#' g_result[[1]]
#' g_result[[2]]
#'
g_regression_fast_const <- function(y, g = 0.5) {

  # accept vector or 1-col matrix/data.frame
  y <- as.matrix(y)
  if (ncol(y) != 1L) stop("y must be a vector or a single-column matrix.")
  dimnames(y) <- NULL

  m <- nrow(y)
  if (m < 3L) stop("Need at least 3 observations (m >= 3).")

  Dilution <- 1

  z <- matrix(1, nrow = m, ncol = 1)
  dimnames(z) <- NULL

  df <- m - 1

  ZtZ <- crossprod(z)              # = m
  ZtZ_inv <- 1 / ZtZ               # faster than solve()
  Zty <- crossprod(z, y)           # sum(y)
  yty <- drop(crossprod(y))

  V <- ZtZ_inv / (1 + g)
  betas <- V %*% Zty

  yPzy <- yty - drop(t(Zty) %*% ZtZ_inv %*% Zty)

  y_m <- mean(y)
  y_center_ss <- yty - m * (y_m^2)

  EX_1 <- (1/(g+1)) * yPzy + (g/(g+1)) * y_center_ss

  co_var <- V * ((EX_1) / (m - 2))
  post_se <- sqrt(diag(co_var))

  log_like <- - ((m-1)/2) * log(EX_1)

  R2 <- 1 - (yPzy / y_center_ss)

  list(
    betas,
    post_se,
    as.numeric(log_like),
    as.numeric(R2),
    as.numeric(df),
    as.numeric(Dilution)
  )
}
