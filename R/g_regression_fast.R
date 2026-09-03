#' Fast regression with g prior
#'
#' The function implements Bayesian regression with g prior (Zellner, 1986)
#'
#' @param data A matrix with data. The first column is interpreted as with the dependent variable, while
#' the remaining columns are interpreted as regressors.
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
#' data <- cbind(y,x1,x2)
#' g_result <- g_regression_fast(data, g = 0.99)
#' g_result[[1]]
#' g_result[[2]]
#'
#' x1 <- rnorm(50, mean = 0, sd = 1)
#' x2 <- rnorm(50, mean = 0, sd = 2)
#' e <- rnorm(50, mean = 0, sd = 0.5)
#' y <- 2 + x1 + 2*x2 + e
#' data <- cbind(y,x1,x2)
#' g_result <- g_regression_fast(data, g = 1.1)
#' g_result[[1]]
#' g_result[[2]]
#'
g_regression_fast <- function(data, g = 0.5) {

  # --- DATA PREPARATION (fast + no names) ---
  data <- as.matrix(data)
  dimnames(data) <- NULL

  m <- nrow(data)
  col <- ncol(data)

  y <- matrix(data[, 1], ncol = 1)
  x <- data[, 2:col, drop = FALSE]

  dimnames(y) <- NULL
  dimnames(x) <- NULL

  r <- col - 1

  # --- Dilution ---
  Dilution <- det(stats::cor(x))

  # --- Add intercept ---
  z <- cbind(1, x)
  dimnames(z) <- NULL

  k <- ncol(z)
  df <- m - k

  # --- Crossproducts ---
  ZtZ <- crossprod(z)
  if (rcond(ZtZ) < .Machine$double.eps) stop("Determinant of Z'Z=0")

  ZtZ_inv <- solve(ZtZ)
  Zty <- crossprod(z, y)
  yty <- drop(crossprod(y))

  # --- Posterior quantities ---
  V <- ZtZ_inv / (1 + g)
  betas <- V %*% Zty

  # y'Pz y
  yPzy <- yty - drop(t(Zty) %*% ZtZ_inv %*% Zty)

  y_m <- mean(y)
  y_center_ss <- yty - m * (y_m^2)

  EX_1 <- (1/(g+1)) * yPzy + (g/(g+1)) * y_center_ss

  se_2 <- EX_1 / m
  co_var <- V * ((m * se_2) / (m - 2))
  post_se <- sqrt(diag(co_var))

  log_like <- (r/2) * log(g/(g+1)) - ((m-1)/2) * log(EX_1)

  SSR <- yPzy
  SST <- y_center_ss
  R2 <- 1 - (SSR / SST)

  list(
    betas,
    post_se,
    as.numeric(log_like),
    as.numeric(R2),
    as.numeric(df),
    as.numeric(Dilution)
  )
}
