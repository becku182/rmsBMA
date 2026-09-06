#' OLS calculation with heteroscedasticity consistent covariance matrix (MacKinnon & White 1985).
#'
#' @param y A vector with the dependent variable.
#' @param x A matrix with with regressors as columns.
#'
#' @return A list with OLS objects: Coefficients, Standard errors, Marginal likelihood, R^2, Degrees of freedom, Determinant of the regressors' matrix.
#' @export
#'
#' @examples
#'
#' x1<-rnorm(10, mean = 0, sd = 1)
#' x2<-rnorm(10, mean = 0, sd = 2)
#' e<-rnorm(10, mean = 0, sd = 0.5)
#' y<-2+x1+2*x2+e
#' x<-cbind(x1,x2)
#' fast_ols_HC(y,x)
#'
fast_ols_HC <- function(y, x){

  # DATA PREPARATION
  y <- as.matrix(y)
  colnames(y) <- NULL
  x <- as.matrix(x)
  colnames(x) <- NULL
  m <- nrow(y)
  r <- ncol(x)

  Diluntion <- det(stats::cor(x))  # (George 2010)

  # add constant
  x <- cbind(1, x)

  # compute (X'X)^-1 once
  XtX_inv <- solve(crossprod(x))

  # OLS coefficients
  betas <- XtX_inv %*% crossprod(x, y)

  # fitted + residuals
  y_hat <- x %*% betas
  res <- y - y_hat
  SSR <- snap_zero(crossprod(res), sum(y^2), m)
  df <- m - r - 1

  # HC1 covariance + SE (ONLY)
  u <- as.numeric(res)
  meat <- crossprod(x, x * (u^2))
  var_B <- (m/df) * XtX_inv %*% meat %*% XtX_inv
  se_B <- sqrt(diag(var_B))

  # R^2, log-like
  y_m <- mean(y)
  SST <- snap_zero(crossprod(y - y_m), sum(y^2), m)
  R2 <- 1 - (SSR/SST)
  log_like <- (-r/2) * log(m) + (-m/2) * log(SSR) #(Leamer, 1978)

  out <- list(betas, se_B, as.numeric(log_like), as.numeric(R2), as.numeric(df), as.numeric(Diluntion))
}
