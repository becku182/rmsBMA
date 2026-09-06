#' OLS calculation with additional objects
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
#' fast_ols(y,x)
#'
fast_ols <- function(y, x){

  y <- as.matrix(y); colnames(y) <- NULL
  x <- as.matrix(x); colnames(x) <- NULL
  m <- nrow(y)
  r <- ncol(x)

  Diluntion <- det(stats::cor(x))  # note: this can be expensive for large r

  # add constant
  x <- cbind(1, x)
  k <- ncol(x)            # = r + 1
  df <- m - k

  # precompute cross-products
  XtX <- crossprod(x)             # X'X
  XtX_inv <- solve(XtX)           # (X'X)^-1
  Xty <- crossprod(x, y)          # X'y

  # OLS coefficients
  betas <- XtX_inv %*% Xty

  # residuals + SSR
  res <- y - x %*% betas
  SSR <- snap_zero(crossprod(res), sum(y^2), m)

  # variance, vcov, se (classical)
  sigma2 <- as.numeric(SSR) / df
  var_B <- sigma2 * XtX_inv
  se_B <- sqrt(diag(var_B))

  # R2
  yc <- y - mean(y)
  SST <- snap_zero(crossprod(yc), sum(y^2), m)
  R2 <- 1 - (SSR / SST)

  # log-likelihood term (as in your code)
  log_like <- (-r/2) * log(m) + (-m/2) * log(SSR)

  list(betas, se_B, as.numeric(log_like), as.numeric(R2), as.numeric(df), as.numeric(Diluntion))
}
