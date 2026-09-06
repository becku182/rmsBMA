#' OLS calculation with additional objects for model with just a constant.
#'
#' @param y A vector with the dependent variable.
#'
#' @return A list with OLS objects: Coefficients, Standard errors, Marginal likelihood, R^2, Degrees of freedom, Determinant of the regressors matrix.
#' @export
#'
#' @examples
#'
#' x1<-rnorm(10, mean = 0, sd = 1)
#' x2<-rnorm(10, mean = 0, sd = 2)
#' e<-rnorm(10, mean = 0, sd = 0.5)
#' y<-2+x1+2*x2+e
#' x<-cbind(x1,x2)
#' fast_ols_const(y)
#'
fast_ols_const <- function(y){

  # DATA PREPARATION
  y <- as.matrix(y) # changing vector of dependent variables into a matrix
  colnames(y) <- NULL # deleting name of the variables from y
  m <- nrow(y) # number of rows in the data

  Diluntion <- 1  # we assign dilution to be used in dilution prior (George 2010)

  x <- matrix(1, nrow = m, ncol = 1) # adding a constant
  betas = solve(t(x)%*%x)%*%t(x)%*%y # we estimate the model parameters - general version
  y_hat = x%*%betas # we obtain the theoretical values
  res = y - y_hat # we calculate the residuals
  SSR = t(res)%*%res # sum of squares of the residuals
  df = m - 1 # we calculate the number of the degrees of freedom
  sigma2 = (t(res)%*%res)/df # we calculate error variance
  sigma = sigma2^(0.5) # we obtain standard error of the regression
  var_B = as.numeric(sigma2)*solve(t(x)%*%x) # we calculate variances of the coefficients - for general version
  se_B = diag(var_B)^(0.5) # we calculate standard errors of the coefficients
  y_m = mean(y) #we calculate the mean value of the dependent variable
  SST = t(y - y_m*matrix(1, nrow = m, 1))%*%(y - y_m*matrix(1, nrow = m, 1)) # total sum of squares of the regression
  # A constant y fits exactly, so SSR and SST are zero in exact arithmetic.
  # Whether the BLAS returns exactly zero or a rounding-level residual decided
  # whether log_like came out Inf or a large finite number, which made results
  # platform-dependent. Snapping both makes the degenerate case identical
  # everywhere: R2 = NaN, log_like = Inf.
  y_scale <- sum(y^2)
  SSR <- snap_zero(SSR, y_scale, m)
  SST <- snap_zero(SST, y_scale, m)

  R2 <- 1 - (SSR/SST) # calculation of R^2
  log_like <- (-m/2) * log(SSR) # value of the likelihood function (Leamer, 1978)

  # PUTTING ALL THE NECESSARY STUF ONE ONE LIST
  out <- list(betas, se_B,as.numeric(log_like), as.numeric(R2), as.numeric(df), as.numeric(Diluntion)) # creates a list of objects needed for modelSpace function
}
