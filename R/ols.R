#' OLS calculation with additional objects
#'
#' @param y A vector with the dependent variable.
#' @param x A matrix with with regressors as columns or 0 for model without any regressors.
#' @param const Binary variable: 1 - include a constant in the estimation, 0 - do not include a constant in the estimation.
#' @param Norm A parameter used to correct likelihood function when it gets to close to zero in the case of high number of observations.
#'
#' @return A list with OLS objects: Coefficients, Standard errors, Marginal likelihood, R^2, Degrees of freedom, Determinant of the regressors's matrix, log(Marginal likelihood).
#' @export
#'
#' @examples
#' x1<-rnorm(10, mean = 0, sd = 1)
#' x2<-rnorm(10, mean = 0, sd = 2)
#' y<-2+x1+2*x2
#' x<-cbind(x1,x2)
#' const<-1
#' ols(y,x,const)
#'
#' x1<-rnorm(10, mean = 0, sd = 1)
#' x2<-rnorm(10, mean = 0, sd = 2)
#' e<-rnorm(10, mean = 0, sd = 0.5)
#' y<-2+x1+2*x2+e
#' x<-cbind(x1,x2)
#' const<-1
#' ols(y,x,const)
#'
ols <- function(y, x, const, Norm = NULL){

  # DATA PREPARATION
  y <- as.matrix(y) # changing vector of dependent variables into a matrix
  colnames(y) <- NULL # deleting name of the variables from y
  x <- as.matrix(x) # changing table of regressors into a matrix
  colnames(x) <- NULL # deleting name of the variables from x
  m <- nrow(y) # number of rows in the data
  r <- ncol(x) # number of regressors in the model

  if (identical(x, 0)){ # when you add "0" as argument model with no variables and a constant will be calculated
    x <- matrix(1, nrow = m, ncol = 1) # creation of the vector of ones
    const = 0}# we make sure another vector of ones is not gonna be added below
  else if (identical(x, matrix(0, nrow = 1, ncol = 1))){
    x <- matrix(1, nrow = m, ncol = 1) # creation of the vector of ones
    const = 0}# we make sure another vector of ones is not gonna be added below

  # we need to calculate Dilution before we add vector of ones to the regressors data matrix
  n <- ncol(as.matrix(x))
  # Dilution must be calculated before adding the vector of ones to data matrix x
  Diluntion <- (det(stats::cor(x)))  # we assign dilution to be used in dilution prior (George 2010)

  # ADDING A CONSTANT by adding a vector of ones
  if (const==1){
    x <- cbind(matrix(1, nrow = m, ncol = 1), x)
  }
  n <- ncol(x) # number of regressors + constant

  betas = solve(t(x)%*%x)%*%t(x)%*%y # we estimate the model parameters - general version
  y_hat = x%*%betas # we obtain the theoretical values
  res = y - y_hat # we calculate the residuals
  SSR = t(res)%*%res # sum of squares of the residuals
  df = m - n # we calculate the number of the degrees of freedom
  sigma2 = (t(res)%*%res)/df # we calculate error variance
  sigma = sigma2^(0.5) # we obtain standard error of the regression
  var_B = as.numeric(sigma2)*solve(t(x)%*%x) # we calculate variances of the coefficients - for general version
  se_B = diag(var_B)^(0.5) # we calculate standard errors of the coefficients
  y_m = mean(y) #we calculate the mean value of the dependent variable
  SST = t(y - y_m*matrix(1, nrow = m, 1))%*%(y - y_m*matrix(1, nrow = m, 1)) # total sum of squares of the regression
  R2 <- 1 - (SSR/SST) # calculation of R^2

  # THERE IS A PROBLEM WITH LIKELIHOOD FUNCTION - NUMBERS TOO CLOSE TO ZERO

  # When the number of observations is relatively low
  if (is.null(Norm)==1){
    like <- (m^(-r/2))*(SSR^(-m/2)) # value of the likelihood function (Leamer, 1978)
  }
  # When the number of observations is relatively high
  if (is.null(Norm)==0){
    like <- (m^(-r/2)*(Norm*SSR)^(-m/2)) # value of the likelihood function with a correcting constant
  }

  # PUTTING ALL THE NECESSARY STUF ONE ONE LIST
  out <- list(betas, se_B,as.numeric(like), as.numeric(R2), as.numeric(df), as.numeric(Diluntion)) # creates a list of objects needed for modelSpace function
}
