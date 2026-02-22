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
#'
#' @return A list with model_space objects: \cr
#' 1. x_names - vector with names of the regressors \cr
#' 2. ols_results - table with the model space - contains ols objects for all the estimated models\cr
#' 3. MS - size of the mode space \cr
#' 4. M - maximum number of regressors in a model \cr
#' 5. K- total number of regressors
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

model_space=function(data, M = NULL, g = "UIP", HC = FALSE){

  # collecting data characteristics
  m <- nrow(data) # number of rows in the data
  n <- ncol(data) # number of columns in the data
  K <- n - 1 # number of regressors

  # What to do if M is not set by the user
  if (is.null(M)){M <- K}

  # collecting names of the used variables
  Var_names <- colnames(data) # names of the variables
  x_names <- Var_names[2:n] # names of the regressors

  if (M > K){# CONDITION about what to do if the user set M that is higher than K (M>K)
    # we tell the user that we are setting M=K
    stop("M>K - maximum number of regressors cannot be bigger than total number of regressors. We set M=K (total number of regressors) and continiue :)")
    M = K # we set M=K
  }# end of the CONDITION about what to do if the user set M that is higher than K (M>K)

  y <- as.matrix(data[,1]) # data on the regressant (dependend variable)
  x <- as.matrix(data[,2:n]) # data on the regressors

  # Total number of models considered
  MS <- sum(choose(K, 0:M))

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


  # NAMES of objects in the ols_results TABLE
  Betas <- matrix(0, nrow = K, ncol = 1)
  SEs <- matrix(0, nrow = K, ncol = 1)
  for (k in 1:K){
    Betas[k,1] = paste0("Coef_", x_names[k])
    SEs[k,1] = paste0("SE_",  x_names[k])
  }
  Betas = rbind("Coef_Const", Betas)
  SEs = rbind("SE_Const", SEs)

  colnames(ols_results) <- rbind(matrix(x_names, nrow = K, ncol = 1),Betas, SEs, matrix(c("log_like", "R^2", "DF", "Dilut"), nrow = 4, ncol =1))

  out<-list(x_names,ols_results,ms,M,K) # we create a model_space object

  return(out)
}
