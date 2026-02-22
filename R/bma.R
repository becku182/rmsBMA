#' Calculation of the bma objects
#'
#' This function calculates bma and related objects for the modelSpace object obtained using model_space function.
#'
#' @param modelSpace Model space object (the result of the model_space function)
#' @param EMS Expected model size for model binomial and binomial-beta model prior.
#' @param dilution Binary parameter: 0 - NO application of a dilution prior; 1 - application of a dilution prior (George 2010).
#' @param dil.Par Parameter associated with dilution prior - the exponent of the determinant (George 2010). Used only if parameter dilution=1.
#' @param Narrative Binary parameter: 0 - NO application of a Narrative dilution prior; 1 - application of a Narrative dilution prior.
#' @param p Parameter or vector that indicates by how much we cut probability of a model with substitutes.
#' @param Nar_vec Vector with information on narrative dilution prior where: 0 - the variable has no substitutes; numbers different than 0 denote consecutive groups of variables considered to be substitutes.
#' @param round Parameter indicating to which place the function should round up the results in final tables.
#'
#' @return A list with Posterior objects: \cr
#' 1. pmp_uniform_table - table with results with PMP under binomial model prior \cr
#' 2. pmp_random_table - table with results with PMP under binomial-beta model prior \cr
#' 3. eba_object - table with results of Extreme Bounds Analysis \cr
#' 4. pms_table - table with prior and posterior model sizes \cr
#' 5. x_names - vector with names of the regressors - to be used by the functions \cr
#' 6. K - total number of regressors \cr
#' 7. MS - size of the mode space \cr
#' 8. EMS - expected model size for binomial and binomial-beta model prior specified by the user (default EMS=K/2) \cr
#' 9. dilution - parameter indicating use of dilution \cr
#' 10. for_jointness - table for jointness function \cr
#' 11. for_best_models - table for best_models function \cr
#' 12. for_model_pmp - table for model_pmp function \cr
#' 13. for_model_sizes - table for model_sizes function \cr
#' 14. alphas - vector with the values of the constant \cr
#' 15. betas_nonzero - matrix with coefficients on regressors \cr
#'
#' @export
#'
#' @examples
#' \donttest{
#' data("Trade_data", package = "rmsBMA")
#' data <- Trade_data
#' modelSpace <- model_space(data, M = 6)
#' bma_list <- bma(modelSpace)
#' bma_list[[1]]
#' }


bma <- function(modelSpace,
                EMS = NULL,
                dilution = 0,
                dil.Par = 0.5,
                Narrative = 0,
                p = NULL,
                Nar_vec = NULL,
                round = 6){

  # Extraction of the elements of the mS object
  K <- modelSpace[[5]][1] # extraction of the total number of regressors
  x_names <- modelSpace[[1]][] # extraction of the regressors names
  MS <- modelSpace[[3]][1] # extraction of the total number of models
  M <- modelSpace[[4]][1] # extraction of the maximum number of regressors in the model
  ols_results<-modelSpace[[2]][] # extraction of the ols the model space

  # Dividing ols results into relevant parts
  Reg_ID <- ols_results[,1:K] # we extract vector indices
  betas <- ols_results[,(K+1):(2*K+1)] # we extract coefficients
  VAR <- ols_results[,(2*K+2):(3*K+2)]^2 # we extract standard errors and change to variance
  log_like <- ols_results[,3*K+3] # we extract value of the log likelihood
  R2 <- ols_results[,3*K+4] # we extract R^2
  DF <- ols_results[,3*K+5] # we extract the number of degrees of freedom
  dilut <- ols_results[,3*K+6] # we extract expression associated with dilution prior

  ##### MODEL PRIORS:
  # binomial (Sala-I-Martin et al. 2004) [uniform];
  # binomial-beta (Ley and Steel 2009) [random]

  # Expected model size
  if (is.null(EMS)){EMS = K / 2}
  if (EMS<0.01|EMS>K){
    message("EMS was changed to K/2 (K - number of regressors)")
    EMS = K / 2}

  uniform_models <- matrix(0, nrow = MS, ncol = 1) # vector to store BINOMIAL probabilities ON MODELS
  uniform_sizes <- matrix(0, nrow = M+1, ncol = 1) # vector to store BINOMIAL probabilities ON MODEL SIZES
  random_models <- matrix(0, nrow = MS, ncol = 1) # vector to store BINOMIAL-BETA probabilities ON MODELS
  random_sizes <- matrix(0, nrow = M+1, ncol = 1) # vector to store BINOMIAL-BETA probabilities ON MODEL SIZES

  r <- matrix(apply(Reg_ID, 1, sum), nrow = MS, ncol = 1)

  for (i in 1:MS){
    uniform_models[i,1] = ((EMS/K)^r[i,1])*(1-EMS/K)^(K-r[i,1])
    random_models[i,1] = gamma(1+r[i,1])*gamma((K-EMS)/EMS+K-r[i,1])
  }

  random_models <-random_models/sum(random_models) # here we do scaling
  uniform_models <-uniform_models/sum(uniform_models) # here we do scaling for case M<K

  # TEST CHECKING IF THE USER HAS CHOSEN JUST ONE PRIOR
  if (dilution==1&Narrative==1){stop("Please choose which diltution prior you want to choose:
        regular (dilution=1 and Narrative=0) or narrative (dilution=0 and Narrative=1).
        YOU CANNOT CHOOSE BOTH!")}

  ###### CONDITION for dilution prior
  if (dilution==1){
    uniform_models <- uniform_models*(dilut)^(dil.Par)
    random_models <- random_models*(dilut)^(dil.Par)
    uniform_models <- uniform_models / sum(uniform_models)
    random_models <- random_models / sum(random_models)
  }
  ##################################

  ######### CONDITION for NARRATIVE dilution prior
  if (Narrative == 1) {
    if (is.null(Nar_vec)) stop("Please provide a vector with NARRATIVE INFORMATION through parameter Nar_vec")
    if (length(Nar_vec) != K) stop("Nar_vec is missspecified: Nar_vec should have K elements")

    # NEW: check p length vs number of groups (only if p is not scalar)
    groups <- sort(unique(as.integer(Nar_vec)[as.integer(Nar_vec) > 0L]))
    G <- length(groups)

    if (length(p) > 1L) {
      if (length(p) != G) {
        stop(sprintf(
          "p is misspecified: Nar_vec implies %d groups (positive IDs), but length(p) = %d. Provide scalar p or a vector of length %d.",
          G, length(p), G
        ))
      }
    }

    NarDilution <- group_dilution(Reg_ID, Nar_vec, p)

    uniform_models <- uniform_models * NarDilution
    random_models  <- random_models  * NarDilution
  }
  ##################################

  ##### POSTERIOR MODEL PROBABILITIES
  log_uniform <- log(uniform_models)
  log_random <- log(random_models)
  post_uniform <- (log_like + log_uniform) - max(log_like + log_uniform)*matrix(1, nrow = MS, ncol = 1)
  post_random <- (log_like + log_random) - max(log_like + log_random)*matrix(1, nrow = MS, ncol = 1)
  pmp_uniform <- exp(post_uniform)
  pmp_random <- exp(post_random)
  pmp_uniform <- pmp_uniform/sum(pmp_uniform)
  pmp_random <- pmp_random/sum(pmp_random)
  ##################################

  ##### FOR model_sizes
  sizes <- matrix(choose(K, 0:K), ncol  = 1)
  ind <- matrix(cumsum(sizes), nrow = K+1, ncol = 1) # we create a vector with the number of models in each model size category
  ind <- ind[1:(M+1),1]

  for_sizes <- cbind(rowSums(Reg_ID),uniform_models,random_models)
  for_sizes <- for_sizes[order(for_sizes[,1]), ]
  uniform_models_ordered <- matrix(for_sizes[,2], nrow = MS, ncol = 1)
  random_models_ordered <- matrix(for_sizes[,3], nrow = MS, ncol = 1)
  pmp_uniform_sizes <- matrix(0, nrow = M+1, ncol = 1)
  pmp_random_sizes <- matrix(0, nrow = M+1, ncol = 1)

  for (i in 1:(M+1)){
    if (i==1){uniform_sizes[i,1] = uniform_models_ordered[1,1]
    random_sizes[i,1] = random_models_ordered[1,1]
    pmp_uniform_sizes[i,1] = pmp_uniform[1,1]
    pmp_random_sizes[i,1] = pmp_random[1,1]} # the case of the model with no regressors
    else{uniform_sizes[i,1] = sum(uniform_models_ordered[(ind[i-1]+1):ind[i],1])
    random_sizes[i,1] = sum(random_models_ordered[(ind[i-1]+1):ind[i],1])
    pmp_uniform_sizes[i,1] = sum(pmp_uniform[(ind[i-1]+1):ind[i],1])
    pmp_random_sizes[i,1] = sum(pmp_random[(ind[i-1]+1):ind[i],1])
    } # the case of models with regressors
  }
  ##########################################################

  for_PM_uniform <- matrix(0, nrow = MS, ncol = K+1)
  for_PM_random <- matrix(0, nrow = MS, ncol = K+1)
  for_var_uniform <- matrix(0, nrow = MS, ncol = K+1)
  for_var_random <- matrix(0, nrow = MS, ncol = K+1)

  for (i in 1:K){
    for_PM_uniform[,i] = betas[,i]*pmp_uniform
    for_PM_random[,i] = betas[,i]*pmp_random
    for_var_uniform[,i] = VAR[,i]*pmp_uniform
    for_var_random[,i] = VAR[,i]*pmp_random
  }

  PM_uniform <- matrix(colSums(for_PM_uniform), nrow = K+1, ncol = 1)
  PM_random <- matrix(colSums(for_PM_random), nrow = K+1, ncol = 1)
  var_uniform <- matrix(colSums(for_var_uniform), nrow = K+1, ncol = 1)
  var_random <- matrix(colSums(for_var_random), nrow = K+1, ncol = 1)

  ind_variables <- matrix(betas, nrow = MS, ncol = K+1)
  PM_dev_uniform_prep <- matrix(0, nrow = MS, ncol = K+1)
  PM_dev_random_prep <- matrix(0, nrow = MS, ncol = K+1)

  for (i in 1:(K+1)){
    PM_dev_uniform_prep[,i] <- ((ind_variables[,i] - PM_uniform[i,1])^2)*pmp_uniform
    PM_dev_random_prep[,i] <- ((ind_variables[,i] - PM_random[i,1])^2)*pmp_random
  }

  PM_dev_uniform <- matrix(colSums(PM_dev_uniform_prep), nrow = K+1, ncol = 1)
  PM_dev_random <- matrix(colSums(PM_dev_random_prep), nrow = K+1, ncol = 1)

  VAR_uniform <- PM_dev_uniform + var_uniform
  VAR_random <- PM_dev_random + var_random

  PSD_uniform <- (PM_dev_uniform + var_uniform)^(0.5)
  PSD_random <- (PM_dev_random + var_random)^(0.5)

  beta_MS <- colSums(Reg_ID)[[1]] # Number of model in which a regressor is present
  alphas <- matrix(betas[,1], nrow = MS, ncol = 1)
  just_betas <- matrix(betas[,2:(K+1)], nrow = MS, ncol = K)
  betas_nonzero <- matrix(0, nrow = beta_MS, ncol = K)
  PM_uniform_nonzero <- matrix(0, nrow = beta_MS, ncol = K)
  PM_random_nonzero <- matrix(0, nrow = beta_MS, ncol = K)
  Positive_betas <- matrix(0, nrow = K, ncol = 1)
  Positive_alpha <- mean(alphas > 0)

  for (i in 1:K){
    k=1
    for (j in 1:MS){
      if (just_betas[j,i]!=0){
        betas_nonzero[k,i] = just_betas[j,i]
        PM_uniform_nonzero[k,i] = pmp_uniform[j,1]
        PM_random_nonzero[k,i] = pmp_random[j,1]
        if (just_betas[j,i]>0){
          Positive_betas[i,1] = 1/beta_MS + Positive_betas[i,1]
        }
        k <- k + 1
      }
    }
  }

  PIP_uniform <- rbind(1, matrix(apply(PM_uniform_nonzero, 2, sum), nrow = K, ncol = 1))
  PIP_random <- rbind(1, matrix(apply(PM_random_nonzero, 2, sum), nrow = K, ncol = 1))
  Positive <- rbind(as.numeric(Positive_alpha), Positive_betas)*100

  con_PM_uniform <-  PM_uniform / PIP_uniform
  con_PM_random <-  PM_random / PIP_random

  con_PSD_uniform <- ((VAR_uniform + PM_uniform^2) / PIP_uniform - con_PM_uniform^2)^(0.5)
  con_PSD_random <- ((VAR_random + PM_random^2) / PIP_random - con_PM_random^2)^(0.5)

  # Calculation of P(+) - posterior probability of positive sign
  out <- plus_pmp_from_pmp(pmp_uniform, pmp_random, betas, VAR, DF, Reg_ID)
  Plus_PMP_uniform <- out$Plus_PMP_uniform
  Plus_PMP_random  <- out$Plus_PMP_random

  # Calculation of PSC - posterior sign certanity - probability the sign in a model will be the same as PM
  out <- posterior_sign_certainty(pmp_uniform, pmp_random, betas, VAR, DF, Reg_ID)
  PSC_uniform <- out$PSC_uniform
  PSC_random  <- out$PSC_random

  uniform_table <- cbind(PIP_uniform,PM_uniform,PSD_uniform,con_PM_uniform,con_PSD_uniform,Plus_PMP_uniform,PSC_uniform)
  random_table <- cbind(PIP_random,PM_random,PSD_random,con_PM_random,con_PSD_random,Plus_PMP_random,PSC_random)
  uniform_table <- round(uniform_table, round)
  random_table <- round(random_table, round)

  bma_names <- c("PIP", "PM", "PSD", "PMcon", "PSDcon", "P(+)", "PSC")
  reg_names <- c("CONST", x_names)

  colnames(uniform_table) <- bma_names
  row.names(uniform_table) <- reg_names
  colnames(random_table) <- bma_names
  row.names(random_table) <- reg_names

  uniform_table[1,1] <- NA
  random_table[1,1] <- NA

  # EXTREME BOUND ANALYSIS
  eba_num <- eba(betas, VAR, Reg_ID)   # numeric matrix
  rownames(eba_num) <- c("CONST", x_names)

  # EBA test (character)
  eba_test <- ifelse(sign(eba_num[, 1]) == sign(eba_num[, 5]), "PASS", "FAIL")

  # Round numeric outputs
  eba_num <- round(eba_num, round)
  Positive <- round(Positive, round)

  # Build a nice-looking data frame
  eba_object <- data.frame(
    Variable    = rownames(eba_num),
    `Lower bound` = eba_num[, 1],
    Minimum       = eba_num[, 2],
    Mean          = eba_num[, 3],
    Maximum       = eba_num[, 4],
    `Upper bound` = eba_num[, 5],
    `EBA test`    = eba_test,
    `%(+)`        = Positive,
    row.names = NULL,
    check.names = FALSE
  )

  ### POSTERIOR MODEL TABLE
  PIPs <- cbind(PIP_uniform, PIP_random)
  prior_MS <- matrix(EMS, nrow = 1, ncol = 2)
  posterior_MS <- matrix((colSums(PIPs) - 1), nrow = 1, ncol = 2)
  pms_table <- round(t(rbind(prior_MS, posterior_MS)), round)
  colnames(pms_table) <- c("Prior model size", "Posterior model size")
  row.names(pms_table) <- c("Binomial", "Binomial-beta")

  # for other functions
  for_jointness <- cbind(Reg_ID, pmp_uniform, pmp_random)
  for_best_models <- cbind(Reg_ID, betas, VAR^(0.5), DF, R2, pmp_uniform, pmp_random)
  for_model_pmp <- cbind(uniform_models, random_models, pmp_uniform, pmp_random)
  for_model_sizes <- cbind(uniform_sizes, random_sizes, pmp_uniform_sizes, pmp_random_sizes)

  bma_list <- list(uniform_table,
                   random_table,
                   eba_object,
                   pms_table,
                   x_names,
                   K,
                   MS,
                   EMS,
                   dilution,
                   for_jointness,
                   for_best_models,
                   for_model_pmp,
                   for_model_sizes,
                   alphas,
                   betas_nonzero)

  names(bma_list) <- c("Table with the binomial model prior results",
                       "Table with the binomial model prior results",
                       "Table with Extreme Bounds Analysis results",
                       "Table with prior and posterior model sizes",
                       "Names of variables",
                       "Number of regressors",
                       "Size of the model space (number of models)",
                       "Expected model size",
                       "Parameter indicating use of dilution",
                       "Table for jointness function",
                       "Table for best_models function",
                       "Table for model_pmp function",
                       "Table for model_sizes function",
                       "Vector with the values of the constant",
                       "Matrix with coefficients on regressors")
  return(bma_list)
}
