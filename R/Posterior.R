#' Calculation of of the posterior objects
#'
#' This function calculates posterior objects for the model space (mS) object obtained using modelSpace function.
#'
#' @param modelSpace Model space (mS) object (the results of the modelSpace function)
#' @param EMS Expected model size for model binomial and binomial-beta model prior. Works only if M=K - Bayesian model averaging on the full model space.
#' @param dilution Binary parameter: 0 - NO application of a dilution prior; 1 - application of a dilution prior (George 2010).
#' @param dil.Par Parameter associated with dilution prior - the exponent of the determinant (George 2010). Used only if parameter dilution=1.
#' @param Narrative Binary parameter: 0 - NO application of a Narrative dilution prior; 1 - application of a Narrative dilution prior.
#' @param p Parameter that indicates by how much we cut probability of a model with substitutes.
#' @param Nar_vec Vector with information on narrative dilution prior where: 0 - the variable has no substitutes; numbers different than 0 denote consecutive groups of variables considered to be substitutes.
#'
#' @return A list with Posterior objects: \cr
#' 1. PMP_uniform_table - table with results with PMP under binomial model prior \cr
#' 2. PMP_random_table - table with results with PMP under binomial-beta model prior \cr
#' 3. EBA - table with results of Extreme Bounds Analysis \cr
#' 4. R2_uniform_table- table with results with R^2 under binomial model model prior \cr
#' 5. R2_random_table - table with results with R^2 under binomial-beta model prior \cr
#' 6. x_names - vector with names of the regressors - to be used by the functions \cr
#' 7. M - maximum number of regressors in a model \cr
#' 8. K - total number of regressors \cr
#' 9. MS - size of the mode space \cr
#' 10. PIPs - table with PIP under different model priors for Jointness function \cr
#' 11. forJointnes - table with model IDs and PMPs for Jointness function \cr
#' 12. forBestModels - table with model IDs, PMPs, coefficients, variances, degrees of freedom, and R^2 for bestModels function \cr
#' 13. sizePriors - table with unifrom and random model priors spread over model sizes for modelSizes function \cr
#' 14. modelPosterior - table with posterior model probabilities for modelSizes function
#'
#' @export
#'
#' @examples
#' x1<-rnorm(20, mean = 0, sd = 1)
#' x2<-rnorm(20, mean = 0, sd = 2)
#' x3<-rnorm(20, mean = 0, sd = 3)
#' x4<-rnorm(20, mean = 0, sd = 1)
#' x5<-rnorm(20, mean = 0, sd = 2)
#' x6<-rnorm(20, mean = 0, sd = 4)
#' e<-rnorm(20, mean = 0, sd = 0.5)
#' y<-2+x1+2*x2+e
#' data<-cbind(y,x1,x2,x3,x4,x5,x6)
#' modelS<-modelSpace(data,M=3)
#' Posterior(modelS)
#'
#' x1<-rnorm(20, mean = 0, sd = 1)
#' x2<-rnorm(20, mean = 0, sd = 2)
#' x3<-rnorm(20, mean = 0, sd = 3)
#' x4<-rnorm(20, mean = 0, sd = 1)
#' x5<-rnorm(20, mean = 0, sd = 2)
#' x6<-rnorm(20, mean = 0, sd = 4)
#' x7<-rnorm(20, mean = 0, sd = 3)
#' x8<-rnorm(20, mean = 0, sd = 1)
#' x9<-rnorm(20, mean = 0, sd = 2)
#' x10<-rnorm(20, mean = 0, sd = 4)
#' e<-rnorm(20, mean = 0, sd = 0.5)
#' y<-2+x3+2*x5+e
#' data<-cbind(y,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10)
#' modelS<-modelSpace(data,M=8)
#' Posterior(modelS)
#'

Posterior=function(modelSpace,EMS=NULL,dilution=0,dil.Par=0.5,Narrative=0,p=0.5,Nar_vec=NULL){

  # Extraction of the elements of the mS object
  K<-modelSpace[[5]][1] # extraction of the total number of regressors
  if (is.null(modelSpace[[1]][])){
    x_names<-matrix(0,nrow=1,ncol=K)
    for (k in 1:K){x_names[1,k]=paste0("k_",k)}
  }else{
    x_names<-modelSpace[[1]][] # extraction of the regressors names from the mS object
  }
  ols_results<-modelSpace[[2]][] # extraction of the ols results (the model space) from the mS object
  MS<-modelSpace[[3]][1] # extraction of the total number of models
  M<-modelSpace[[4]][1] # extraction of the maximum number of regressors in the model

  # Dividing ols results into relevant parts
  Reg_ID<-ols_results[,1:M] # we extract vector indices
  betas<-ols_results[,(M+1):(2*M+1)] # we extract coefficients
  VAR<-ols_results[,(2*M+2):(3*M+2)]^2 # we extract standard errors and change to variance (VAR)
  Like<-ols_results[,3*M+3] # we extract value of the likelihood function (expression proportional to likelihood)
  R2<-ols_results[,3*M+4] # we extract R^2
  DF<-ols_results[,3*M+5] # we extract the number of degrees of freedom
  dilut<-ols_results[,3*M+6] # we extract expression associated with dilution prior

  #### Model Priors
  # We Make separate case for:
  # 1) M<K - here we just assume that all models or model sizes have the same probability
  # 2) M=K - here we introduce binomial and binomial-beta priors with user specified EMS -expected model size
  if (M<K){# CONDITION for Model Priors: how to calculate priors when M<K
    uniform_models<-(1/MS)*matrix(1,nrow=MS,ncol=1) # we create a vector of uniform model priors (ON MODELS)
    uniform_sizes<-matrix(0,nrow=M+1,ncol=1) # we create a vector to store probabilities (ON MODEL SIZES) for model sizes under uniform model prior
    random_models<-matrix(0,nrow=MS,ncol=1) # we create a vector to store probabilities (ON MODELS) for the case of equal probabilities in all model sizes
    random_sizes<-(1/(M+1))*matrix(1,nrow=M+1,ncol=1) # we create a vector of model size priors (ON MODEL SIZES) for the case of equal probabilities in all model sizes

    sizes<-matrix(0,nrow=M+1,ncol=1) # we create vector to store number of models in a given model size

    for (k in 0:M){# at this LOOP we add all combinations of regressors up models with M variables
      sizes[k+1,1]<-choose(K,k) # number of models of the size k out of K regressors
    } # this sum adds up all the models for each possible model size

    ind<-cumsum(sizes) # we create a vector with the number of models in each model size category

    h=1 # here we create artificial variable to perform vector combination through rbind function
    for (i in 1:(M+1)){ # this LOOP creates random prior for individual models and uniform prior for model sizes
      q=matrix(1,nrow=sizes[i,1],ncol=1)*(1/(M+1))*(1/sizes[i,1]) # we create vectors of probabilities for models of a given size
      h=rbind(h,q) # we bind vectors with model random prior probabilities
      if (i==1){uniform_sizes[i,1]=(1/MS)} # we collect probabilities for different model sizes (UNIFORM prior): the case of the model with no regressors
      else{uniform_sizes[i,1]=sum(uniform_models[(ind[i-1]+1):ind[i],1])} # we collect probabilities for different model sizes (UNIFORM prior): the case of models with regressors
    } # end of the LOOP creates random prior for individual models and uniform prior for model sizes
    random_models<-matrix(h[-1,],nrow=MS,ncol=1) # we create a vector to store probabilities (ON MODELS) for the case of equal probabilities in all model sizes

  } else if (M==K) {# CONDITION for Model Priors: how to calculate priors when M=K
    if (is.null(EMS)){EMS=K/2}
    uniform_models<-matrix(0,nrow=MS,ncol=1) # we create a vector to store BINOMIAL probabilities ON MODELS
    uniform_sizes<-matrix(0,nrow=M+1,ncol=1) # we create a vector to store BINOMIAL probabilities ON MODEL SIZES
    random_models<-matrix(0,nrow=MS,ncol=1) # we create a vector to store BINOMIAL-BETA probabilities ON MODELS
    random_sizes<-matrix(0,nrow=M+1,ncol=1) # we create a vector to store BINOMIAL-BETA probabilities ON MODEL SIZES

    km <- matrix(apply(Reg_ID, 1, function(row) sum(row != 0)),nrow=MS,ncol=1)
    for (m in 1:MS){# At this LOOP we calculate expressions proportional to prior model probabilities:
      uniform_models[m,1]=((EMS/K)^km[m,1])*(1-EMS/K)^(K-km[m,1])
      random_models[m,1]=gamma(1+km[m,1])*gamma((K-EMS)/EMS+K-km[m,1])}
    random_models<-random_models/sum(random_models) # here we do scaling

    sizes<-matrix(0,nrow=M+1,ncol=1) # we create vector to store number of models in a given model size

    for (k in 0:M){# at this LOOP we add all combinations of regressors up models with M variables
      sizes[k+1,1]<-choose(K,k) # number of models of the size k out of K regressors
    } # this sum adds up all the models for each possible model size

    ind<-matrix(cumsum(sizes),nrow=M+1,ncol=1) # we create a vector with the number of models in each model size category

    for (i in 1:(M+1)){
      if (i==1){uniform_sizes[i,1]=uniform_models[1,1]
      random_sizes[i,1]=random_models[1,1]} # we collect probabilities for different model sizes (UNIFORM prior): the case of the model with no regressors
      else{uniform_sizes[i,1]=sum(uniform_models[(ind[i-1]+1):ind[i],1])
      random_sizes[i,1]=sum(random_models[(ind[i-1]+1):ind[i],1])
      } # we collect probabilities for different model sizes (UNIFORM prior): the case of models with regressors
    }

  }

  # Posterior model probabilities (PMP) and R^2 weights of individual models
  PMP_uniform<-as.matrix((uniform_models*Like)/sum(uniform_models*Like)) # calculation of PMPs based on uniform prior
  PMP_R2_uniform<-as.matrix((uniform_models*R2)/sum(uniform_models*R2)) # calculation of R^2 weights based on uniform prior
  PMP_random<-as.matrix((random_models*Like)/sum(random_models*Like))  # calculation of PMPs based on random prior
  PMP_R2_random<-as.matrix((random_models*R2)/sum(random_models*R2)) # calculation of R^2 weights based on random prior

  # Creating products of PMPs and betas and stds^2=VAR for calculation of posterior statistics
  betas_PMP_uniform<-matrix(0,nrow=MS,ncol=K+1) # matrix to store products of betas and PMP_uniform
  betas_R2_uniform<-matrix(0,nrow=MS,ncol=K+1) # matrix to store products of betas and PMP_R2_uniform
  betas_PMP_random<-matrix(0,nrow=MS,ncol=K+1) # matrix to store products of betas and PMP_random
  betas_R2_random<-matrix(0,nrow=MS,ncol=K+1) # matrix to store products of betas and PMP_R2_random
  VAR_PMP_uniform<-matrix(0,nrow=MS,ncol=K+1) # matrix to store products of VAR and PMP_uniform
  VAR_R2_uniform<-matrix(0,nrow=MS,ncol=K+1) # matrix to store products of VAR and PMP_R2_uniform
  VAR_PMP_random<-matrix(0,nrow=MS,ncol=K+1) # matrix to store products of VAR and PMP_random
  VAR_R2_random<-matrix(0,nrow=MS,ncol=K+1) # matrix to store products of VAR and PMP_R2_random

  # Multiplication of betas and VARS by PMPs
  for (j in 1:(M+1)){ # at this LOOP we multiply columns of betas and VAR by posterior model measures
    betas_PMP_uniform[,j]=PMP_uniform*betas[,j] # matrix to store products of betas and PMP_uniform
    betas_R2_uniform[,j]=PMP_R2_uniform*betas[,j] # matrix to store products of betas and R2_uniform
    betas_PMP_random[,j]=PMP_random*betas[,j] # matrix to store products of betas and PMP_random
    betas_R2_random[,j]=PMP_R2_random*betas[,j] # matrix to store products of betas and R2_random
    VAR_PMP_uniform[,j]=PMP_uniform*VAR[,j] # matrix to store products of VAR and PMP_uniform
    VAR_R2_uniform[,j]=PMP_R2_uniform*VAR[,j] # matrix to store products of VAR and R2_uniform
    VAR_PMP_random[,j]=PMP_random*VAR[,j] # matrix to store products of VAR and PMP_random
    VAR_R2_random[,j]=PMP_R2_random*VAR[,j] # matrix to store products of VAR and R2_random
  } # end of the LOOP at which we multiply columns of betas and VAR by posterior model measures

  # POSTERIOR INCLUSION PROBABILITIES (PIP) AND POSTERIOR MEANS (PM)
  PM_PMP_uniform<-matrix(0,nrow=K,ncol=1) # matrix to store PMs under uniform model prior
  PM_R2_uniform<-matrix(0,nrow=K,ncol=1) # matrix to store PMs obtained with R^2 under uniform model prior
  PM_PMP_random<-matrix(0,nrow=K,ncol=1) # matrix to store PMs under random model prior
  PM_R2_random<-matrix(0,nrow=K,ncol=1) # matrix to store PMs obtained with R^2 under random model prior
  PIP_PMP_uniform<-matrix(0,nrow=K,ncol=1) # matrix to store PIPs under uniform model prior
  PIP_R2_uniform<-matrix(0,nrow=K,ncol=1) # matrix to store PIPs obtained with R^2 under uniform model prior
  PIP_PMP_random<-matrix(0,nrow=K,ncol=1) # matrix to store PIPs under random model prior
  PIP_R2_random<-matrix(0,nrow=K,ncol=1) # matrix to store PIPs obtained with R^2 under random model prior
  Plus_PMP_uniform<-matrix(0,nrow=K,ncol=1) # matrix to store P(+) under uniform model prior
  Plus_R2_uniform<-matrix(0,nrow=K,ncol=1) # matrix to store P(+) obtained with R^2 under uniform model prior
  Plus_PMP_random<-matrix(0,nrow=K,ncol=1) # matrix to store P(+) under random model prior
  Plus_R2_random<-matrix(0,nrow=K,ncol=1) # matrix to store P(+) obtained with R^2 under random model prior

  Plus_prep<-matrix(0,nrow=MS,ncol=K) # matrix to store ones (1s) representing positive value of the coefficient
  Denominator<-matrix(0,nrow=K,ncol=1) # matrix that counts total number of models with a given variable
  beta_k<-matrix(0,nrow=MS,ncol=K) # matrix to store all the estimated coefficients for regressor k
  se<-matrix(0,nrow=MS,ncol=K) # matrix to store all the estimated standard errors for regressor k

  # Calculation of posterior inclusion probabilities (PIP), posterior means (PM), and preparations for EBA
  for (k in 1:K){ # at this LOOP we move along regressors to collect elements to appropriate sums
    for (i in 1:MS){ # at this LOOP we move along all models to find the models with a given regressor k
      for (t in 1:M){ # at this LOOP we go through Reg_ID to find regressor k index
        if (Reg_ID[i,t]==k){# CONDITION for finding regressor k
          # PREPARATION OF THE OBJECTS FOR EXTREME BOUNDS ANALYSIS (EBA)
          beta_k[i,k]=betas[i,t+1] # we collect all the estimated coefficients for regressor k
          se[i,k]=(VAR[i,t+1])^(0.5) # we collect all the estimated standard errors for regressor k
          Denominator[k,1]=Denominator[k,1]+1 # we get the number of models with a given regressor - the same for every regressors
          # POSTERIOR INCLUSION PROBABILITIES (PIP)
          PIP_PMP_uniform[k,1]=PMP_uniform[i,1]+PIP_PMP_uniform[k,1] # we sum PMPs for models including a given regressor
          PIP_R2_uniform[k,1]=PMP_R2_uniform[i,1]+PIP_R2_uniform[k,1] # we sum PMPs for models including a given regressor
          PIP_PMP_random[k,1]=PMP_random[i,1]+PIP_PMP_random[k,1] # we sum PMPs for models including a given regressor
          PIP_R2_random[k,1]=PMP_R2_random[i,1]+PIP_R2_random[k,1] # we sum PMPs for models including a given regressor
          # POSTERIOR MEANS
          PM_PMP_uniform[k,1]=betas_PMP_uniform[i,t+1]+PM_PMP_uniform[k,1] # we sum products of PMPs with a given coefficient
          PM_R2_uniform[k,1]=betas_R2_uniform[i,t+1]+PM_R2_uniform[k,1] # we sum products of PMPs with a given coefficient
          PM_PMP_random[k,1]=betas_PMP_random[i,t+1]+PM_PMP_random[k,1] # we sum products of PMPs with a given coefficient
          PM_R2_random[k,1]=betas_R2_random[i,t+1]+PM_R2_random[k,1] # we sum products of PMPs with a given coefficient
          # POSTERIOR PROBABILITY OF A POSTERIOR SIGN P(+)
          # here we calculate posterior probability of a positive sign P(+) for the relevant expression
          # see Doppelhofer and Weeks (2009) p. 216 for the details of the expression
          Plus_PMP_uniform[k,1]=Plus_PMP_uniform[k,1]+PMP_uniform[i,1]*stats::pt(betas[i,t+1]/(VAR[i,t+1])^0.5,df =DF[i])
          Plus_R2_uniform[k,1]=Plus_R2_uniform[k,1]+PMP_R2_uniform[i,1]*stats::pt(betas[i,t+1]/(VAR[i,t+1])^0.5,df =DF[i])
          Plus_PMP_random[k,1]=Plus_PMP_random[k,1]+PMP_random[i,1]*stats::pt(betas[i,t+1]/(VAR[i,t+1])^0.5,df =DF[i])
          Plus_R2_random[k,1]=Plus_R2_random[k,1]+PMP_R2_random[i,1]*stats::pt(betas[i,t+1]/(VAR[i,t+1])^0.5,df =DF[i])
          # % OF POSITVE BETAS
          if (betas[i,t+1]>0){# CONDITION for counting models with POSITIVE coefficients
            Plus_prep[i,k]=1 # we set 1 to Plus_prep for models in which regressors have positive coefficients
            # it is latter used to calculate %(+) - percentage of positive coefficients
          } # end of the CONDITION for counting models with POSITIVE coefficients
        } # the end of the CONDITION for finding regressor k
      } # the end of the LOOP we go through Reg_ID to find regressor k index
    } # end of the LOOP that moves along all models to find the models with a given regressor k
  } # end of the LOOP that moves along regressors to collect elements to appropriate sums

  # Distinguishing between the case of positive PM and negative PM - needed for calculation of P(+)
  # for details see Doppelhofer and Weeks (2009) p. 216
  for (k in 1:K){ # at this LOOP we change the expression for P(+) for each regressor with negative PM
    if (PM_PMP_uniform[k,1]<0){Plus_PMP_uniform[k,1]=1-Plus_PMP_uniform[k,1]} # case of PMP and uniform prior
    if (PM_R2_uniform[k,1]<0){Plus_R2_uniform[k,1]=1-Plus_R2_uniform[k,1]} # case of R2 and uniform prior
    if (PM_PMP_random[k,1]<0){Plus_PMP_random[k,1]=1-Plus_PMP_random[k,1]} # case of PMP and random prior
    if (PM_R2_random[k,1]<0){Plus_R2_random[k,1]=1-Plus_R2_random[k,1]} # case of R2 and random prior
  } # the end of the LOOP at which change the expression for P(+) for each regressor with negative PM

  # Plus denotes percentage of positive betas
  Plus<-colSums(Plus_prep)/Denominator

  # POSTERIOR STANDARD DEVIATION (PSD) PREP
  VAR_PMP_prep_uniform<-matrix(0,nrow=K,ncol=1)
  VAR_R2_prep_uniform<-matrix(0,nrow=K,ncol=1)
  VAR_PMP_prep_random<-matrix(0,nrow=K,ncol=1)
  VAR_R2_prep_random<-matrix(0,nrow=K,ncol=1)

  # Calculations of posterior standard deviations (PSD)
  for (k in 1:K){ # at this LOOP we move along regressors to collect elements to appropriate sums
    for (i in 1:MS){ # at this LOOP we move along all models to find the models with a given regressor k
      for (t in 1:M){ # at this LOOP we go through Reg_ID to find regressor k index
        if (Reg_ID[i,t]==k){# CONDITION for finding regressor k
          # Below we add variance of regressor k multiplied by posterior model probability to the square of the difference
          # between coefficient estimated in model i and Posterior mean of regressor k also multiplied by posterior model
          # probability. As a result we get posterior variances we are going to use to calculate posterior standard
          # deviations (PSD) under the assumption of uniform model prior and random model prior (using PMP and R2 in both cases)
          VAR_PMP_prep_uniform[k,1]=VAR_PMP_prep_uniform[k,1]+((betas[i,t+1]-PM_PMP_uniform[k,1])^2)*PMP_uniform[i,1]
          VAR_R2_prep_uniform[k,1]=VAR_R2_prep_uniform[k,1]+((betas[i,t+1]-PM_R2_uniform[k,1])^2)*PMP_R2_uniform[i,1]
          VAR_PMP_prep_random[k,1]=VAR_PMP_prep_random[k,1]+((betas[i,t+1]-PM_PMP_random[k,1])^2)*PMP_random[i,1]
          VAR_R2_prep_random[k,1]=VAR_R2_prep_random[k,1]+((betas[i,t+1]-PM_R2_random[k,1])^2)*PMP_R2_random[i,1]
        } # the end of the CONDITION for finding regressor k
      } # the end of the LOOP we go through Reg_ID to find regressor k index
    } # end of the LOOP that moves along all models to find the models with a given regressor k
  } # end of the LOOP that moves along regressors to collect elements to appropriate sums

  # Calculation of posterior standard deviations (PSD)
  PSD_PMP_uniform<-VAR_PMP_prep_uniform^0.5
  PSD_R2_uniform<-VAR_R2_prep_uniform^0.5
  PSD_PMP_random<-VAR_PMP_prep_random^0.5
  PSD_R2_random<-VAR_R2_prep_random^0.5

  # EXTREME BOUDS ANALYSIS (EBA)
  beta_max<-matrix(0,nrow=K,ncol=1) # matrix to store the highest values of the coefficients
  beta_min<-matrix(0,nrow=K,ncol=1) # matrix to store the lowest values of the coefficients
  beta_mean<-matrix(0,nrow=K,ncol=1) # matrix to store the mean values of the coefficients
  upper<-matrix(0,nrow=K,ncol=1) # matrix to store upper bounds
  lower<-matrix(0,nrow=K,ncol=1) # matrix to store lower bounds

  for (k in 1:K){ # at this LOOP we calculate upper and lower bounds
    beta_max[k,1]=max(beta_k[,k]) # we find the highest values of the coefficients
    beta_min[k,1]=min(beta_k[,k]) # we find the lowest values of the coefficients
    beta_mean[k,1]=mean(beta_k[,k]) # we find mean values of the coefficients
    u<-which(beta_k[,k]==max(beta_k[,k])) # finding the index of the highest value of the coefficient
    if (length(u)>1){u<-1} # CONDITION TO avoid the problem of repeating values - zeros
    upper[k,1]<-beta_k[u,k]+2*se[u,k] # calculation of the upper bound
    l<-which(beta_k[,k]==min(beta_k[,k])) # finding the index of the lowest value of the coefficient
    if (length(l)>1){l<-1}  # CONDITION TO avoid the problem of repeating values - zeros
    lower[k,1]<-beta_k[l,k]-2*se[u,k] # calculation of the lower bound
  } # end of the LOOP where we calculate upper and lower bounds

  # Calculation of the posterior objects (PIP, PM, and PSD) for the constant

  const_betas_PMP_uniform<-PMP_uniform*betas[,1] # we create a vector of products of constants and PMP_uniform
  const_betas_R2_uniform<-PMP_R2_uniform*betas[,1] # we create a vector of products of constants and R2_uniform
  const_betas_PMP_random<-PMP_random*betas[,1] # we create a vector of products of constants and PMP_random
  const_betas_R2_random<-PMP_R2_random*betas[,1] # we create a vector of products of constants and R2_random
  const_VAR_PMP_uniform<-PMP_uniform*VAR[,1] # we create a vector of products of VAR of constants and PMP_uniform
  const_VAR_R2_uniform<-PMP_R2_uniform*VAR[,1] # we create a vector of products of VAR of constants and R2_uniform
  const_VAR_PMP_random<-PMP_random*VAR[,1] # we create a vector of products of VAR of constants and PMP_random
  const_VAR_R2_random<-PMP_R2_random*VAR[,1] # we create a vector of products of VAR of constants and R2_random

  # Calculation of posterior means (PM)
  const_PM_PMP_uniform<-0 # object to store PMs under uniform model prior
  const_PM_R2_uniform<-0 # object to store PMs obtained with R^2 under uniform model prior
  const_PM_PMP_random<-0 # object to store PMs under random model prior
  const_PM_R2_random<-0 # object to store PMs obtained with R^2 under random model
  Plus_const<-0
  Plus_const_PMP_uniform<-0 # object to store P(+) under uniform model prior
  Plus_const_R2_uniform<-0 # object to store P(+) obtained with R^2 under uniform model prior
  Plus_const_PMP_random<-0 # object to store P(+) under random model prior
  Plus_const_R2_random<-0 # object to store P(+) obtained with R^2 under random model

  for (i in 1:MS){ # at this LOOP we move along all models
    # POSTERIOR MEANS
    const_PM_PMP_uniform<-const_betas_PMP_uniform[i,1]+const_PM_PMP_uniform # we sum products of PMPs with a given coefficient
    const_PM_R2_uniform<-const_betas_R2_uniform[i,1]+const_PM_R2_uniform # we sum products of PMPs with a given coefficient
    const_PM_PMP_random<-const_betas_PMP_random[i,1]+const_PM_PMP_random # we sum products of PMPs with a given coefficient
    const_PM_R2_random<-const_betas_R2_random[i,1]+const_PM_R2_random # we sum products of PMPs with a given coefficient
    # POSTERIOR PROBABILITY OF A POSTERIOR SIGN P(+)
    # here we calculate posterior probability of a positive sign P(+) for the relevant expression
    # see Doppelhofer and Weeks (2009) p. 216 for the details of the expression
    Plus_const_PMP_uniform<-Plus_const_PMP_uniform+PMP_uniform[i,1]*stats::pt(betas[i,1]/(VAR[i,1])^0.5,df =DF[i])
    Plus_const_R2_uniform<-Plus_const_R2_uniform+PMP_R2_uniform[i,1]*stats::pt(betas[i,1]/(VAR[i,1])^0.5,df =DF[i])
    Plus_const_PMP_random<-Plus_const_PMP_random+PMP_random[i,1]*stats::pt(betas[i,1]/(VAR[i,1])^0.5,df =DF[i])
    Plus_const_R2_random<-Plus_const_R2_random+PMP_R2_random[i,1]*stats::pt(betas[i,1]/(VAR[i,1])^0.5,df =DF[i])
    # % OF POSITVE BETAS
    if (betas[i,1]>0){ # CONDITION for counting models with POSITIVE coefficients
      Plus_const<-(1/MS)+Plus_const # at this step we add percentage associated with
    } # the end of the CONDITION for counting models with POSITIVE coefficients
  } # end of the LOOP that moves along all models

  # Distinguishing between the case of positive PM and negative PM for the constant - needed for calculation of P(+)
  if (const_PM_PMP_uniform<0){Plus_const_PMP_uniform=1-Plus_const_PMP_uniform} # case of PMP and uniform prior
  if (const_PM_R2_uniform<0){Plus_const_R2_uniform=1-Plus_const_R2_uniform} # case of R2 and uniform prior
  if (const_PM_PMP_random<0){Plus_const_PMP_random=1-Plus_const_PMP_random} # case of PMP and random prior
  if (const_PM_R2_random<0){Plus_const_R2_random=1-Plus_const_R2_random} # case of R2 and random prior

  # Calculations of posterior standard deviations (PSD)
  const_VAR_PMP_prep_uniform<-0 # object to store PSDs under uniform model prior
  const_VAR_R2_prep_uniform<-0 # object to store PSDs obtained with R^2 under uniform model prior
  const_VAR_PMP_prep_random<-0 # object to store PSDs under random model prior
  const_VAR_R2_prep_random<-0 # object to store PSDs obtained with R^2 under random model prior
  for (i in 1:MS){ # at this LOOP we move along all models
    # Below we add variance of the constant multiplied by posterior model probability to the square of the difference
    # between the constant estimated in model i and Posterior mean of the constant k also multiplied by posterior model
    # probability. As a result we get posterior variances we are doing to use to calculate posterior standard
    # deviation (PSD) under the assumption of uniform model prior and random model prior (using PMP and R2 in both cases)
    const_VAR_PMP_prep_uniform<-const_VAR_PMP_prep_uniform+((betas[i,1]-const_PM_PMP_uniform)^2)*PMP_uniform[i,1]
    const_VAR_R2_prep_uniform<-const_VAR_R2_prep_uniform+((betas[i,1]-const_PM_R2_uniform)^2)*PMP_R2_uniform[i,1]
    const_VAR_PMP_prep_random<-const_VAR_PMP_prep_random+((betas[i,1]-const_PM_PMP_random)^2)*PMP_random[i,1]
    const_VAR_R2_prep_random<-const_VAR_R2_prep_random+((betas[i,1]-const_PM_R2_random)^2)*PMP_R2_random[i,1]
  } # end of the LOOP that moves along all models

  # Calculation of posterior standard deviations (PSD)
  const_PSD_PMP_uniform<-const_VAR_PMP_prep_uniform^0.5
  const_PSD_R2_uniform<-const_VAR_R2_prep_uniform^0.5
  const_PSD_PMP_random<-const_VAR_PMP_prep_random^0.5
  const_PSD_R2_random<-const_VAR_R2_prep_random^0.5

  # Extreme Bound Analysis (EBA)
  const_max<-max(betas[,1]) # finding the highest constant
  const_min<-min(betas[,1]) # finding the lowest constant
  const_mean<-mean(betas[,1]) # finding mean constant
  u<-which(betas[,1]==max(betas[,1])) # finding the highest value of the coefficient
  if (length(u)>1){u<-1} # CONDITION TO avoid the problem of repeating values - zeros
  const_upper<-betas[u,1]+2*se[u,1] # calculation of the upper bound
  l<-which(betas[,1]==min(betas[,1])) # finding the lowest value of the coefficient
  if (length(l)>1){l<-1}  # CONDITION TO avoid the problem of repeating values - zeros
  const_lower<-betas[l,1]-2*se[l,1] # calculation of the lower bound

  # At this step we add results for the constant to posterior statistics
  PM_PMP_uniform<-as.matrix(as.numeric(rbind(const_PM_PMP_uniform,PM_PMP_uniform)))
  PM_R2_uniform<-as.matrix(as.numeric(rbind(const_PM_R2_uniform,PM_R2_uniform)))
  PM_PMP_random<-as.matrix(as.numeric(rbind(const_PM_PMP_random,PM_PMP_random)))
  PM_R2_random<-as.matrix(as.numeric(rbind(const_PM_R2_random,PM_R2_random)))
  PIP_PMP_uniform<-as.matrix(as.numeric(rbind(1,PIP_PMP_uniform)))
  PIP_R2_uniform<-as.matrix(as.numeric(rbind(1,PIP_R2_uniform)))
  PIP_PMP_random<-as.matrix(as.numeric(rbind(1,PIP_PMP_random)))
  PIP_R2_random<-as.matrix(as.numeric(rbind(1,PIP_R2_random)))
  PSD_PMP_uniform<-as.matrix(as.numeric(rbind(const_PSD_PMP_uniform,PSD_PMP_uniform)))
  PSD_R2_uniform<-as.matrix(as.numeric(rbind(const_PSD_R2_uniform,PSD_R2_uniform)))
  PSD_PMP_random<-as.matrix(as.numeric(rbind(const_PSD_PMP_random,PSD_PMP_random)))
  PSD_R2_random<-as.matrix(as.numeric(rbind(const_PSD_R2_random,PSD_R2_random)))
  lower<-as.matrix(as.numeric(rbind(const_lower,lower)))
  beta_min<-as.matrix(as.numeric(rbind(const_min,beta_min)))
  beta_mean<-as.matrix(as.numeric(rbind(const_mean,beta_mean)))
  beta_max<-as.matrix(as.numeric(rbind(const_max,beta_max)))
  upper<-as.matrix(as.numeric(rbind(const_upper,upper)))
  Plus<-as.matrix(as.numeric(rbind(Plus_const,Plus)))
  Plus_PMP_uniform<-as.matrix(as.numeric(rbind(Plus_const_PMP_uniform,Plus_PMP_uniform)))
  Plus_R2_uniform<-as.matrix(as.numeric(rbind(Plus_const_R2_uniform,Plus_R2_uniform)))
  Plus_PMP_random<-as.matrix(as.numeric(rbind(Plus_const_PMP_random,Plus_PMP_random)))
  Plus_R2_random<-as.matrix(as.numeric(rbind(Plus_const_R2_random,Plus_R2_random)))
  c_name="Constant"
  x_names<-c(c_name,x_names)

  # Calculation of posterior mean to posterior standard deviation (PM/PSD) ratios
  Ratio_PMP_uniform<-PM_PMP_uniform/PSD_PMP_uniform
  Ratio_R2_uniform<-PM_R2_uniform/PSD_R2_uniform
  Ratio_PMP_random<-PM_PMP_random/PSD_PMP_random
  Ratio_R2_random<-PM_R2_random/PSD_R2_random

  # Objects to be available to the user
  PMP_uniform_table<-cbind(PIP_PMP_uniform,PM_PMP_uniform,PSD_PMP_uniform,Ratio_PMP_uniform,Plus_PMP_uniform)
  rownames(PMP_uniform_table)<-x_names
  colnames(PMP_uniform_table)<-c("PIP","PM","PSD","PM/PSD","P(+)")
  R2_uniform_table<-cbind(PIP_R2_uniform,PM_R2_uniform,PSD_R2_uniform,Ratio_R2_uniform,Plus_R2_uniform)
  rownames(R2_uniform_table)<-x_names
  colnames(R2_uniform_table)<-c("PIP","PM","PSD","PM/PSD","P(+)")
  PMP_random_table<-cbind(PIP_PMP_random,PM_PMP_random,PSD_PMP_random,Ratio_PMP_random,Plus_PMP_random)
  rownames(PMP_random_table)<-x_names
  colnames(PMP_random_table)<-c("PIP","PM","PSD","PM/PSD","P(+)")
  R2_random_table<-cbind(PIP_R2_random,PM_R2_random,PSD_R2_random,Ratio_R2_random,Plus_R2_random)
  rownames(R2_random_table)<-x_names
  colnames(R2_random_table)<-c("PIP","PM","PSD","PM/PSD","P(+)")
  EBA<-cbind(lower,beta_min,beta_mean,beta_max,upper,Plus)
  rownames(EBA)<-x_names
  colnames(EBA)<-c("Lower bound","Minimum","Mean","Maximum","Upper bound","%(+)")

  # objects to be available for other functions
  PIPs<-cbind(PIP_PMP_uniform,PIP_PMP_random,PIP_R2_uniform,PIP_R2_random)# Table with PIP under different model priors for Jointness function: PIPs
  forJointness<-cbind(Reg_ID,PMP_uniform,PMP_random,PMP_R2_uniform,PMP_R2_random)# Table with model IDs and PMPs for Jointness function: forJointnes
  forBestModels<-cbind(Reg_ID,PMP_uniform,PMP_random,PMP_R2_uniform,PMP_R2_random,betas,VAR,R2,DF)# Table with model IDs, PMPs, coefficients, variances, degrees of freedom, and R^2 for bestModels function: forBestModels
  sizePriors<-cbind(uniform_sizes,random_sizes)# Table with unifrom and random model priors spread over model sizes
  modelPosterior<-cbind(PMP_uniform,PMP_random,PMP_R2_uniform,PMP_R2_random)# Table with posterior model probabilities
  # creation of the list of Posterior objects (Post objects)
  out<-list(PMP_uniform_table,PMP_random_table,EBA,R2_uniform_table,R2_random_table,
            x_names,M,K,MS,PIPs,forJointness,forBestModels,sizePriors,modelPosterior) # we create a Posterior object (Post object) - a list with:

  # WE NEED TO ADD Nar_vec AT THE END OF THE LIST

  return(out)

}# THE END of the Posterior FUNCTION
