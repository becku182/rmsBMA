#' Calculation of of the jointness measures
#'
#' This function calculates four types of the jointness measures based on the posterior model probabilities (or R^2 measures) calculated using binomial and binomial-beta model prior. The four measures are: \cr
#' 1) HCGHM - for Hofmarcher et al. (2018) measure; \cr
#' 2) LS - for Ley & Steel (2007) measure; \cr
#' 3) DW - for Doppelhofer & Weeks (2009) measure; \cr
#' 4) PPI - for posterior probability of including both variables. \cr
#' The measures will appear in a table above or below the diagonal. \cr
#' \cr
#' REFERENCES \cr
#' Doppelhofer G, Weeks M (2009) Jointness of growth determinants. Journal of Applied Econometrics., 24(2), 209-244. doi: 10.1002/jae.1046 \cr
#' Hofmarcher P, Crespo Cuaresma J, Grün B, Humer S, Moser M (2018) Bivariate jointness measures in Bayesian Model Averaging: Solving the conundrum. Journal of Macroeconomics, 57, 150-165. doi: 10.1016/j.jmacro.2018.05.005 \cr
#' Ley E, Steel M (2007) Jointness in Bayesian variable selection with applications to growth regression. Journal of Macroeconomics, 29(3), 476-493. doi: 10.1016/j.jmacro.2006.12.002
#'
#' @param Post Posterior (Post) object (the result of the Posterior function)
#' @param above Parameter taking one of the four values (1,2,3,4) that indicates what model prior is to be used in the calculation of the measure displayed ABOVE the diagonal.The options are: \cr
#' 1 - for PMP from binomial model prior; \cr
#' 2 - for PMP form binomial-beta model prior; \cr
#' 3 - for R2 from binomial model prior; \cr
#' 4 - for R2 form binomial-beta model prior.
#' @param below Parameter taking one of the four values (1,2,3,4) that indicates what model prior is to be used in the calculation of the measure displayed BELOW the diagonal. The options are: \cr
#' 1 - for PMP from binomial model prior; \cr
#' 2 - for PMP form binomial-beta model prior; \cr
#' 3 - for R2 from binomial model prior; \cr
#' 4 - for R2 form binomial-beta model prior.
#' @param measure Parameter for choosing the measure of jointness: \cr
#' HCGHM - for Hofmarcher et al. (2018) measure; \cr
#' LS - for Ley & Steel (2007) measure; \cr
#' DW - for Doppelhofer & Weeks (2009) measure; \cr
#' PPI - for posterior probability of including both variables.
#' @param rho The parameter "rho" (\eqn{\rho}) to be used in HCGHM jointness measure (default app=0.5). Works only if HCGHM measure is chosen (Hofmarcher et al. 2018).
#' @param app Parameter indicating the decimal place to which the jointness measures should be rounded (default app=3).
#'
#' @return A table with jointness measures for all the pairs of regressors used in the analysis. Parameter "above" indicates what model prior is used for the values ABOVE the diagonal, and parameter "below" indicates what model prior is used for the values BELOW the diagonal.
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
#' Post<-Posterior(modelS)
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

#THE Jointness FUNCTION starts here
Jointness=function(Post,above=1,below=2,measure="HCGHM",rho=0.5,app=3){

  # Extraction of the elements of the mS object
  x_names<-Post[[6]]# we extract vector with names of the regressors from Posterior object (Post object)
  x_names<-x_names[-1]
  M<-Post[[7]]# we extract M - maximum number of regressors in a model from Posterior object (Post object)
  K<-Post[[8]]# we extract K - total number of regressors from Posterior object (Post object)
  MS<-Post[[9]]# we extract MS - size of the mode space from Posterior object (Post object)
  PIPs<-Post[[10]] # we extract matrix with PIPs under different model prior assumptions
  PIPs<-PIPs[-1,] # we delete a row for constants

  # from Posterior object (Post object)
  forJointness<-Post[[12]] # we extract forJointness matrix from Posterior object (Post object)
  Reg_ID<-forJointness[,1:M] # we extract regressor ID
  PMP_uniform<-as.matrix(forJointness[,M+1]) # we extract PMP under unifrom prior
  PMP_random<-as.matrix(forJointness[,M+2]) # we extract PMP under unifrom prior
  PMP_R2_uniform<-as.matrix(forJointness[,M+3]) # we extract PMP under unifrom prior
  PMP_R2_random<-as.matrix(forJointness[,M+4]) # we extract PMP under unifrom prior

  # Information about pairs of regressors: number of pairs, list of all possible regressors pairs
  c=choose(K,2) # number of pairs e.i. jointness measures
  Pairs<-as.matrix(utils::combn(1:K,2)) # a list of all possible pairs of regressors

  # introducing notation a matrices to store posterior objects
  PMP_uniform_a_b<-matrix(0,nrow=c,ncol=1) # a_b - P(a and b) for PMP under uniform prior
  PMP_uniform_Na_b<-matrix(0,nrow=c,ncol=1) # Na_b - P(NOT a and b) for PMP under uniform prior
  PMP_uniform_a_Nb<-matrix(0,nrow=c,ncol=1) # a_Nb - P(a and NOT b) for PMP under uniform prior
  PMP_uniform_Na_Nb<-matrix(0,nrow=c,ncol=1) # Na_Nb - P(NOT a and NOT b) for PMP under uniform prior
  PMP_random_a_b<-matrix(0,nrow=c,ncol=1) # a_b - P(a and b) for PMP under random prior
  PMP_random_Na_b<-matrix(0,nrow=c,ncol=1) # Na_b - P(NOT a and b) for PMP under random prior
  PMP_random_a_Nb<-matrix(0,nrow=c,ncol=1) # a_Nb - P(a and NOT b) for PMP under random prior
  PMP_random_Na_Nb<-matrix(0,nrow=c,ncol=1) # Na_Nb - P(NOT a and NOT b) for PMP under random prior
  R2_uniform_a_b<-matrix(0,nrow=c,ncol=1) # a_b - P(a and b) for R2 under uniform prior
  R2_uniform_Na_b<-matrix(0,nrow=c,ncol=1) # Na_b - P(NOT a and b) for R2 under uniform prior
  R2_uniform_a_Nb<-matrix(0,nrow=c,ncol=1) # a_Nb - P(a and NOT b) for R2 under uniform prior
  R2_uniform_Na_Nb<-matrix(0,nrow=c,ncol=1) # Na_Nb - P(NOT a and NOT b) for R2 under uniform prior
  R2_random_a_b<-matrix(0,nrow=c,ncol=1) # a_b - P(a and b) for R2 under random prior
  R2_random_Na_b<-matrix(0,nrow=c,ncol=1) # Na_b - P(NOT a and b) for R2 under random prior
  R2_random_a_Nb<-matrix(0,nrow=c,ncol=1) # a_Nb - P(a and NOT b) for R2 under random prior
  R2_random_Na_Nb<-matrix(0,nrow=c,ncol=1) # Na_Nb - P(NOT a and NOT b) for R2 under random prior

  for (j in 1:c){# this LOOP goes trough all the regressors pairs
    a<-Pairs[1,j] # ID of the first regressor
    b<-Pairs[2,j] # ID of the second regressor
    for (i in 1:MS){# this LOOP goes through all the models
      aIN<-0 # setting the initial value before the LOOP
      bIN<-0 # setting the initial value before the LOOP
      for (t in 1:M){# this LOOP goes trough model IDs
        if (Reg_ID[i,t]==a){aIN<-1} # CONDITION for presence of the regressors a in the model
        if (Reg_ID[i,t]==b){bIN<-1} # CONDITION for presence of the regressors b in the model
      }# the end of the LOOP that goes trough model IDs
      if (aIN==1&bIN==1){# CONDNION for both variables being included in the model
        PMP_uniform_a_b[j,1]=PMP_uniform_a_b[j,1]+PMP_uniform[i,1]
        PMP_random_a_b[j,1]=PMP_random_a_b[j,1]+PMP_random[i,1]
        R2_uniform_a_b[j,1]=R2_uniform_a_b[j,1]+PMP_R2_uniform[i,1]
        R2_random_a_b[j,1]=R2_random_a_b[j,1]+PMP_R2_random[i,1]
      }# the end of the CONDNION for both variables being included in the model
      else if(aIN==0&bIN==1) {# CONDNION for only regressor b being included in the model
        PMP_uniform_Na_b[j,1]=PMP_uniform_Na_b[j,1]+PMP_uniform[i,1]
        PMP_random_Na_b[j,1]=PMP_random_Na_b[j,1]+PMP_random[i,1]
        R2_uniform_Na_b[j,1]=R2_uniform_Na_b[j,1]+PMP_R2_uniform[i,1]
        R2_random_Na_b[j,1]=R2_random_Na_b[j,1]+PMP_R2_random[i,1]
      }# the end of the CONDNION for only regressor b being included in the model
      else if(aIN==1&bIN==0) {# CONDNION for only regressor a being included in the model
        PMP_uniform_a_Nb[j,1]=PMP_uniform_a_Nb[j,1]+PMP_uniform[i,1]
        PMP_random_a_Nb[j,1]=PMP_random_a_Nb[j,1]+PMP_random[i,1]
        R2_uniform_a_Nb[j,1]=R2_uniform_a_Nb[j,1]+PMP_R2_uniform[i,1]
        R2_random_a_Nb[j,1]=R2_random_a_Nb[j,1]+PMP_R2_random[i,1]
      }# the end of the CONDNION for only regressor a being included in the model
      else if(aIN==0&bIN==0) {# CONDNION for both variables being excluded from the model
        PMP_uniform_Na_Nb[j,1]=PMP_uniform_Na_Nb[j,1]+PMP_uniform[i,1]
        PMP_random_Na_Nb[j,1]=PMP_random_Na_Nb[j,1]+PMP_random[i,1]
        R2_uniform_Na_Nb[j,1]=R2_uniform_Na_Nb[j,1]+PMP_R2_uniform[i,1]
        R2_random_Na_Nb[j,1]=R2_random_Na_Nb[j,1]+PMP_R2_random[i,1]
      }# the end of the CONDNION for both variables being excluded from the model
    }# the end of the LOOP that goes trough all the models
  }# the end of the LOOP that goes trough all the regressors pairs

  # MEASURES OF JOINTNESS
  PMP_uniform_HCGHM_m<-matrix(0,nrow=c,ncol=1) # Hofmarcher et al. (2018)
  PMP_uniform_LS_m<-matrix(0,nrow=c,ncol=1) # Ley & Steel (2007)
  PMP_uniform_DW_m<-matrix(0,nrow=c,ncol=1) # Doppelhofer & Weeks (2009)
  PMP_random_HCGHM_m<-matrix(0,nrow=c,ncol=1) # Hofmarcher et al. (2018)
  PMP_random_LS_m<-matrix(0,nrow=c,ncol=1) # Ley & Steel (2007)
  PMP_random_DW_m<-matrix(0,nrow=c,ncol=1) # Doppelhofer & Weeks (2009)
  R2_uniform_HCGHM_m<-matrix(0,nrow=c,ncol=1) # Hofmarcher et al. (2018)
  R2_uniform_LS_m<-matrix(0,nrow=c,ncol=1) # Ley & Steel (2007)
  R2_uniform_DW_m<-matrix(0,nrow=c,ncol=1) # Doppelhofer & Weeks (2009)
  R2_random_HCGHM_m<-matrix(0,nrow=c,ncol=1) # Hofmarcher et al. (2018)
  R2_random_LS_m<-matrix(0,nrow=c,ncol=1) # Ley & Steel (2007)
  R2_random_DW_m<-matrix(0,nrow=c,ncol=1) # Doppelhofer & Weeks (2009)

  # this LOOP calculates 3 measures of jointness
  for (j in 1:c){# this LOOP goes trough all the regressors pairs
    if (identical(measure,"HCGHM")){ # CONDITION for the HCGHM measure
      PMP_uniform_HCGHM_m[j,1]=((PMP_uniform_a_b[j,1]+rho)*(PMP_uniform_Na_Nb[j,1]+rho)-(PMP_uniform_Na_b[j,1]+rho)*(PMP_uniform_a_Nb[j,1]+rho))/((PMP_uniform_a_b[j,1]+rho)*(PMP_uniform_Na_Nb[j,1]+rho)+(PMP_uniform_Na_b[j,1]+rho)*(PMP_uniform_a_Nb[j,1]+rho)-rho)
      PMP_random_HCGHM_m[j,1]=((PMP_random_a_b[j,1]+rho)*(PMP_random_Na_Nb[j,1]+rho)-(PMP_random_Na_b[j,1]+rho)*(PMP_random_a_Nb[j,1]+rho))/((PMP_random_a_b[j,1]+rho)*(PMP_random_Na_Nb[j,1]+rho)+(PMP_random_Na_b[j,1]+rho)*(PMP_random_a_Nb[j,1]+rho)-rho)
      R2_uniform_HCGHM_m[j,1]=((R2_uniform_a_b[j,1]+rho)*(R2_uniform_Na_Nb[j,1]+rho)-(R2_uniform_Na_b[j,1]+rho)*(R2_uniform_a_Nb[j,1]+rho))/((R2_uniform_a_b[j,1]+rho)*(R2_uniform_Na_Nb[j,1]+rho)+(R2_uniform_Na_b[j,1]+rho)*(R2_uniform_a_Nb[j,1]+rho)-rho)
      R2_random_HCGHM_m[j,1]=((R2_random_a_b[j,1]+rho)*(R2_random_Na_Nb[j,1]+rho)-(R2_random_Na_b[j,1]+rho)*(R2_random_a_Nb[j,1]+rho))/((R2_random_a_b[j,1]+rho)*(R2_random_Na_Nb[j,1]+rho)+(R2_random_Na_b[j,1]+rho)*(R2_random_a_Nb[j,1]+rho)-rho)
    } # the end of the CONDITION for the HCGHM measure
    else if (identical(measure,"LS")){ # CONDITION for the LS measure
      PMP_uniform_LS_m[j,1]=PMP_uniform_a_b[j,1]/(PMP_uniform_Na_b[j,1]+PMP_uniform_a_Nb[j,1])
      PMP_random_LS_m[j,1]=PMP_random_a_b[j,1]/(PMP_random_Na_b[j,1]+PMP_random_a_Nb[j,1])
      R2_uniform_LS_m[j,1]=R2_uniform_a_b[j,1]/(R2_uniform_Na_b[j,1]+R2_uniform_a_Nb[j,1])
      R2_random_LS_m[j,1]=R2_random_a_b[j,1]/(R2_random_Na_b[j,1]+R2_random_a_Nb[j,1])
    } # the end of the CONDITION for the LS measure
    else if (identical(measure,"DW")){ # CONDITION for the DW measure
      PMP_uniform_DW_m[j,1]=log((PMP_uniform_a_b[j,1]/PMP_uniform_Na_b[j,1])*(PMP_uniform_Na_Nb[j,1]/PMP_uniform_a_Nb[j,1]))
      PMP_random_DW_m[j,1]=log((PMP_random_a_b[j,1]/PMP_random_Na_b[j,1])*(PMP_random_Na_Nb[j,1]/PMP_random_a_Nb[j,1]))
      R2_uniform_DW_m[j,1]=log((R2_uniform_a_b[j,1]/R2_uniform_Na_b[j,1])*(R2_uniform_Na_Nb[j,1]/PMP_uniform_a_Nb[j,1]))
      R2_random_DW_m[j,1]=log((R2_random_a_b[j,1]/R2_random_Na_b[j,1])*(R2_random_Na_Nb[j,1]/PMP_random_a_Nb[j,1]))
    } # the end of the CONDITION for the DW measure
  }# the end of the LOOP that goes trough all the regressors pairs

  # here we check which combination of model priors and jointness measures has the user chosen
  if (measure=="HCGHM"&above==1){first<-PMP_uniform_HCGHM_m} #above the diagonal
  if (measure=="HCGHM"&above==2){first<-PMP_random_HCGHM_m} #above the diagonal
  if (measure=="HCGHM"&above==3){first<-R2_uniform_HCGHM_m} #above the diagonal
  if (measure=="HCGHM"&above==4){first<-R2_random_HCGHM_m} #above the diagonal
  if (measure=="LS"&above==1){first<-PMP_uniform_LS_m} #above the diagonal
  if (measure=="LS"&above==2){first<-PMP_random_LS_m} #above the diagonal
  if (measure=="LS"&above==3){first<-R2_uniform_LS_m} #above the diagonal
  if (measure=="LS"&above==4){first<-R2_random_LS_m} #above the diagonal
  if (measure=="DW"&above==1){first<-PMP_uniform_DW_m} #above the diagonal
  if (measure=="DW"&above==2){first<-PMP_random_DW_m} #above the diagonal
  if (measure=="DW"&above==3){first<-R2_uniform_DW_m} #above the diagonal
  if (measure=="DW"&above==4){first<-R2_random_DW_m} #above the diagonal
  if (measure=="PPI"&above==1){first<-PMP_uniform_a_b} #above the diagonal
  if (measure=="PPI"&above==2){first<-PMP_random_a_b} #above the diagonal
  if (measure=="PPI"&above==3){first<-R2_uniform_a_b} #above the diagonal
  if (measure=="PPI"&above==4){first<-R2_random_a_b} #above the diagonal
  if (measure=="HCGHM"&below==1){second<-PMP_uniform_HCGHM_m} #below the diagonal
  if (measure=="HCGHM"&below==2){second<-PMP_random_HCGHM_m} #below the diagonal
  if (measure=="HCGHM"&below==3){second<-R2_uniform_HCGHM_m} #below the diagonal
  if (measure=="HCGHM"&below==4){second<-R2_random_HCGHM_m} #below the diagonal
  if (measure=="LS"&below==1){second<-PMP_uniform_LS_m} #below the diagonal
  if (measure=="LS"&below==2){second<-PMP_random_LS_m} #below the diagonal
  if (measure=="LS"&below==3){second<-R2_uniform_LS_m} #below the diagonal
  if (measure=="LS"&below==4){second<-R2_random_LS_m} #below the diagonal
  if (measure=="DW"&below==1){second<-PMP_uniform_DW_m} #below the diagonal
  if (measure=="DW"&below==2){second<-PMP_random_DW_m} #below the diagonal
  if (measure=="DW"&below==3){second<-R2_uniform_DW_m} #below the diagonal
  if (measure=="DW"&below==4){second<-R2_random_DW_m} #below the diagonal
  if (measure=="PPI"&below==1){second<-PMP_uniform_a_b} #below the diagonal
  if (measure=="PPI"&below==2){second<-PMP_random_a_b} #below the diagonal
  if (measure=="PPI"&below==3){second<-R2_uniform_a_b} #below the diagonal
  if (measure=="PPI"&below==4){second<-R2_random_a_b} #below the diagonal

  # We prepare a table to store the jointness results
  JointnessTable<-matrix(1,nrow=K,ncol=K)

  for (j in 1:c){# this LOOP goes trough all the regressors pairs
    # ABOVE THE DIAGONAL
    JointnessTable[Pairs[1,j],Pairs[2,j]]=first[j,1]
    # BELOW THE DIAGONAL
    JointnessTable[Pairs[2,j],Pairs[1,j]]=second[j,1]
  }# the end of the LOOP that goes trough all the regressors pairs

  # Rounding up the numbers in a table
  JointnessTable<-round(JointnessTable,digits=app)

  # Adding names to rows and columns
  rownames(JointnessTable)<-x_names # we add regressor names to the rows
  colnames(JointnessTable)<-x_names # we add regressor names to the columns

  # creation of the Jointness object (joint object)
  return(JointnessTable)
}
