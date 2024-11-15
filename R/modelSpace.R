#' Calculation of the model space
#'
#'This function calculates all possible models with M regressors that can be constructed out of K regressors. The main object of this function (ols_results) is a table with ols objects for all the estimated models.
#'
#' @param data Data set to work with. The first column is the data for the dependent variable, and the other columns is the data for the regressors.
#' @param M Maximum number of regressor in the estimated models.
#'
#' @return A list with modelSpace objects: \cr
#' 1. x_names - vector with names of the regressors \cr
#' 2. ols_results - table with the model space - contains ols objects for all the estimated models\cr
#' 3. ms - size of the mode space (the number of the last estimated model) \cr
#' 4. M - maximum number of regressors in a model \cr
#' 5. K- total number of regressors
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
#' modelSpace(data,M=3)
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
#' modelSpace(data,M=8)
#'

modelSpace=function(data,M){
  # collecting data characteristics
  m<-nrow(data) # number of rows in the data
  n<-ncol(data) # number of columns in the data
  K<-n-1 # number of regressors

  # collecting names of the used variables
  Var_names<-colnames(data) # names of the variables
  x_names<-Var_names[2:n] # names of the regressors

  if (M>K){# CONDITION about what to do if the user set M that is higher than K (M>K)
    # we tell the user that we are setting M=K
    message("M>K - maximum number of regressors cannot be bigger than total number of regressors. We set M=K and continiue :)")
    M=K # we set M=K
  }# end of the CONDITION about what to do if the user set M that is higher than K (M>K)

  y<-as.matrix(data[,1]) # data on the regressant (dependend variable)
  x<-as.matrix(data[,2:n]) # data on the regressors

  MS=0 # MS a variable representing the size of the model space (total number of models)
  for (k in 0:M){# at this LOOP we add all combinations of regressors up to models with M variables
    c=choose(K,k) # number of models of the size k out of K regressors
    MS=MS+c} # this sum adds up all the models for each possible model size

  # we build a table for all ols statistics
  ols_results=matrix(0,nrow=MS,ncol=3*M+6)

  # ms - model index
  ms=0 #starting value for counting the number of the model

  # k indicates the number of variables in the model
  for (k in 0:M){ # at this LOOP we create all possible model sizes
    c=choose(K,k) # number of models of the size k out of K regressors
    if (k==0){ # CONDITION for the special case of a model with no variables and a constant
      ols1_model<-ols(y=y,x=0,const=1) # estimation of the model with a constant and no regressors
      ms=ms+1 #we change the number of the model
      ols_results[ms,M+1]=as.numeric(ols1_model[1]) # extraction of the coefficients
      ols_results[ms,2*M+2]=as.numeric(ols1_model[2]) # extraction of the standard errors
      ols_results[ms,3*M+3]=as.numeric(ols1_model[3]) # here we extract value of the Likelihood function
      ols_results[ms,3*M+4]=as.numeric(ols1_model[4]) # here we extract R2
      ols_results[ms,3*M+5]=as.numeric(ols1_model[5]) # here we extract the number of degrees of freedom
      ols_results[ms,3*M+6]=as.numeric(ols1_model[6]) # here we extract information for dilution prior
    } # end of the CONDITION for the special case of a model with no variables
    else{# CONDITION for the general case - model with regressors
      models<-as.matrix(utils::combn(1:K,k)) #generates all possible models of the size k out of K regressors
      #i indicates the number of the model from the group of the models with k our of K regressors
      for (i in 1:c){ # at this LOOP we are estimating the individual models with with one or more regressors (k>0)
        mod<-models[,i] # we extract the indices of the regressors to be used in a given model ms
        x_ms=matrix(0,nrow=m,ncol=1) # we crate an artificial vector so we can perform the binding in the next LOOP
        for (t in 1:k){ # at this LOOP we collect the regressors for the model t
          h<-mod[t]  # here we extract an index of the regressor used for estimation
          x_mi=x[,h] # we extract regressor
          x_ms=cbind(x_ms,x_mi) #we bind all k regressors together
        } # end of the LOOP that collects the regressors for the model t
        x_ms<-x_ms[,-1] # we delete an artificial vector from the regressor matrix
        ms=ms+1 # we update the index of the model
        model_ms<-ols(y,x_ms,const=1) #estimation of the model ms
        ols_results[ms,3*M+3]=as.numeric(model_ms[3]) #here we extract value of the Likelihood function
        ols_results[ms,3*M+4]=as.numeric(model_ms[4]) #here we extract R2
        ols_results[ms,3*M+5]=as.numeric(model_ms[5]) #here we extract the number of degrees of freedom
        ols_results[ms,3*M+6]=as.numeric(model_ms[6]) # here we extract information for dilution prior
        for (p in 1:(k+1)){# LOOP performs extraction of the coefficients and the standard errors from the model
          ols_results[ms,M+p]=as.numeric(model_ms[[1]][[p]]) # extraction of the coefficients
          ols_results[ms,2*M+1+p]=as.numeric(model_ms[[2]][[p]]) # extraction of the standard errors
          if (p<k+1){ols_results[ms,p]=mod[p]} # here we extract indices of the used regressors
        }# end of the LOOP performs extraction of the coefficients and the standard errors from the model
      } # end of the LOOP that estimates the individual models with with one or more regressors (k>=0)
    } # end of the CONDITION for the general case - model with regressors
  } # end of the LOOP that creates all model sizes

  # NAMES of objects in the ols_results TABLE
  reg_presence<-matrix(0,nrow=1,ncol=M)
  Betas<-matrix(0,nrow=1,ncol=M)
  SEs<-matrix(0,nrow=1,ncol=M)
  for (m in 1:M){
    reg_presence[1,m]=paste0("Reg_",m)
    Betas[1,m]=paste0("Coef_",m)
    SEs[1,m]=paste0("SE_",m)}
  Betas=cbind("Coef_Const",Betas)
  SEs=cbind("SE_Const",SEs)

  ols_names<-cbind(reg_presence,Betas,SEs,matrix(c("like","R^2","DF","Dilut"),nrow=1,ncol=4))
  colnames(ols_results)<-cbind(reg_presence,Betas,SEs,matrix(c("like","R^2","DF","Dilut"),nrow=1,ncol=4))

  out<-list(x_names,ols_results,ms,M,K) # we create a modelSpace object (mS object)
  return(out)
}
