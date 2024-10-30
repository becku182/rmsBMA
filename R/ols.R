## 2. ols - function that calculates all relevant ols statistics (coefficients, standard errors, Marginal likelihood, R^2, degrees of freedom, determinant of the regressors's matrix,log(Marginal likelihood)) #########

### PARAMETERS OF THE FUNCTION ###
## y - dependent variables ##
## x - matrix of regressors ##
## const should we include a constants

ols<-function(y,x,const){ # function that provides ols estimates and additional statistics
  # HERE SHOULD BE A TEST IF X AND Y ARE OF THE SAME LENGTH

  # DATA PREPARATION
  y<-as.matrix(y) # changing vector of dependent variables into a matrix
  colnames(y)<-NULL # deleting name of the variables from y
  x<-as.matrix(x) # changing table of regressors into a matrix
  colnames(x)<-NULL # deleting name of the variables from x
  m<-nrow(y) # number of rows in the data
  r<-ncol(x) # number of regressors in the model

  if (identical(x,0)){ # when you add "0" as argument model with no variables and a constant will be calculated
    x<-matrix(1,nrow=m,ncol=1) # creation of the vector of ones
    const=0}# we make sure another vector of ones is not gonna be added below
  else if (identical(x,matrix(0,nrow=1,ncol=1))){
    x<-matrix(1,nrow=m,ncol=1) # creation of the vector of ones
    const=0}# we make sure another vector of ones is not gonna be added below

  if (det(t(x)%*%x)==0){stop("Determinant of X'X=0")} # test if X'X can be inverted: Warning that X'X cannot be inverted

  # we need to calculate Dilution before we add vector of ones to the regressors data matrix
  n<-ncol(as.matrix(x))
  # Dilution must be calculated before adding the vector of ones to data matrix x
  Diluntion<-(det(cor(x)))  # we assign dilution to be used in dilution prior (George 2010)

  # ADDING A CONSTANT by adding a vector of ones
  if (const==1){
    ones<-matrix(1,nrow=m,1)
    x<-cbind(ones,x)
  }
  n<-ncol(x) # number of regressors + constant

  betas=solve(t(x)%*%x)%*%t(x)%*%y # we estimate the model parameters - general version
  y_hat=x%*%betas # we obtain the theoretical values
  res=y-y_hat # we calculate the residuals
  SSR=t(res)%*%res # sum of squares of the residuals
  df=m-n # we calculate the number of the degrees of freedom
  sigma2=(t(res)%*%res)/df # we calculate error variance
  sigma=sigma2^(0.5) # we obtain standard error of the regression
  var_B=as.numeric(sigma2)*solve(t(x)%*%x) # we calculate variances of the coefficients - for general version
  se_B=diag(var_B)^(0.5) # we calculate standard errors of the coefficients
  y_m=mean(y) #we calculate the mean value of the dependent variable
  SST=t(y-y_m*matrix(1,nrow=m,1))%*%(y-y_m*matrix(1,nrow=m,1)) # total sum of squares of the regression
  R2<-1-(SSR/SST) # calculation of R^2

  # THERE IS A PROBLEM WITH LIKELIHOOD FUNCTION - NUMBERS TOO CLOSE TO ZERO
  like<-(m^(-r/2))*(SSR^(-m/2)) # value of the likelihood function (Leamer, 1978)
  loglike<-(-r/2)*log(m)+(-m/2)*log(SSR) # ln of "like" above

  # PUTTING ALL THE NECESSARY STUF ONE ONE LIST
  out <- list(betas,se_B,as.numeric(like),as.numeric(R2),as.numeric(df),as.numeric(Diluntion),as.numeric(loglike)) # creates a list of objects needed for modelSpace function
}# THE END OF THE ols FUNCTION
