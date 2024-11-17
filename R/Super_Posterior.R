#' Calculation of of the posterior objects
#'
#' This function calculates posterior objects for the data set. \cr
#' If the data is in the panel form the function assumes it has the following structure\cr
#' \cr
#' section_1  year_1   y x1 x2 x3 ....\cr
#' section_2  year_1   y x1 x2 x3 ....\cr
#' section_3  year_1   y x1 x2 x3 ....\cr
#' ........\cr
#' section_n  year_1   y x1 x2 x3 ....\cr
#' section_1  year_2   y x1 x2 x3 ....\cr
#' section_2  year_2   y x1 x2 x3 ....\cr
#' section_3  year_2   y x1 x2 x3 ....\cr
#' ........\cr
#' section_n  year_2  y x1 x2 x3 ....\cr
#' ........\cr
#' section_n  year_(T-1)   y x1 x2 x3 ....\cr
#' section_1  year_T       y x1 x2 x3 ....\cr
#' section_2  year_T       y x1 x2 x3 ....\cr
#' section_3  year_T       y x1 x2 x3 ....\cr
#' ........\cr
#'section_n  year_T       y x1 x2 x3 ....\cr
#'
#' @param data Data set to work with. The first column is the data for the dependent variable, and the other columns is the data for the regressors.
#' @param M Maximum number of regressor in the estimated models.
#' @param EMS Expected model size for model binomial and binomial-beta model prior. Works only if M=K - Bayesian model averaging on the full model space.
#' @param dilution Binary parameter: 0 - NO application of a dilution prior; 1 - application of a dilution prior (George 2010).
#' @param dil.Par Parameter associated with dilution prior - the exponent of the determinant (George 2010). Used only if parameter dilution=1.
#' @param Narrative Binary parameter: 0 - NO application of a Narrative dilution prior; 1 - application of a Narrative dilution prior.
#' @param p Parameter that indicates by how much we cut probability of a model with substitutes.
#' @param Nar_vec Vector with information on narrative dilution prior where: 0 - the variable has no substitutes; numbers different than 0 denote consecutive groups of variables considered to be substitutes.
#' @param FE Binary variable: 1 - include fixed effect, 0 - do not include fixed effects.
#' @param Time The number of time periods - works only if FE=1.
#' @param Section The number of cross-sections - works only if EF=1.
#' @param Time_FE Binary variable: 1 - include time fixed effect, 0 - do not include time fixed effects. Works only if EF=1.
#' @param Section_FE Binary variable: 1 - include cross-section fixed effect, 0 - do not include cross-section fixed effects. Works only if EF=1.
#' @param STD Binary variable: 1 - standardize the data set, 0 - do not standardize the data set. By standardization we mean subracion of amean and division  by standard deviation of each variable.
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
#' 13. sizePriors - table with uniform and random model priors spread over model sizes for modelSizes function \cr
#' 14. modelPosterior - table with posterior model probabilities for modelSizes function \cr
#' 15. EMS - expected model size for binomial and binomial=beta model prior specified by the user (default EMS=K/2). Works only if M=K \cr
#' 16. beta_k - table with coefficients from all the estimated models \cr
#' 17. NarDilution - vector with factors that multiply model priors obtained from Narrative approach (appears only when parameter Narrative=1)
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
#' M<-3
#' Super_Posterior(data,M)
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
#' M<-8
#' Super_Posterior(data,M)
#'
#' x1<-rnorm(20, mean = 0, sd = 1)
#' x2<-rnorm(20, mean = 0, sd = 2)
#' x3<-rnorm(20, mean = 0, sd = 3)
#' x4<-rnorm(20, mean = 0, sd = 1)
#' x5<-rnorm(20, mean = 0, sd = 2)
#' x6<-rnorm(20, mean = 0, sd = 4)
#' e<-rnorm(20, mean = 0, sd = 0.5)
#' y<-2+x1+2*x2+e
#' data<-cbind(y,x1,x2,x3,x4,x5,x6)
#' colnames(data)<-c("y","x1","x2","x3","x4","x5","x6")
#' M<-6
#' Super_Post<-Super_Posterior(data,M,dilution=1,dil.Par=0.5)
#'
#' x1<-rnorm(50, mean = 0, sd = 5)
#' x2<-rnorm(50, mean = 0, sd = 2)
#' x3<-rnorm(50, mean = 0, sd = 7)
#' x4<-rnorm(50, mean = 0, sd = 1)
#' x5<-rnorm(50, mean = 0, sd = 3)
#' x6<-rnorm(50, mean = 0, sd = 4)
#' e<-rnorm(50, mean = 0, sd = 20)
#' y<-2+x1+2*x2+e
#' data<-cbind(y,x1,x2,x3,x4,x5,x6)
#' colnames(data)<-c("y","x1","x2","x3","x4","x5","x6")
#' M<-5
#' Super_Post<-Super_Posterior(data,M,dilution=1,dil.Par=0.5)
#'
#' x1<-rnorm(20, mean = 0, sd = 1)
#' x2<-rnorm(20, mean = 0, sd = 2)
#' x3<-rnorm(20, mean = 0, sd = 3)
#' x4<-rnorm(20, mean = 0, sd = 1)
#' x5<-rnorm(20, mean = 0, sd = 2)
#' x6<-rnorm(20, mean = 0, sd = 4)
#' e<-rnorm(20, mean = 0, sd = 0.5)
#' y<-2+x1+2*x2+e
#' data<-cbind(y,x1,x2,x3,x4,x5,x6)
#' colnames(data)<-c("y","x1","x2","x3","x4","x5","x6")
#' M<-6
#' Nar_vec<-as.matrix(c(0,1,1,1,2,2))
#' Super_Post<-Super_Posterior(data,M,Narrative=0,p=0.5,Nar_vec=Nar_vec)
#'
#' x1<-rnorm(50, mean = 0, sd = 5)
#' x2<-rnorm(50, mean = 0, sd = 2)
#' x3<-rnorm(50, mean = 0, sd = 7)
#' x4<-rnorm(50, mean = 0, sd = 1)
#' x5<-rnorm(50, mean = 0, sd = 3)
#' x6<-rnorm(50, mean = 0, sd = 4)
#' e<-rnorm(50, mean = 0, sd = 20)
#' y<-2+x1+2*x2+e
#' data<-cbind(y,x1,x2,x3,x4,x5,x6)
#' colnames(data)<-c("y","x1","x2","x3","x4","x5","x6")
#' M<-5
#' Nar_vec<-as.matrix(c(0,1,1,1,2,2))
#' Super_Post<-Super_Posterior(data,M,Narrative=0,p=0.5,Nar_vec=Nar_vec)
#'

Super_Posterior=function(data,M,EMS=NULL,dilution=0,dil.Par=0.5,Narrative=0,p=0.5,Nar_vec=NULL,FE=0,Time=0,Section=0,Time_FE=0,Section_FE=0,STD=0){

  # data_prep
  data<-data_prep(data,FE=FE,Time=Time,Section=Section,Time_FE=Time_FE,Section_FE=Section_FE,STD=STD)

  # ModelSpace function
  modelSpace<-modelSpace(data,M)

  # Posterior function
  Post<-Posterior(modelSpace,EMS=EMS,dilution=dilution,dil.Par=dil.Par,Narrative=Narrative,p=p,Nar_vec=Nar_vec)

  return(Post)
}# THE END of the Super_Posterior FUNCTION
