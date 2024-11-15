#' Table with the best models according to one of the posterior criteria
#'
#' This function creates a ranking of best models according to one of the possible criterion (PMP under binomial model prior, PMP under binomial-beta model prior, R^2 under binomial model prior, R^2 under binomial-beta model prior).
#' The function gives two types of tables in three different formats: inclusion table (where 1 indicates presence of the regressor in the model and 0 indicates that the variable is excluded from the model) and estimation results table (it displays the best models and estimation output for those models: point estimates, standard errors, significance level, and R^2).
#'
#' @param Post Posterior (Post) object (the result of the Posterior function)
#' @param criterion The criterion that will be used for a basis of the model ranking: \cr
#' 1 - PMP, binomial model prior \cr
#' 2 - PMP, binomial-beta model prior \cr
#' 3 - R2, binomial model prior \cr
#' 4 - R2, binomial-beta model prior
#' @param best The number of the best models to be considered
#' @param estimate A parameter with values TRUE or FALSE indicating which table should be displayed at when the function finishes calculations. \cr
#' TRUE - table with estimation to the results \cr
#' FALSE - table with the inclusion of regressors in the best models \cr
#' @param showConst A parameter with values TRUE or FALSE indicating whether a constant should be included in the output tables. \cr
#' TRUE - the constant is included in the output tables \cr
#' FALSE - the constant is not included in the output tables
#'
#' @return A list with BestModels objects (BS objects): \cr
#' 1. knitr_kable table with inclusion of the regressors in the best models (the best for the display on the console - up to 11 models) \cr
#' 2. knitr_kable table with estimation output in the best models (the best for the display on the console - up to 6 models) \cr
#' 3. matrix with inclusion of the regressors in the best models \cr
#' 4. matrix with estimation output in the best models \cr
#' 5. gTree table with inclusion of the regressors in the best models (displayed as a plot). Use grid::grid.draw() to display.\cr
#' 6. gTree table with estimation output in the best models (displayed as a plot). Use grid::grid.draw() to display.
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
#' mSpace<-modelSpace(data,M)
#' Post<-Posterior(mSpace,dilution=1,dil.Par=0.5)
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
#' mSpace<-modelSpace(data,M)
#' Post<-Posterior(mSpace,dilution=1,dil.Par=0.5)
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
#' mSpace<-modelSpace(data,M)
#' Nar_vec<-as.matrix(c(0,1,1,1,2,2))
#' Post<-Posterior(mSpace,Narrative=0,p=0.5,Nar_vec=Nar_vec)
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
#' mSpace<-modelSpace(data,M)
#' Nar_vec<-as.matrix(c(0,1,1,1,2,2))
#' Post<-Posterior(mSpace,Narrative=0,p=0.5,Nar_vec=Nar_vec)
#'

BestModels=function(Post,criterion=1,best=5,estimate=TRUE,showConst=TRUE){

  #Extraction of the information form the Posterior object (Post object)
  x_names<-Post[[6]]# we extract vector with names of the regressors from Posterior object (Post object)
  x_names<-x_names[-1]# we delete the name of the dependent variable
  M<-Post[[7]]# we extract M - maximum number of regressors in a model from Posterior object (Post object)
  K<-Post[[8]]# we extract K - total number of regressors from Posterior object (Post object)
  MS<-Post[[9]]# we extract MS - size of the mode space from Posterior object (Post object)
  forBestModels<-Post[[12]] # Extraction of all the information

  # Division of the information into the relevant categories
  Reg_ID<-as.matrix(forBestModels[,1:M]) # extraction of the model IDs
  PMPs<-as.matrix(forBestModels[,(M+1):(M+4)]) # extraction of the model weights: PMP_uniform,PMP_random,R2_uniform,R2_random
  betas<-as.matrix(forBestModels[,(M+5):(2*M+5)]) # extraction of the parameter coefficients
  stds<-as.matrix((forBestModels[,(2*M+6):(3*M+6)])^0.5) # extraction of the parameter variances and turning them into standard errors
  R2<-as.matrix(forBestModels[,3*M+7]) # extraction of R^2
  DF<-as.matrix(forBestModels[,3*M+8]) # extraction of degrees of freedom (DF)

  if (best>MS){# CONDITION about what to do if the user set best that is higher than MS
    # we tell the user that we are setting best=MS
    message("best>MS - number of best models cannot be bigger than the total number of models. We set best=MS and continiue :)")
    best=MS # we set best=MS
  }# end of the CONDITION about what to do if the user set best that is higher than MS

  # check for the criterion chosen by the user
  if (criterion==1){ranking<-PMPs[,1]} # PMP,uniform
  if (criterion==2){ranking<-PMPs[,2]} # PMP,random
  if (criterion==3){ranking<-PMPs[,3]} # R2,uniform
  if (criterion==4){ranking<-PMPs[,4]} # R2,random

  Ranking<-cbind(ranking,Reg_ID,betas,stds,R2,DF) # we add ranking criterion based on the users choice

  # here we order the models according to PMP criterion
  Ranking<-Ranking[order(Ranking[,1],decreasing=T),] # ordering of the models
  Best_models<-Ranking[1:best,2:(M+1)] # extraction of the matrix with ordered model IDs
  Ranks<-round(Ranking[1:best,1], digits = 3) # PMPs (R^2s) of the first 'best' models
  bestBetas<-Ranking[1:best,(M+2):(2*M+2)] # extraction of the coefficients
  bestSTDs<-Ranking[1:best,(2*M+3):(3*M+3)] # extraction of the standard errors
  bestR2<-as.matrix(Ranking[1:best,(3*M+4)]) # extraction of the standard errors
  bestDF<-as.matrix(Ranking[1:best,(3*M+5)]) # extraction of the degrees of freedom (DF)

  forRegressors<-matrix(0,nrow=best,ncol=K) # matrix for regressors indices

  for (i in 1:best){# at this LOOP we go through all the models up to "best" chosen by the user
    for (k in 1:M){# at this LOOP we go through all regressors ID
      if (Best_models[i,k]!=0){# CONDITION for the presence of a regressor in a model
        forRegressors[i,Best_models[i,k]]=1 #
      }# the end of the CONDITION for the presence of a regressor in a model
    }# the end of the LOOP at which we go through all regressors ID
  }# the end of the LOOP at which we go through all the models up to "best" chosen by the user

  final_names<-c(x_names,"PMP")
  Rank<-matrix(0,nrow=best,ncol=1)# matrix to store model ranks
  for (i in 1:best){# at this LOOP we go through all the models up to "best" chosen by the user
    Rank[i,1]=paste0("'No. ",i,"'")
  }# the end of the LOOP at which we go through all the models up to "best" chosen by the user
  table<-cbind(forRegressors,Ranks) # we add a column with PMPs to the models
  colnames(table)<-final_names # we add regressor names to the columns
  rownames(table)<-Rank # we add regressor names to the rows
  bestTable<-t(table) # we make transposition to make the table look better
  matbestTable<-bestTable
  bestTable<-as.data.frame(bestTable)

    bestCoefs<-matrix(0,nrow=K,ncol=best)# matrix to store COEFFICIENTS of the best models
    bestSE<-matrix(0,nrow=K,ncol=best)# matrix to store STANDARD ERRORS of the best models
    bestPvalues<-matrix(1,nrow=K,ncol=best)# matrix to store P-VALUES of the best models

    ### CONDITION checking if the user wants to see CONSTANTS in the  estimation output
    if (showConst==TRUE){# CONDITION for the extraction of statistics for constants form the best models

      const_Coefs<-matrix(0,nrow=1,ncol=best)# matrix to store POINT ESTIMATES of the constants
      const_SE<-matrix(0,nrow=1,ncol=best)# matrix to store STANDARD ERRORS of the constants
      const_P_values<-matrix(1,nrow=1,ncol=best)# matrix to store P-VALUES of the constants

      # here we extract all the statistics for the best CONSTANTS
      for (j in 1:best){# at this LOOP we go through all the "best" models
        const_Coefs[1,j]=bestBetas[j,1]# extraction of coefficients of the constants
        const_SE[1,j]=bestSTDs[j,1]# extraction of standard errors of the constants
        const_P_values[1,j]=2*stats::pt(abs(bestBetas[j,1]/bestSTDs[j,1]),df=DF[j,1],lower.tail=FALSE)# extraction of p-values of the constants
      }# the end of the LOOP at which we go through all the "best" models
    }# the end of the CONDITION for the extraction of statistics for constants form the best models

    # here we extract all the statistics for the best REGRESSORS
    for (k in 1:K){# at this LOOP we go though all the regressors
      for (j in 1:best){# at this LOOP we go through all the "best" models
        for (i in 1:M){# at this LOOP we go though all the regressors in a model "j"
          if (Best_models[j,i]==k){# CONDITION checking if a given regressor is in a model
            bestCoefs[k,j]=bestBetas[j,1+i]# extraction of coefficients of the regressors
            bestSE[k,j]=bestSTDs[j,1+i]# extraction of standard errors of the regressors
            bestPvalues[k,j]=2*stats::pt(abs(bestBetas[j,1+i]/bestSTDs[j,1+i]),df=DF[j,1],lower.tail=FALSE)# extraction of p-values of the regressors
          }# the end of the CONDITION checking if a given regressor is in a model
        }# end of the LOOP at which we go though all the regressors in a model "j"
      }# the end of the LOOP at which we go through all the "best" models
    }# the end of the LOOP we go though all the regressors

    const=0

    if (showConst==TRUE){# CONDITION checking if constants should be included in the output - the case of adding the results for constant terms
      bestCoefs<-rbind(const_Coefs,bestCoefs) # adding coefficients for constants and regressors
      bestSE<-rbind(const_SE,bestSE) # adding standard errors for constants and regressors
      bestPvalues<-rbind(const_P_values,bestPvalues)# adding p-values for constants and regressors
      const=1
    }# the end of the CONDITION checking if constants should be included in the output

    ## Here we assign asterisks (*) to coefficients based on p-values
    Asterisks<-matrix(0,nrow=K+const,ncol=best) # matrix to store asterisks (*)

    for (i in 1:(K+const)){# LOOP that goes through all regressors (+ constant)
      for (j in 1:best){# LOOP that goes through all the best models
        # NOT statistically significant, p-value: above 0.1
        if (bestPvalues[i,j]>=0.1){# CONDITION: not statistically significant
          Asterisks[i,j]=""
          # statistically significant, p-value: below 0.1 (*)
        }else if (bestPvalues[i,j]<0.1&bestPvalues[i,j]>=0.05){# CONDITION: statistically significant at 0.9 level
          Asterisks[i,j]="*"
          # statistically significant, p-value: below 0.05 (**)
        }else if (bestPvalues[i,j]<0.05&bestPvalues[i,j]>=0.01){# CONDITION: statistically significant at 0.95 level
          Asterisks[i,j]="**"
          # statistically significant, p-value: below 0.01 (***)
        }else if  (bestPvalues[i,j]<0.01){# CONDITION: statistically significant at 0.99 level
          Asterisks[i,j]="***"}
      }# the end of the LOOP that goes through all the best models
    }# the end of the LOOP that goes through all regressors (+ constant)

    finalTable<-matrix(0,nrow=K+const,ncol=best) # matrix to store final tables with estimation output

    for (i in 1:(K+const)){# LOOP that goes through all regressors (+ constant)
      for (j in 1:best){# LOOP that goes through all the best models
        if (bestCoefs[i,j]==0){# CONDITION for the EXCLUSION of the variable from a model
          finalTable[i,j]=""
        }# the end of the CONDITION for the EXCLUSION of the variable from a model
        else {# CONDITION for the INCLUSION of the variable from a model
          finalTable[i,j]=paste0(round(bestCoefs[i,j],digits=2)," (",round(bestSE[i,j],digits=2),")",Asterisks[i,j])
        }# the end of the CONDITION for the INCLUSION of the variable from a model
      }# the end of the LOOP that goes through all the best models
    }# the end of the LOOP that goes through all regressors (+ constant)

    ## Adding constant to variables if the user chooses to do so
    if (showConst==TRUE){# CONDITION for adding the CONSTANT
      Estimation_names<-c("Constant",x_names) # names of variables + constant at the begining
    }else{ # CONDITION for NOT adding the CONSTANT
      Estimation_names<-x_names # names of variables WITHOUT w constant
    }#the end of the CONDITIONS for EXCLUDING/INCLUDING a CONSTANT

    finalTable<-cbind(Estimation_names,finalTable) # we add regressor names to the columns
    Rank2<-rbind("Model",Rank) # we add "Model" to column names
    finalR2<-cbind("R^2",t(round(bestR2,digits=2)))# Creation of a row with R^2
    finalTable<-rbind(finalTable,finalR2)# we add row with r^2 to the final table
    colnames(finalTable)<-Rank2 # we add model ranks
    matfinalTable<-finalTable

  # ADDITIONAL OUTPUT
  if (estimate==FALSE){# CONDITION for display of a INCLUSION table
    BestTable<-gridExtra::grid.table(bestTable) # here is a more fancy version of the INCLUSION table (works with not many models)
  }# the end of the CONDITION for display of a INCLUSION table
  if (estimate==TRUE){# CONDITION for display of a ESTIMATION table
    FinalTable<-gridExtra::grid.table(finalTable) # here is a more fancy version of the ESTIMATION table (works with not many models)
  }# the end of the CONDITION for display of a ESTIMATION table
  FinalTable<-grid::grid.grabExpr(gridExtra::grid.table(finalTable))
  BestTable<-grid::grid.grabExpr(gridExtra::grid.table(bestTable))
    finalTable<-as.data.frame(finalTable)
  bestTable<-knitr::kable(bestTable, row.names = TRUE,align = "c")
  finalTable<-knitr::kable(finalTable, row.names = FALSE,align = "c")
  out<-list(bestTable,finalTable,matbestTable,matfinalTable,BestTable,FinalTable)
  return(out)
}
