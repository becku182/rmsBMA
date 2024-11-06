#' Introduction of time and section fixed effects and data standardization.
#'
#' If the data is in the panel form we assume it has the following structure\cr
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
#' @param data A data file.
#' @param FE Binary variable: 1 - include fixed effect, 0 - do not include fixed effects.
#' @param Time The number of time periods - works only if FE=1.
#' @param Section The number of cross-sections - works only if EF=1.
#' @param Time_FE Binary variable: 1 - include time fixed effect, 0 - do not include time fixed effects. Works only if EF=1.
#' @param Section_FE Binary variable: 1 - include cross-section fixed effect, 0 - do not include cross-section fixed effects. Works only if EF=1.
#' @param STD Binary variable: 1 - standardize the data set, 0 - do not standardize the data set. By standardization we mean subracion of amean and division  by standard deviation of each variable.
#'
#' @return Formatted data set.
#' @export
#'
#' @examples
#' y <- matrix(1:20,nrow=20,ncol=1)
#' x1 <- matrix(21:40,nrow=20,ncol=1)
#' x2 <- matrix(41:60,nrow=20,ncol=1)
#' data <- cbind(y,x1,x2)
#' new_data <- data_prep(data,FE=1,Time=5,Section=4,Time_FE=1,Section_FE=1,STD=0)
#'
#' y <- rnorm(20, mean = 0, sd = 1)
#' x1 <- rnorm(20, mean = 0, sd = 1)
#' x2 <- rnorm(20, mean = 0, sd = 1)
#' data <- cbind(y,x1,x2)
#' new_data <- data_prep(data,FE=1,Time=5,Section=4,Time_FE=1,Section_FE=1,STD=1)
#'

data_prep <- function(data,FE=0,Time=0,Section=0,Time_FE=0,Section_FE=0,STD=0){# BEGINING OF THE data_prep function

  Var_names <- names(data) # names of the variables
  colnames(data) <- NULL
  data <- as.matrix(data)

  # collecting data characteristics
  m <- nrow(data) # number of rows in the data
  n <- ncol(data) # number of columns in the data

  #### CONDITION for FIXED EFFECTS
  if (FE==1){# CONDITION checking if the user wants to use fixed effects

    ##### CONDITION checking if the user specified Section/Time fixed effects
    if (Section_FE==0&Time_FE==0){# CONDITION checking if the user specified CROSS-SECTION and TIME fixed effects
      stop("Please specify if you want to use Cross-section (Section_FE) or/and Time (Time_FE) fixed effects. If you do NOT want to use fixed effects please set FE=0")
    }# the end of the CONDITION checking if the user specified CROSS-SECTION and TIME fixed effects

    ###### CONDITION that checks if "Section" and "Time" provied by the user match the number of observations in the dataset
    if (m!=Time*Section){# CONDITION about what to do if the total number of observations is different than
      # product of the number of cross-sections (Section) and the number of time periods (Time)
      # if the CONDITION is met the function stops and provides
      # a user with a messege on the reason why the function stops
      stop("total number of observations in not equal to the product of cross-sections and periods (Section*Time)")
    }# the end of the CONDITION about what to do if the total number of observations

    Section_ID <- kronecker(matrix(1,nrow=Time,ncol=1),matrix(1:Section,nrow=Section,ncol=1))
    Time_ID <- kronecker(matrix(1:Time,nrow=Time,ncol=1),matrix(1,nrow=Section,ncol=1))
    ID <- cbind(Section_ID,Time_ID)

    FEdata <- cbind(ID,data)
    FEdata

    ########### TIME FIXED EFFECTS
    if (Time_FE==1){# CONDITION that checks if the user wanted TIME fixed effects
      For_TFE <- FEdata[order(FEdata[, 2]),]
      For_TFE2 <- For_TFE[,3:(n+2)]
      TFE_ID <- For_TFE[,1:2]
      For_S_means <- diag(m)-kronecker(matrix(1,nrow=Time,ncol=Time),(1/Time)*diag(Section))
      TFE <- round(For_S_means%*%For_TFE2,11)
      After <- cbind(TFE_ID,TFE)
      FEdata <- After[order(After[, 2]),]
    }# the end of the CONDITION that checks if the user wanted TIME fixed effects

    ########### CROSS-SECTION FIXED EFFECTS (covers the case of TIME and COUNTRY fixed effects)
    if (Section_FE==1){# CONDITION that checks if the user wanted SECTION fixed effects
      For_SFE <- FEdata[order(FEdata[, 1]),]
      For_SFE2 <- For_SFE[,3:(n+2)]
      SFE_ID <- For_SFE[,1:2]
      For_T_means <- diag(m)-kronecker(matrix(1,nrow=Section,ncol=Section),(1/Section)*diag(Time))
      SFE <- round(For_T_means%*%For_SFE2,11)
      After <- cbind(SFE_ID,SFE)
      FEdata <- After[order(After[, 2]),]
    }# the end of the CONDITION that checks if the user wanted CROSS-SECTION fixed effects

    data <- round(FEdata[,3:(n+2)],11)

  }# the end of the CONDITION checking if the user wants to use fixed effects



  ########## STANDARDIZATION OF THE DATA
  if (STD==1){# CONDITION checking if the user wants to STANDARDIZE the data

    # TEST IF data is not a matrix of zeros
    if (all(data == 0)){
      stop("Fixed effects left the matrix of zeros. Standardization cannot be perforemd")
    }

    STDmeans <- apply(data, 2, mean) # calculation of the column means
    STDstds <- apply(data, 2, stats::sd) # calculation of the column standard deviations
    STDdata <- matrix(0,nrow=m,ncol=n) # Matrix to store standardized data

    for (i in 1:n){# at this LOOP we go though all the variables
      STDdata[1:m,i]=(data[1:m,i]-STDmeans[i])/STDstds[i]
    }# the end of the LOOP at which we go though all the variables

    # Here we provide the data for further calculations
    data <- STDdata # standardized data (STDdata) is now a data file for further use

  }# the end of the CONDITION checking if the user wants to STANDARDIZE the data

  colnames(data) <- Var_names
  return(data)

}# THE END OF THE data_prep function
