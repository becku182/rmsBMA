#' Introduction of time and section fixed effects and data standardization.
#'
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
#' @param data A data file.
#' @param FE Binary variable: TRUE - include fixed effect, FALSE - do not include fixed effects.
#' @param Time The number of time periods - works only if FE=1.
#' @param Section The number of cross-sections - works only if EF=1.
#' @param Time_FE Binary variable: 1 - include time fixed effect, 0 - do not include time fixed effects. Works only if EF=1.
#' @param Section_FE Binary variable: 1 - include cross-section fixed effect, 0 - do not include cross-section fixed effects. Works only if EF=1.
#' @param STD Binary variable: 1 - standardize the data set, 0 - do not standardize the data set. By standardization we mean subtraction of a mean and division  by standard deviation of each variable.
#'
#' @return Formatted data set.
#'
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

data_prep <- function(data,FE=FALSE,Time=0,Section=0,Time_FE=0,Section_FE=0,STD=0){

  Var_names <- colnames(data)
  colnames(data) <- NULL
  data <- as.matrix(data)

  m <- nrow(data)
  n <- ncol(data)

  if (FE==TRUE){

    if (Section_FE==0 & Time_FE==0){
      stop("Please specify if you want to use Cross-section (Section_FE) or/and Time (Time_FE) fixed effects. If you do NOT want to use fixed effects please set FE=0")
    }

    if (m != Time*Section){
      stop("total number of observations in not equal to the product of cross-sections and periods (Section*Time)")
    }

    Section_ID <- kronecker(matrix(1,nrow=Time,ncol=1),matrix(1:Section,nrow=Section,ncol=1))
    Time_ID <- kronecker(matrix(1:Time,nrow=Time,ncol=1),matrix(1,nrow=Section,ncol=1))
    ID <- cbind(Section_ID,Time_ID)

    FEdata <- cbind(ID,data)
    # removed: FEdata

    if (Section_FE==1){
      For_TFE <- FEdata[order(FEdata[, 2]),]
      For_TFE2 <- For_TFE[,3:(n+2)]
      TFE_ID <- For_TFE[,1:2]
      For_S_means <- diag(m)-kronecker(matrix(1,nrow=Time,ncol=Time),(1/Time)*diag(Section))
      For_S_means <- Matrix::Matrix(For_S_means, sparse=TRUE)
      TFE <- as.matrix(round(For_S_means %*% For_TFE2, 11))  # <- force base matrix
      After <- cbind(TFE_ID, TFE)
      FEdata <- After[order(After[, 2]),]
    }

    if (Time_FE==1){
      For_SFE <- FEdata[order(FEdata[, 1]),]
      For_SFE2 <- For_SFE[,3:(n+2)]
      SFE_ID <- For_SFE[,1:2]
      For_T_means <- diag(m)-kronecker(matrix(1,nrow=Section,ncol=Section),(1/Section)*diag(Time))
      For_T_means <- Matrix::Matrix(For_T_means, sparse=TRUE)
      SFE <- as.matrix(round(For_T_means %*% For_SFE2, 11))  # <- force base matrix
      After <- cbind(SFE_ID, SFE)
      FEdata <- After[order(After[, 2]),]
    }

    data <- round(as.matrix(FEdata[,3:(n+2)]), 11)  # <- ensure base matrix here too
  }

  if (STD==1){

    if (all(data == 0)){
      stop("Fixed effects left the matrix of zeros. Standardization cannot be perforemd")
    }

    STDmeans <- apply(data, 2, mean)
    STDstds <- apply(data, 2, stats::sd)
    STDdata <- matrix(0,nrow=m,ncol=n)

    for (i in 1:n){
      STDdata[1:m,i] <- (data[1:m,i]-STDmeans[i]) / STDstds[i]
    }

    data <- STDdata
  }

  colnames(data) <- Var_names
  return(data)
}

