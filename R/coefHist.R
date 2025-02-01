#' Graphs of the distribution of the coefficients over the model space
#'
#' This function draws graphs of the distribution (in the form of histogram or kernel density) of the coefficients for all the considered regressors over the part of the model space that includes this regressors (half of the model space).
#'
#' @name coefHist
#'
#' @param Post Posterior (Post) object (the result of the Posterior function)
#' @param BW Parameter indicating what method should be chosen to find bin widths for the histograms: \cr
#' 1) "FD" Freedman-Diaconis method \cr
#' 2) "SC" Scott method \cr
#' 3) "vec" user specified bin widths provided through a vector (parameter: binW)
#' @param binW A vector with bin widths to be used to construct histograms for the regressors. The vector must be of the size equal to total number of regressors. The vector with bin widths is used only if parameter BW="vec".
#' @param BN Parameter taking the values(default: BN=0): \cr
#' 1 - the histogram will be build based on the number of bins specified by the user through parameter num. If BN=1, the function ignores parameters BW. \cr
#' 0 - the histogram will be build in line with parameter BW
#' @param num A vector with the numbers of bins used to be used to construct histograms for the regressors. The vector must be of the size equal to total number of regressors. The vector with bin widths is used only if parameter BN=1.
#' @param kernel A parameter taking the values (default: kernel=0):\cr
#' 1 - the function will build graphs using kernel density for the distribution of coefficients (with kernel=1, the function ignores parameters BW and BN) \cr
#' 0 - the function will build regular histogram density for the distribution of coefficients
#'
#' @return A list with the graphs of the distribution of coefficients for all the considered regressors.
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
#' histPlots<-coefHist(Post)
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
#' histPlots<-coefHist(Post)
#'

utils::globalVariables(".data")

coefHist=function(Post,BW="FD",binW=NULL,BN=0,num=NULL,kernel=0){
  x_names<-Post[[6]] # names of the regressors
  x_names<-x_names[-1]
  K<-Post[[8]] # number of regressors
  MS<-Post[[9]] # total number of models
  betas<-Post[[16]] # tables with coefficients on all regressors for all the models


  # Adding colnames and changing to dataframe
  colnames(betas)<-x_names
  betas<-as.data.frame(betas)

  histPlots<-list() # Opening a list  for the histogram plots

  # CONDITION for using kernel or regular histogram
  if (kernel==1){
    for (i in 1:K){
      histPlots[[i]]<-invisible(ggplot2::ggplot(betas, ggplot2::aes(x = .data[[x_names[i]]])) +
                                  ggplot2::geom_density(fill = "skyblue", alpha = 0.7) +
                                  ggplot2::labs(
                                    title = paste("Distribiution of", x_names[i], "coefficients"),
                                    x = paste0("Coefficients on ",x_names[i]),
                                    y = "Frequency") +
                                  ggplot2::theme_minimal(base_size = 12) +
                                  ggplot2::theme(
                                    plot.title = ggplot2::element_text(size = 12, hjust = 0.5, face = "bold"),
                                    axis.title = ggplot2::element_text(face = "bold")))
      names(histPlots)[[i]] <-x_names[[i]]
    }

  }else{# REGULAR HISTOGRAM BELOW

    # CONDITION for graphs plotted with binwidth:
    if (BN==0){### Rules for bin width
      for (i in 1:K){ # at this LOOP we go through all the regressors
        # 1) Freedman-Diaconis (FD)
        if (BW=="FD"){BW<-(stats::IQR(betas[,i])*2)/sqrt(length(betas[,i]))}
        # 2) Scott (SC)
        if (BW=="SC"){BW<-(stats::sd(betas[,i])*3.5)/(length(betas[,i])^(1/3))}
        # 3) Binwidth sizes
        if(BW=="vec"){
          if (is.null(binW)){stop("Please provide a vector with bin width sizes through parameter binW")}
          if (length(binW)!=K){stop("binW is missspecified: binW should have K elements")}
          BW<-binW[i]
        }
        # Building the plot
        histPlots[[i]]<-invisible(ggplot2::ggplot(betas, ggplot2::aes(x = .data[[x_names[i]]])) +
                                    ggplot2::geom_histogram(binwidth=BW, fill = "skyblue", color = "white", alpha = 0.8) +
                                    ggplot2::labs(
                                      title = paste("Distribiution of", x_names[i], "coefficients"),
                                      x = paste0("Coefficients on ",x_names[i]),
                                      y = "Frequency") +
                                    ggplot2::theme_minimal(base_size = 12) +
                                    ggplot2::theme(
                                      plot.title = ggplot2::element_text(size = 12, hjust = 0.5, face = "bold"),
                                      axis.title = ggplot2::element_text(face = "bold")))
        names(histPlots)[[i]] <-x_names[[i]]
      }
    }

    # CONDITION for graphs plotted with bins - through setting the number of bins:
    if (BN==1){### Rules for bin width
      for (i in 1:K){ # at this LOOP we go through all the regressors
        if (is.null(num)){stop("Please provide a vector with number of bins through parameter num")}
        if (length(num)!=K){stop("num is missspecified: num should have K elements")}

        # Building the plot
        histPlots[[i]]<-invisible(ggplot2::ggplot(betas, ggplot2::aes(x = .data[[x_names[i]]])) +
                                    ggplot2::geom_histogram(bins=num[i], fill = "skyblue", color = "white", alpha = 0.8) +
                                    ggplot2::labs(
                                      title = paste("Distribiution of", x_names[i], "coefficients"),
                                      x = paste0("Coefficients on ",x_names[i]),
                                      y = "Frequency") +
                                    ggplot2::theme_minimal(base_size = 12) +
                                    ggplot2::theme(
                                      plot.title = ggplot2::element_text(size = 12, hjust = 0.5, face = "bold"),
                                      axis.title = ggplot2::element_text(face = "bold")))
        names(histPlots)[[i]] <-x_names[[i]]
      }
    }
  }
  return(histPlots)
}
