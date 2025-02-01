#' Graphs of the prior and posterior model probabilities for the best individual models
#'
#' This function draws four graphs of prior and posterior model probabilities for the best individual models: \cr
#' a) The results with binomial model prior (based on PMP - posterior model probability) \cr
#' b) The results with binomial-beta model prior (based on PMP - posterior model probability) \cr
#' c) The results with binomial model prior based on R^2 \cr
#' d) The results with binomial-beta model prior based on R^2 \cr
#' Models on the graph are ordered according to their posterior model probability.
#'
#'
#' @param Post Posterior (Post) object (the result of the Posterior function)
#' @param Top The number of the best model to be placed on the graphs
#'
#' @return Four graphs with prior and posterior model probabilities for individual models:\cr
#' 1) The results with binomial model prior (based on PMP - posterior model probability) \cr
#' 2) The results with binomial-beta model prior (based on PMP - posterior model probability) \cr
#' 3) The results with binomial model prior based on R^2 \cr
#' 4) The results with binomial-beta model prior based on R^2 \cr
#' 5) On graph combining all the four aforementioned graphs
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
#' PMPgraphs<-modelPMP(Post)
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
#' PMPgraphs<-modelPMP(Post,Top=7)
#'
#'@name modelPMP

utils::globalVariables(c("ID", "Value", "Probability"))

modelPMP=function(Post,Top=4){

  K<-Post[[8]] # Total number of regressors
  MS<-Post[[9]] # Number of models in the models space

  if (Top>MS){# CONDITION about what to do if the user sets Top that is higher than MS (Top>MS)
    # we tell the user that we are setting Top=K
    message("The number of the best models (Top) cannot be higher than the total number of models. We set Top=4 (total number of regressors) and continiue :)")
    MS=K # we set M=K
  }# end of the CONDITION about what to do if the user set M that is higher than K (M>K)

  # Collecting information from the Post object
  M<-Post[[7]] # Maximum number of regressors in a model
  MS<-Post[[9]] # Number of models in the models space
  PMPs<-Post[[12]][,(M+1):(M+4)] # PMP_uniform,PMP_random,PMP_R2_uniform,PMP_R2_random
  Priors<-Post[[17]] # Priors: uniform and random

  # Objects to store posteriors and priors
  PMP_uniform<-cbind(PMPs[,1],Priors[,1])
  PMP_random<-cbind(PMPs[,2],Priors[,2])
  R2_uniform<-cbind(PMPs[,3],Priors[,1])
  R2_random<-cbind(PMPs[,4],Priors[,2])

  # Ordering of the models according to posterior criterion
  PMP_uniform<-PMP_uniform[order(PMP_uniform[,1],decreasing=T),]
  PMP_random<-PMP_random[order(PMP_random[,1],decreasing=T),]
  R2_uniform<-R2_uniform[order(R2_uniform[,1],decreasing=T),]
  R2_random<-R2_random[order(R2_random[,1],decreasing=T),]

  ranking<-matrix(1:MS,nrow=MS,ncol=1)

  # Adding a ranking number
  PMP_uniform<-cbind(ranking[1:Top,],PMP_uniform[1:Top,])
  PMP_random<-cbind(ranking[1:Top,],PMP_random[1:Top,])
  R2_uniform<-cbind(ranking[1:Top,],R2_uniform[1:Top,])
  R2_random<-cbind(ranking[1:Top,],R2_random[1:Top,])

  IDnames<-cbind("ID","Posterior","Prior") #creating names of the variables to be used by 'tidyverse' package

  colnames(PMP_uniform)<-IDnames # we add names to the columns
  colnames(PMP_random)<-IDnames # we add names to the columns
  colnames(R2_uniform)<-IDnames # we add names to the columns
  colnames(R2_random)<-IDnames # we add names to the columns

  forGraph1<-as.data.frame(PMP_uniform) # changing the table to data frame for ggplot
  forGraph2<-as.data.frame(PMP_random) # changing the table to data frame for ggplot
  forGraph3<-as.data.frame(R2_uniform) # changing the table to data frame for ggplot
  forGraph4<-as.data.frame(R2_random) # changing the table to data frame for ggplot

  ## Preparation of the Figures with ggplot
  # for graph with uniform model prior and likelihood based posterior
  forGraph1 <-tidyr::gather(forGraph1,key = "Probability", value = "Value", -ID)
  # for graph with random model prior and likelihood based posterior
  forGraph2 <-tidyr::gather(forGraph2,key = "Probability", value = "Value", -ID)
  # for graph with uniform model prior and R^2 based posterior
  forGraph3 <-tidyr::gather(forGraph3,key = "Probability", value = "Value", -ID)
  # for graph with random model prior and R^2 based posterior
  forGraph4 <-tidyr::gather(forGraph4,key = "Probability", value = "Value", -ID)

  ## Preparation of the Figures with ggplot
  # for graph with uniform model prior and likelihood based posterior
  Graph1<-ggplot2::ggplot(forGraph1, ggplot2::aes(x = ID, y = Value)) +
    ggplot2::geom_line(ggplot2::aes(color = Probability, linetype = Probability)) +
    ggplot2::scale_color_manual(values = c("darkred", "steelblue"))+
    ggplot2::ylab("Prior, Posterior") + ggplot2::xlab("Model number in the raniking")
  # for graph with random model prior and likelihood based posterior
  Graph2<-ggplot2::ggplot(forGraph2, ggplot2::aes(x = ID, y = Value)) +
    ggplot2::geom_line(ggplot2::aes(color = Probability, linetype = Probability)) +
    ggplot2::scale_color_manual(values = c("darkred", "steelblue"))+
    ggplot2::ylab("Prior, Posterior") + ggplot2::xlab("Model number in the raniking")
  # for graph with uniform model prior and R^2 based posterior
  Graph3<-ggplot2::ggplot(forGraph3, ggplot2::aes(x = ID, y = Value)) +
    ggplot2::geom_line(ggplot2::aes(color = Probability, linetype = Probability)) +
    ggplot2::scale_color_manual(values = c("darkred", "steelblue"))+
    ggplot2::ylab("Prior, Posterior") + ggplot2::xlab("Model number in the raniking")
  # for graph with random model prior and R^2 based posterior
  Graph4<-ggplot2::ggplot(forGraph4, ggplot2::aes(x = ID, y = Value)) +
    ggplot2::geom_line(ggplot2::aes(color = Probability, linetype = Probability)) +
    ggplot2::scale_color_manual(values = c("darkred", "steelblue"))+
    ggplot2::ylab("Prior, Posterior") + ggplot2::xlab("Model number in the raniking")

  ## Preparation of the data for BIG COMBINED GRAPH
  # for graph with uniform model prior and likelihood based posterior
  Graph1_2<-ggplot2::ggplot(forGraph1, ggplot2::aes(x = ID, y = Value)) +
    ggplot2::geom_line(ggplot2::aes(color = Probability, linetype = Probability)) +
    ggplot2::scale_color_manual(values = c("darkred", "steelblue"))+
    ggplot2::ylab("Prior, Posterior") + ggplot2::xlab("Model number in the ranking")+ggplot2::ggtitle("Results with binomial model prior")
  # for graph with random model prior and likelihood based posterior
  Graph2_2<-ggplot2::ggplot(forGraph2, ggplot2::aes(x = ID, y = Value)) +
    ggplot2::geom_line(ggplot2::aes(color = Probability, linetype = Probability)) +
    ggplot2::scale_color_manual(values = c("darkred", "steelblue"))+
    ggplot2::ylab("Prior, Posterior") + ggplot2::xlab("Model number in the ranking")+ggplot2::ggtitle("Results with binomial-beta model prior")
  # for graph with unifrom model prior and R^2 based posterior
  Graph3_2<-ggplot2::ggplot(forGraph3, ggplot2::aes(x = ID, y = Value)) +
    ggplot2::geom_line(ggplot2::aes(color = Probability, linetype = Probability)) +
    ggplot2::scale_color_manual(values = c("darkred", "steelblue"))+
    ggplot2::ylab("Prior, Posterior") + ggplot2::xlab("Model number in the ranking")+ggplot2::ggtitle(bquote("Results with bimonial model prior and" ~ R^2))
  # for graph with random model prior and R^2 based posterior
  Graph4_2<-ggplot2::ggplot(forGraph4, ggplot2::aes(x = ID, y = Value)) +
    ggplot2::geom_line(ggplot2::aes(color = Probability, linetype = Probability)) +
    ggplot2::scale_color_manual(values = c("darkred", "steelblue"))+
    ggplot2::ylab("Prior, Posterior") + ggplot2::xlab("Model number in the ranking")+ggplot2::ggtitle(bquote("Results with binomial-beta model prior and" ~ R^2))

  # Putting together the last plot
  Finalplot<-ggpubr::ggarrange(Graph1_2,Graph2_2,Graph3_2,Graph4_2,
                               labels = c("a)", "b)", "c)","d)"),
                               ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")

  print(Finalplot)
  # creation of the modelSizes object (BS object)
  out<-list(Graph1,Graph2,Graph3,Graph4,Finalplot)
  return(out)
}
