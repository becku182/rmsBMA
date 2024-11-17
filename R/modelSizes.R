#' Graphs of the prior and posterior model probabilities
#'
#' This function draws four graphs of prior and posterior model probabilities: \cr
#' a) The results with binomial model prior (based on PMP - posterior model probability) \cr
#' b) The results with binomial-beta model prior (based on PMP - posterior model probability) \cr
#' c) The results with binomial model prior based on R^2 \cr
#' d) The results with binomial-beta model prior based on R^2
#'
#' @param Post Posterior (Post) object (the result of the Posterior function)
#'
#' @return Four graphs with prior and posterior model probabilities:\cr
#' 1) The results with binomial model prior (based on PMP - posterior model probability) \cr
#' 2) The results with binomial-beta model prior (based on PMP - posterior model probability) \cr
#' 3) The results with binomial model prior based on R^2 \cr
#' 4) The results with binomial-beta model prior based on R^2 \cr
#' 5) On graph combining all the four aformentioned graphs
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
#' mSizes<-modelSizes(Post)
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
#' mSizes<-modelSizes(Post)
#'
#'@name modelSizes

utils::globalVariables(c("ID", "Value", "Probability"))

modelSizes=function(Post){

  M<-Post[[7]]# we extract M - maximum number of regressors in a model from Posterior object (Post object)
  K<-Post[[8]]# we extract K - total number of regressors from Posterior object (Post object)
  MS<-Post[[9]]# we extract MS - size of the mode space from Posterior object (Post object)
  sizePriors<-Post[[13]]# we extract sizePriors - table with unifrom and random model priors spread over model sizes from Posterior object (Post object)
  modelPosterior<-Post[[14]]# we extract modelPosterior - table with posterior model probabilities from Posterior object (Post object)

  sizes<-matrix(0,nrow=M+1,ncol=1) #we create vector to store number of models in a given model size

  for (k in 0:M){# at this LOOP we add all combinations of regressors up models with M variables
    sizes[k+1,1]<-choose(K,k) # number of models of the size k out of K regressors
  } # this sum adds up all the models for each possible model size

  ind<-cumsum(sizes) # we create a vector with the number of models in each model size category

  Posterior_sizes<-matrix(0,nrow=M,ncol=4) # matrix to store posterior probabilities over model sizes

  No_regressors<-matrix(0,nrow=1,ncol=4)
  No_regressors[1,1:4]=modelPosterior[1,1:4]  # insertion of posteriors for model with no regressors

  for (i in 1:M){# at this LOOP we go through all the model sizes
    Posterior_sizes[i,1]=sum(modelPosterior[(ind[i]+1):(ind[i+1]),1])# posterior under uniform model prior
    Posterior_sizes[i,2]=sum(modelPosterior[(ind[i]+1):(ind[i+1]),2])# posterior under random model prior
    Posterior_sizes[i,3]=sum(modelPosterior[(ind[i]+1):(ind[i+1]),3])# R2 posterior under uniform model prior
    Posterior_sizes[i,4]=sum(modelPosterior[(ind[i]+1):(ind[i+1]),4])# R2 posterior under random model prior
  }# the end of the LOOP at which we go through all the model sizes

  Posterior_sizes<-rbind(No_regressors,Posterior_sizes) # insertions of posteriors for model with no regressors

  # Preparation of the tables for graphs
  forGraph1<-cbind(0:M,sizePriors[,1],Posterior_sizes[,1])# for graph with uniform model prior and likelihood based posterior
  forGraph2<-cbind(0:M,sizePriors[,2],Posterior_sizes[,2])# for graph with random model prior and likelihood based posterior
  forGraph3<-cbind(0:M,sizePriors[,1],Posterior_sizes[,3])# for graph with uniform model prior and R^2 based posterior
  forGraph4<-cbind(0:M,sizePriors[,2],Posterior_sizes[,4])# for graph with random model prior and R^2 based posterior

  IDnames<-cbind("ID","Prior","Posterior") #creating names of the variables to be used by 'tidyverse' package

  colnames(forGraph1)<-IDnames # we add names to the columns
  colnames(forGraph2)<-IDnames # we add names to the columns
  colnames(forGraph3)<-IDnames # we add names to the columns
  colnames(forGraph4)<-IDnames # we add names to the columns

  forGraph1<-as.data.frame(forGraph1) # changing the table to data frame for ggplot
  forGraph2<-as.data.frame(forGraph2) # changing the table to data frame for ggplot
  forGraph3<-as.data.frame(forGraph3) # changing the table to data frame for ggplot
  forGraph4<-as.data.frame(forGraph4) # changing the table to data frame for ggplot

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
    ggplot2::ylab("Prior, Posterior") + ggplot2::xlab("Model size (number of regressors)")
  # for graph with random model prior and likelihood based posterior
  Graph2<-ggplot2::ggplot(forGraph2, ggplot2::aes(x = ID, y = Value)) +
    ggplot2::geom_line(ggplot2::aes(color = Probability, linetype = Probability)) +
    ggplot2::scale_color_manual(values = c("darkred", "steelblue"))+
    ggplot2::ylab("Prior, Posterior") + ggplot2::xlab("Model size (number of regressors)")
  # for graph with uniform model prior and R^2 based posterior
  Graph3<-ggplot2::ggplot(forGraph3, ggplot2::aes(x = ID, y = Value)) +
    ggplot2::geom_line(ggplot2::aes(color = Probability, linetype = Probability)) +
    ggplot2::scale_color_manual(values = c("darkred", "steelblue"))+
    ggplot2::ylab("Prior, Posterior") + ggplot2::xlab("Model size (number of regressors)")
  # for graph with random model prior and R^2 based posterior
  Graph4<-ggplot2::ggplot(forGraph4, ggplot2::aes(x = ID, y = Value)) +
    ggplot2::geom_line(ggplot2::aes(color = Probability, linetype = Probability)) +
    ggplot2::scale_color_manual(values = c("darkred", "steelblue"))+
    ggplot2::ylab("Prior, Posterior") + ggplot2::xlab("Model size (number of regressors)")

  ## Preparation of the data for BIG COMBINED GRAPH
  # for graph with uniform model prior and likelihood based posterior
  Graph1_2<-ggplot2::ggplot(forGraph1, ggplot2::aes(x = ID, y = Value)) +
    ggplot2::geom_line(ggplot2::aes(color = Probability, linetype = Probability)) +
    ggplot2::scale_color_manual(values = c("darkred", "steelblue"))+
    ggplot2::ylab("Prior, Posterior") + ggplot2::xlab("Model size (number of regressors)")+ggplot2::ggtitle("Results with binomial model prior")
  # for graph with random model prior and likelihood based posterior
  Graph2_2<-ggplot2::ggplot(forGraph2, ggplot2::aes(x = ID, y = Value)) +
    ggplot2::geom_line(ggplot2::aes(color = Probability, linetype = Probability)) +
    ggplot2::scale_color_manual(values = c("darkred", "steelblue"))+
    ggplot2::ylab("Prior, Posterior") + ggplot2::xlab("Model size (number of regressors)")+ggplot2::ggtitle("Results with binomial-beta model prior")
  # for graph with unifrom model prior and R^2 based posterior
  Graph3_2<-ggplot2::ggplot(forGraph3, ggplot2::aes(x = ID, y = Value)) +
    ggplot2::geom_line(ggplot2::aes(color = Probability, linetype = Probability)) +
    ggplot2::scale_color_manual(values = c("darkred", "steelblue"))+
    ggplot2::ylab("Prior, Posterior") + ggplot2::xlab("Model size (number of regressors)")+ggplot2::ggtitle(bquote("Results with bimonial model prior and" ~ R^2))
  # for graph with random model prior and R^2 based posterior
  Graph4_2<-ggplot2::ggplot(forGraph4, ggplot2::aes(x = ID, y = Value)) +
    ggplot2::geom_line(ggplot2::aes(color = Probability, linetype = Probability)) +
    ggplot2::scale_color_manual(values = c("darkred", "steelblue"))+
    ggplot2::ylab("Prior, Posterior") + ggplot2::xlab("Model size (number of regressors)")+ggplot2::ggtitle(bquote("Results with binomial-beta model prior and" ~ R^2))

  # Putting together the last plot
  Finalplot<-ggpubr::ggarrange(Graph1_2,Graph2_2,Graph3_2,Graph4_2,
                               labels = c("a)", "b)", "c)","d)"),
                               ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")

  print(Finalplot)
  # creation of the modelSizes object (BS object)
  out<-list(Graph1,Graph2,Graph3,Graph4,Finalplot)
}
