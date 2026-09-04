#' Graphs of the prior and posterior model probabilities of the model sizes
#'
#' This function draws four graphs of prior and posterior model probabilities: \cr
#' a) The results with binomial model prior (based on PMP - posterior model probability) \cr
#' b) The results with binomial-beta model prior (based on PMP - posterior model probability)
#'
#' @param bma_list bma_list object (the result of the bma function)
#'
#' @return A list with three graphs with prior and posterior model probabilities for model sizes:\cr
#' 1) The results with binomial model prior (based on PMP - posterior model probability) \cr
#' 2) The results with binomial-beta model prior (based on PMP - posterior model probability) \cr
#' 3) One graph combining all the aforementioned graphs
#'
#' @export
#'
#' @examples
#' \donttest{
#' data("Trade_data", package = "rmsBMA")
#' data <- Trade_data[,1:10]
#' modelSpace <- model_space(data, M = 9, g = "UIP")
#' bma_list <- bma(modelSpace)
#' sizes <- model_sizes(bma_list)
#' sizes[[1]]
#' }
#'
#' @name model_sizes

utils::globalVariables(c("ID", "Value", "Probability"))

model_sizes <- function(bma_list){

  EMS <- bma_list[[8]] # expected model size
  uniform_posterior <- bma_list[[13]][,3]
  M <- length(uniform_posterior)
  random_posterior <- bma_list[[13]][,4]
  uniform_prior <- bma_list[[13]][,1]
  random_prior <- bma_list[[13]][,2]
  dilution <- bma_list[[9]] # 0 - no dilution prior, 1 - dilution prior

  # Preparation of the tables for graphs
  forGraph1 <- cbind(0:(M-1), uniform_prior, uniform_posterior)
  forGraph2 <- cbind(0:(M-1), random_prior, random_posterior)

  IDnames <- cbind("ID", "Prior", "Posterior") # names of the variables to be used by 'tidyverse'

  colnames(forGraph1) <- IDnames
  colnames(forGraph2) <- IDnames

  forGraph1 <- as.data.frame(forGraph1)
  forGraph2 <- as.data.frame(forGraph2)

  ## Preparation of the Figures with ggplot
  forGraph1 <- tidyr::gather(forGraph1, key = "Probability", value = "Value", -ID)
  forGraph2 <- tidyr::gather(forGraph2, key = "Probability", value = "Value", -ID)

  ## Preparation of the Figures with ggplot

  Graph1 <- ggplot2::ggplot(forGraph1, ggplot2::aes(x = ID, y = Value)) +
    ggplot2::geom_line(ggplot2::aes(color = Probability, linetype = Probability)) +
    ggplot2::scale_color_manual(values = c("darkred", "steelblue")) +
    ggplot2::ylab("Prior, Posterior") + ggplot2::xlab("Model size (number of regressors)")

  Graph2 <- ggplot2::ggplot(forGraph2, ggplot2::aes(x = ID, y = Value)) +
    ggplot2::geom_line(ggplot2::aes(color = Probability, linetype = Probability)) +
    ggplot2::scale_color_manual(values = c("darkred", "steelblue")) +
    ggplot2::ylab("Prior, Posterior") + ggplot2::xlab("Model size (number of regressors)")

  ## Preparation of the data for BIG COMBINED GRAPH
  if (dilution==0){
    Graph1_2 <- ggplot2::ggplot(forGraph1, ggplot2::aes(x = ID, y = Value)) +
      ggplot2::geom_line(ggplot2::aes(color = Probability, linetype = Probability)) +
      ggplot2::scale_color_manual(values = c("darkred", "steelblue")) +
      ggplot2::ylab("Prior, Posterior") + ggplot2::xlab("Model size (number of regressors)") +
      ggplot2::ggtitle(paste0("Results with binomial model prior (EMS = ", EMS, ")"))

    Graph2_2 <- ggplot2::ggplot(forGraph2, ggplot2::aes(x = ID, y = Value)) +
      ggplot2::geom_line(ggplot2::aes(color = Probability, linetype = Probability)) +
      ggplot2::scale_color_manual(values = c("darkred", "steelblue")) +
      ggplot2::ylab("Prior, Posterior") + ggplot2::xlab("Model size (number of regressors)") +
      ggplot2::ggtitle(paste0("Results with binomial-beta model prior (EMS = ", EMS, ")"))
  }

  if (dilution==1){
    Graph1_2 <- ggplot2::ggplot(forGraph1, ggplot2::aes(x = ID, y = Value)) +
      ggplot2::geom_line(ggplot2::aes(color = Probability, linetype = Probability)) +
      ggplot2::scale_color_manual(values = c("darkred", "steelblue")) +
      ggplot2::ylab("Prior, Posterior") + ggplot2::xlab("Model size (number of regressors)") +
      ggplot2::ggtitle(paste0("Results with diluted binomial model prior (EMS = ", EMS, ")"))

    Graph2_2 <- ggplot2::ggplot(forGraph2, ggplot2::aes(x = ID, y = Value)) +
      ggplot2::geom_line(ggplot2::aes(color = Probability, linetype = Probability)) +
      ggplot2::scale_color_manual(values = c("darkred", "steelblue")) +
      ggplot2::ylab("Prior, Posterior") + ggplot2::xlab("Model size (number of regressors)") +
      ggplot2::ggtitle(paste0("Results with diluted binomial-beta model prior (EMS = ", EMS, ")"))
  }

  # Putting together the last plot
  Finalplot <- ggpubr::ggarrange(Graph1_2,Graph2_2,
                                 labels = c("a)", "b)"),
                                 ncol = 1, nrow = 2, common.legend = TRUE, legend = "bottom")

  print(Finalplot)

  out <- list(Graph1, Graph2, Finalplot)
  return(out)
}
