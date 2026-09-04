#' Graphs of the prior and posterior model probabilities for the best individual models
#'
#' This function draws four graphs of prior and posterior model probabilities for the best individual models: \cr
#' a) The results with binomial model prior (based on PMP - posterior model probability) \cr
#' b) The results with binomial-beta model prior (based on PMP - posterior model probability) \cr
#' Models on the graph are ordered according to their posterior model probability.
#'
#'
#' @param bma_list bma_list object (the result of the bma function)
#' @param top The number of the best model to be placed on the graphs
#' @param type Character, either \code{"line"} (the default) for the prior and
#' posterior drawn as lines against the model ranking, or \code{"histogram"}
#' for them drawn as side-by-side bars. The bar form is best suited to a small
#' \code{top}, where individual models can still be told apart; with a large
#' \code{top} the line form remains readable. The two are otherwise identical.
#'
#' @return A list with three graphs with prior and posterior model probabilities for individual models:\cr
#' 1) The results with binomial model prior (based on PMP - posterior model probability) \cr
#' 2) The results with binomial-beta model prior (based on PMP - posterior model probability) \cr
#' 3) On graph combining the aforementioned graphs
#'
#' @export
#'
#' @examples
#' \donttest{
#' data("Trade_data", package = "rmsBMA")
#' data <- Trade_data[,1:10]
#' modelSpace <- model_space(data, M = 9, g = "UIP")
#' bma_list <- bma(modelSpace)
#' model_pmps <- model_pmp(bma_list, top = 100)
#' model_pmps[[1]]
#' model_bars <- model_pmp(bma_list, top = 10, type = "histogram")
#' model_bars[[1]]
#'}
#'@name model_pmp

utils::globalVariables(c("ID", "Value", "Probability"))

model_pmp <- function(bma_list, top = NULL, type = c("line", "histogram")){

type <- match.arg(type)

# Collecting information from the bma_list
R <- bma_list[[6]] # total number of regressors
M <- bma_list[[7]] # size of the model space
EMS <- bma_list[[8]] # expected model size


PMPs <- bma_list[[12]][,3:4] # PMP_uniform, PMP_random
Priors <- bma_list[[12]][,1:2] # Priors: uniform and random
dilution <- bma_list[[9]] # 0 - no dilution prior, 1 - dilution prior

if (is.null(top)){
  top <- M
}

if (top>M){# CONDITION about what to do if the user sets top that is higher than M
  message("The number of the best models (top) cannot be higher than the total ",
          "number of models. Setting top = ", M, ", the size of the model space.")
  top = M
}

# Objects to store posteriors and priors
PMP_uniform <- cbind(PMPs[,1], Priors[,1])
PMP_random <- cbind(PMPs[,2], Priors[,2])

# Ordering of the models according to posterior criterion
PMP_uniform <- PMP_uniform[order(PMP_uniform[,1], decreasing=T),]
PMP_random <- PMP_random[order(PMP_random[,1], decreasing=T),]

ranking <- matrix(1:M, nrow = M, ncol = 1)

# Adding a ranking number
PMP_uniform <- cbind(ranking[1:top,], PMP_uniform[1:top,])
PMP_random <- cbind(ranking[1:top,], PMP_random[1:top,])

IDnames <- cbind("ID","Posterior","Prior") # names of the variables to be used by 'tidyverse'

colnames(PMP_uniform) <- IDnames
colnames(PMP_random) <- IDnames

forGraph1 <- as.data.frame(PMP_uniform)
forGraph2 <- as.data.frame(PMP_random)

## Preparation of the Figures with ggplot
forGraph1 <- tidyr::gather(forGraph1, key = "Probability", value = "Value", -ID)

forGraph2 <- tidyr::gather(forGraph2, key = "Probability", value = "Value", -ID)

# One builder for both display types, both priors and both dilution settings,
# so that they cannot drift apart. The line form is unchanged apart from the
# axis label, which previously read "raniking" on the two untitled graphs.
pmp_plot <- function(df, title = NULL) {
  p <- ggplot2::ggplot(df, ggplot2::aes(x = ID, y = Value))
  if (type == "line") {
    p <- p +
      ggplot2::geom_line(ggplot2::aes(color = Probability, linetype = Probability)) +
      ggplot2::scale_color_manual(values = c("darkred", "steelblue"))
  } else {
    p <- p +
      ggplot2::geom_col(ggplot2::aes(fill = Probability),
                        position = ggplot2::position_dodge(width = 0.75),
                        width = 0.7) +
      ggplot2::scale_fill_manual(values = c("darkred", "steelblue"))
    # label every rank while there are few enough for that to stay legible
    if (top <= 30) p <- p + ggplot2::scale_x_continuous(breaks = 1:top)
  }
  p <- p +
    ggplot2::ylab("Prior, Posterior") +
    ggplot2::xlab("Model number in the ranking")
  if (!is.null(title)) p <- p + ggplot2::ggtitle(title)
  p
}

Graph1 <- pmp_plot(forGraph1)
Graph2 <- pmp_plot(forGraph2)

## Titled versions for the BIG COMBINED GRAPH
dil_label <- if (identical(as.numeric(dilution), 1)) "diluted " else ""
Graph1_2 <- pmp_plot(forGraph1, paste0("Results with ", dil_label,
                                       "binomial model prior (EMS = ", EMS, ")"))
Graph2_2 <- pmp_plot(forGraph2, paste0("Results with ", dil_label,
                                       "binomial-beta model prior (EMS = ", EMS, ")"))

# Putting together the last plot
Finalplot <- ggpubr::ggarrange(Graph1_2,Graph2_2,
                               labels = c("a)", "b)"),
                               ncol = 1, nrow = 2, common.legend = TRUE, legend = "bottom")

print(Finalplot)

out <- list(Graph1, Graph2, Finalplot)
return(out)
}
