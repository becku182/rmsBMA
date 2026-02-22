#' Table with the best models according to one of the posterior criteria
#'
#' This function creates a ranking of best models according to one of the possible criterion (PMP under binomial model prior, PMP under binomial-beta model prior, R^2 under binomial model prior, R^2 under binomial-beta model prior).
#' The function gives two types of tables in three different formats: inclusion table (where 1 indicates presence of the regressor in the model and 0 indicates that the variable is excluded from the model) and estimation results table (it displays the best models and estimation output for those models: point estimates, standard errors, significance level, and R^2).
#'
#' @param bma_list bma object (the result of the bma function)
#' @param criterion The criterion that will be used for a basis of the model ranking: \cr
#' 1 - binomial model prior \cr
#' 2 - binomial-beta model prior
#' @param best The number of the best models to be considered
#' @param round Parameter indicating the decimal place to which number in the tables should be rounded (default round = 3)
#' @param estimate A parameter with values TRUE or FALSE indicating which table should be displayed when
#' TRUE - table with estimation to the results \cr
#' FALSE - table with the inclusion of regressors in the best models
#'
#' @return A list with best_models objects: \cr
#' 1. matrix with inclusion of the regressors in the best models \cr
#' 2. matrix with estimation output in the best models with regular standard errors \cr
#' 3. knitr_kable table with inclusion of the regressors in the best models (the best for the display on the console - up to 11 models) \cr
#' 4. knitr_kable table with estimation output in the best models with regular standard errors (the best for the display on the console - up to 6 models) \cr
#' 5. gTree table with inclusion of the regressors in the best models (displayed as a plot). Use grid::grid.draw() to display.\cr
#' 6. gTree table with estimation output in the best models with regular standard errors (displayed as a plot). Use grid::grid.draw() to display.
#'
#' @export
#'
#' @examples
#' \donttest{
#' data <- Trade_data[,1:10]
#' modelSpace <- model_space(data, M = 9, g = "UIP")
#' bma_list <- bma(modelSpace)
#' models <- best_models(bma_list, best = 3)
#' models[[4]]
#' }
#'@name best_models

best_models <- function(bma_list, criterion = 1, best = 5, round = 3, estimate = TRUE){

  R <- bma_list[[6]] # number of regressors from bma object
  K <- R+1 # number of variables
  reg_names <- matrix(c("Const",bma_list[[5]]), nrow = K, ncol = 1) # vector with names of the regressors from bma object
  M <- bma_list[[7]] # size of the mode space from bma object
  info <- bma_list[[11]][,1:(R+2*K)]
  PMP_unifrom <- matrix(bma_list[[11]][,R+2*K+3], nrow = M, ncol = 1)
  PMP_random <- matrix(bma_list[[11]][,R+2*K+4], nrow = M, ncol = 1)
  d_free <- matrix(bma_list[[11]][,R+2*K+1], nrow = M, ncol = 1)
  R2 <- matrix(bma_list[[11]][,R+2*K+2], nrow = M, ncol = 1)

  if (best>M){
    message("best > M - number of best models cannot be bigger than the total number of models. We set best = M and continiue :)")
    best = M
  }

  # check for the criterion chosen by the user
  if (criterion==1){ranking <- PMP_unifrom}
  if (criterion==2){ranking <- PMP_random}

  Ranking<-cbind(ranking,info,d_free,R2) # we add ranking criterion based on the users choice

  # ordering the models according to PMP criterion
  Ranking <- Ranking[order(Ranking[,1],decreasing=T),] # ordering of the models

  Best_models <- Ranking[1:best, 2:(R+1)] # model IDs
  Ranks <- matrix(round(Ranking[1:best, 1], digits = 3), nrow = best, ncol = 1) # PMPs of the first 'best' models
  bestBetas <- Ranking[1:best, (R+2):(R+K+1)] # coefficients
  bestBetas <- t(round(bestBetas, round))
  bestBetas[bestBetas == 0] <- NA
  bestSTDs <- Ranking[1:best, (R+K+2):(R+2*K+1)] # standard errors
  bestSTDs <- t(round(bestSTDs, round))
  bestBetas[bestBetas == 0] <- NA

  best_d_free <- matrix( Ranking[1:best, R+2*K+2], nrow = best, ncol = 1)
  best_R2 <- matrix( Ranking[1:best, R+2*K+3], nrow = best, ncol = 1)
  best_R2 <- t(round(best_R2, round))

  inclusion_table <- t(cbind(matrix(1, nrow = best, ncol = 1), Best_models, Ranks, t(best_R2)))
  row.names(inclusion_table) <- rbind(reg_names, "PMP", "R^2")

  names <- matrix(0, nrow = best, ncol = 1)

  for (i in 1:best){
    names[i,1] = paste0("'No. ",i,"'")
  }

  colnames(inclusion_table) <- names

  models_std <- matrix(0, nrow = K, ncol = best)
  p_values <- matrix(0, nrow = K, ncol = best)
  asterisks <- matrix(0, nrow = K, ncol = best)

  for (i in 1:K){
    for (j in 1:best){
      if (!is.na(bestBetas[i,j])){
        models_std[i,j] = paste0(bestBetas[i,j]," (",bestSTDs[i,j],")")
        p_values[i,j] = 2*stats::pt(abs(bestBetas[i,j]/bestSTDs[i,j]), df = best_d_free[j,1], lower.tail = FALSE)

        if (is.na(p_values[i,j]) || p_values[i,j] >= 0.1){
          asterisks[i,j] = NA
        } else if (p_values[i,j] >= 0.05){
          asterisks[i,j] = "*"
        } else if (p_values[i,j] >= 0.01){
          asterisks[i,j] = "**"
        } else {
          asterisks[i,j]="***"
        }

      } else{
        models_std[i,j] = NA
        p_values[i,j] = NA
        asterisks[i,j] = NA
      }
    }
  }

  for (i in 1:K){
    for (j in 1:best){
      if (!is.na(asterisks[i,j])){
        models_std[i,j] = paste0(models_std[i,j],asterisks[i,j])
      }
    }
  }

  models_std <- rbind(models_std, t(Ranks), best_R2)

  colnames(models_std) <- names
  row.names(models_std) <- rbind(reg_names,"PMP","R^2")

  inclusion_2 <- knitr::kable(inclusion_table, row.names = TRUE, align = "c")
  models_std_2 <- knitr::kable(models_std, row.names = TRUE, align = "c")
  inclusion_3 <- grid::grid.grabExpr(gridExtra::grid.table(inclusion_table))
  models_std_3 <- grid::grid.grabExpr(gridExtra::grid.table(models_std))

  out <- list(inclusion_table, models_std, inclusion_2,
              models_std_2, inclusion_3, models_std_3)

  if (estimate==FALSE){
    gridExtra::grid.table(inclusion_table)
  }
  if (estimate==TRUE){
      gridExtra::grid.table(models_std)
  }
  return(out)
}
