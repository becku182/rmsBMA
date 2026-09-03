#' Regressors choice based on the inclusion vector
#'
#' The function creates a matrix with regressors that should be included in a specific model
#' based on inclusion vector
#'
#' @param x matrix with regressors
#' @param model_row inclusion vector - row of a model space matrix
#'
#' @return A matrix with regressors to be used in a specific model
#'
#'
#' @keywords internal
subset_design <- function(x, model_row) {
  idx <- which(model_row != 0L)   # or == 1L
  x[, idx, drop = FALSE]
}
