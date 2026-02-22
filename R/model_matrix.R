#' Matrix reflecting model space
#'
#' The function creates a matrix with ones indicating inclusion and zeros indicating inclusion
#' of a regressor in a model
#'
#' @param K total number of regressors
#' @param M maximum number of regressor in a model
#'
#' @return A matrix with ones indicating inclusion and zeros indicating inclusion
#' of a regressor in a model
#'

model_matrix <- function(K, M) {
  stopifnot(M <= K)

  models <- vector("list", sum(choose(K, 0:M)))
  idx <- 1

  # k = 0 (null model)
  models[[idx]] <- integer(K)
  idx <- idx + 1

  # k = 1,...,M
  for (k in 1:M) {
    cmb <- utils::combn(K, k)
    for (j in seq_len(ncol(cmb))) {
      v <- integer(K)
      v[cmb[, j]] <- 1
      models[[idx]] <- v
      idx <- idx + 1
    }
  }

  do.call(rbind, models)
}
