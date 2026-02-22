#' Calculation of of the jointness measures
#'
#' This function calculates four types of the jointness measures based on the posterior model probabilities calculated using binomial and binomial-beta model prior. The four measures are: \cr
#' 1) HCGHM - for Hofmarcher et al. (2018) measure; \cr
#' 2) LS - for Ley & Steel (2007) measure; \cr
#' 3) DW - for Doppelhofer & Weeks (2009) measure; \cr
#' 4) PPI - for posterior probability of including both variables. \cr
#' The measures under binomial model prior will appear in a table above the diagonal, and the measure calculated under binomial-beta model prior below the diagonal. \cr
#' \cr
#' REFERENCES \cr
#' Doppelhofer G, Weeks M (2009) Jointness of growth determinants. Journal of Applied Econometrics., 24(2), 209-244. doi: 10.1002/jae.1046 \cr
#' Hofmarcher P, Crespo Cuaresma J, Grün B, Humer S, Moser M (2018) Bivariate jointness measures in Bayesian Model Averaging: Solving the conundrum. Journal of Macroeconomics, 57, 150-165. doi: 10.1016/j.jmacro.2018.05.005 \cr
#' Ley E, Steel M (2007) Jointness in Bayesian variable selection with applications to growth regression. Journal of Macroeconomics, 29(3), 476-493. doi: 10.1016/j.jmacro.2006.12.002
#'
#' @param bma_list bma object (the result of the bma function)
#' @param measure Parameter for choosing the measure of jointness: \cr
#' HCGHM - for Hofmarcher et al. (2018) measure; \cr
#' LS - for Ley & Steel (2007) measure; \cr
#' DW - for Doppelhofer & Weeks (2009) measure; \cr
#' PPI - for posterior probability of including both variables.
#' @param rho The parameter "rho" (\eqn{\rho}) to be used in HCGHM jointness measure (default rho = 0.5). Works only if HCGHM measure is chosen (Hofmarcher et al. 2018).
#' @param round Parameter indicating the decimal place to which the jointness measures should be rounded (default round = 3).
#'
#' @return A table with jointness measures for all the pairs of regressors used in the analysis. The results obtained with the binomial model prior are above the diagonal, while the ones obtained with the binomial-beta prior are below.
#'
#' @export
#'
#' @examples
#'
#' \donttest{
#' data("Trade_data", package = "rmsBMA")
#' data <- Trade_data[,1:10]
#' modelSpace <- model_space(data, M = 9, g = "UIP")
#' bma_list <- bma(modelSpace)
#' jointness_table <- jointness(bma_list)
#' }

jointness <- function(bma_list, measure = "HCGHM", rho = 0.5, round = 3){

  reg_names <- bma_list[[5]]
  R <- bma_list[[6]]
  M <- bma_list[[7]]
  forJointness <- bma_list[[10]]

  # Inclusion matrix (0/1): M x R
  Z <- as.matrix(forJointness[, 1:R])
  # weights: length M
  wU <- as.numeric(forJointness[, R+1])  # uniform prior PMP
  wR <- as.numeric(forJointness[, R+2])  # random prior PMP

  # Weighted joint inclusion matrices: R x R
  Pab_U <- crossprod(Z * wU, Z)
  Pab_R <- crossprod(Z * wR, Z)

  # Marginal inclusion probs: length R
  Pa_U <- as.numeric(crossprod(wU, Z))
  Pa_R <- as.numeric(crossprod(wR, Z))

  # Build component matrices (R x R)
  # (diag elements are fine; we'll set them to NA at end)
  Na_b_U  <- matrix(Pa_U, nrow = R, ncol = R, byrow = TRUE) - Pab_U
  a_Nb_U  <- matrix(Pa_U, nrow = R, ncol = R, byrow = FALSE) - Pab_U
  Na_Nb_U <- 1 - matrix(Pa_U, nrow = R, ncol = R, byrow = FALSE) -
    matrix(Pa_U, nrow = R, ncol = R, byrow = TRUE) + Pab_U

  Na_b_R  <- matrix(Pa_R, nrow = R, ncol = R, byrow = TRUE) - Pab_R
  a_Nb_R  <- matrix(Pa_R, nrow = R, ncol = R, byrow = FALSE) - Pab_R
  Na_Nb_R <- 1 - matrix(Pa_R, nrow = R, ncol = R, byrow = FALSE) -
    matrix(Pa_R, nrow = R, ncol = R, byrow = TRUE) + Pab_R

  # Compute the measure matrices
  if (measure == "HCGHM") {
    num_U <- (Pab_U + rho) * (Na_Nb_U + rho) - (Na_b_U + rho) * (a_Nb_U + rho)
    den_U <- (Pab_U + rho) * (Na_Nb_U + rho) + (Na_b_U + rho) * (a_Nb_U + rho) - rho
    first  <- num_U / den_U

    num_R <- (Pab_R + rho) * (Na_Nb_R + rho) - (Na_b_R + rho) * (a_Nb_R + rho)
    den_R <- (Pab_R + rho) * (Na_Nb_R + rho) + (Na_b_R + rho) * (a_Nb_R + rho) - rho
    second <- num_R / den_R

  } else if (measure == "LS") {
    first  <- Pab_U / (Na_b_U + a_Nb_U)
    second <- Pab_R / (Na_b_R + a_Nb_R)

  } else if (measure == "DW") {
    first  <- log((Pab_U / Na_b_U) * (Na_Nb_U / a_Nb_U))
    second <- log((Pab_R / Na_b_R) * (Na_Nb_R / a_Nb_R))

  } else if (measure == "PPI") {
    first  <- Pab_U
    second <- Pab_R

  } else {
    stop("Unknown measure. Use one of: 'HCGHM', 'LS', 'DW', 'PPI'.")
  }

  # Assemble table: above diag = uniform, below diag = random
  out <- first
  out[lower.tri(out)] <- second[lower.tri(second)]

  diag(out) <- NA_real_
  out <- round(out, digits = round)

  rownames(out) <- reg_names
  colnames(out) <- reg_names
  out
}
