#' Posterior probability of a positive coefficient sign (P(+))
#'
#' Computes posterior probabilities of a positive coefficient sign, P(+),
#' for the intercept and each regressor by averaging model-specific probabilities
#' across the model space, weighted by posterior model probabilities.
#'
#' For a given model \eqn{i} and coefficient \eqn{j}, the contribution is
#' \deqn{
#'   p(M_i \mid y) \cdot F_t\!\left(
#'     \frac{\beta_{ij}}{\sqrt{\mathrm{VAR}_{ij}}}; \mathrm{DF}_i
#'   \right),
#' }
#' where \eqn{F_t(\cdot;\mathrm{DF}_i)} is the CDF of the Student-\eqn{t}
#' distribution with \eqn{\mathrm{DF}_i} degrees of freedom.
#'
#' The intercept is included in all models, while each regressor contributes
#' only in those models in which it is included, as indicated by the model
#' inclusion matrix \code{Reg_ID}.
#'
#' @param pmp_uniform Numeric vector of length \code{MS} containing posterior
#'   model probabilities under a uniform model prior.
#' @param pmp_random Numeric vector of length \code{MS} containing posterior
#'   model probabilities under a random model prior.
#' @param betas Numeric matrix of dimension \code{MS x (K+1)} containing
#'   estimated coefficients for each model. Column 1 corresponds to the intercept,
#'   columns 2 to \code{K+1} correspond to regressors.
#' @param VAR Numeric matrix of dimension \code{MS x (K+1)} containing variances
#'   of the coefficient estimates. Must have the same dimensions as \code{betas}.
#' @param DF Numeric vector of length \code{MS} containing the degrees of freedom
#'   associated with each model.
#' @param Reg_ID Numeric or integer matrix of dimension \code{MS x K} indicating
#'   regressor inclusion. Entry \code{Reg_ID[i, k] = 1} if regressor \code{k}
#'   is included in model \code{i}, and 0 otherwise.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{Plus_PMP_uniform}{A \code{(K+1) x 1} numeric matrix containing
#'     posterior probabilities of a positive coefficient sign under the
#'     uniform model prior. The first row corresponds to the intercept.}
#'   \item{Plus_PMP_random}{A \code{(K+1) x 1} numeric matrix containing
#'     posterior probabilities of a positive coefficient sign under the
#'     random model prior. The first row corresponds to the intercept.}
#' }
#'
#' @details
#' The posterior probability of a positive sign for coefficient \eqn{j} is
#' computed as
#' \deqn{
#'   P(\beta_j > 0 \mid y)
#'   =
#'   \sum_{i \in \mathcal{M}_j}
#'   p(M_i \mid y)\,
#'   F_t\!\left(
#'     \frac{\beta_{ij}}{\sqrt{\mathrm{VAR}_{ij}}}; \mathrm{DF}_i
#'   \right),
#' }
#' where \eqn{\mathcal{M}_j} denotes the set of models that include regressor
#' \eqn{j}. For the intercept, \eqn{\mathcal{M}_j} contains all models.
#'
#' This definition follows the sign-probability interpretation in
#' \cite{Doppelhofer and Weeks (2009)}.
#'
#' @references
#' Doppelhofer, G. and Weeks, M. (2009).
#' \emph{Jointness of growth determinants}.
#' Journal of Applied Econometrics, 24(2), 209--244.
#'
#'
#' @export
plus_pmp_from_pmp <- function(pmp_uniform, pmp_random, betas, VAR, DF, Reg_ID) {
  MS <- nrow(Reg_ID)
  K  <- ncol(Reg_ID)

  stopifnot(
    nrow(betas) == MS,
    ncol(betas) == K + 1,
    all(dim(VAR) == dim(betas)),
    length(DF) == MS,
    length(pmp_uniform) == MS,
    length(pmp_random)  == MS
  )

  w_u <- as.numeric(pmp_uniform)
  w_r <- as.numeric(pmp_random)
  DF  <- as.numeric(DF)

  out_u <- numeric(K + 1)
  out_r <- numeric(K + 1)

  # ---- Intercept (always included) ----
  t0 <- betas[, 1] / sqrt(VAR[, 1])
  p0 <- stats::pt(t0, df = DF)

  out_u[1] <- sum(w_u * p0)
  out_r[1] <- sum(w_r * p0)

  # ---- Slopes: compute only where included ----
  B <- betas[, 2:(K + 1), drop = FALSE]  # MS x K
  V <- VAR[,   2:(K + 1), drop = FALSE]  # MS x K

  incl <- (Reg_ID == 1L)
  ok   <- incl & is.finite(B) & is.finite(V) & (V > 0)

  # Contribution matrix, zero by default (excluded regressors contribute 0)
  contrib_u <- matrix(0, nrow = MS, ncol = K)
  contrib_r <- matrix(0, nrow = MS, ncol = K)

  if (any(ok)) {
    t1 <- B[ok] / sqrt(V[ok])
    df_ok <- DF[row(B)[ok]]
    p1 <- stats::pt(t1, df = df_ok)

    # fill only included entries
    contrib_u[ok] <- w_u[row(B)[ok]] * p1
    contrib_r[ok] <- w_r[row(B)[ok]] * p1
  }

  out_u[2:(K + 1)] <- colSums(contrib_u)
  out_r[2:(K + 1)] <- colSums(contrib_r)

  list(
    Plus_PMP_uniform = matrix(out_u, ncol = 1),
    Plus_PMP_random  = matrix(out_r, ncol = 1)
  )
}
