#' Posterior sign certainty probability
#'
#' Computes the posterior probability that a regression coefficient has the
#' same sign as its posterior mean, following the sign-certainty measure used
#' in Sala-i-Martin, Doppelhofer and Miller (2004) and Doppelhofer and Weeks (2009).
#'
#' For each coefficient \eqn{j}, define
#' \deqn{
#' S_j = \sum_{i=1}^{MS} p(M_i \mid y)\,
#' F_t\!\left(
#'   \frac{\beta_{ij}}{\sqrt{\mathrm{VAR}_{ij}}}; \mathrm{DF}_i
#' \right),
#' }
#' where \eqn{F_t(\cdot;\mathrm{DF}_i)} is the CDF of the Student-\eqn{t}
#' distribution with \eqn{\mathrm{DF}_i} degrees of freedom.
#'
#' The posterior sign certainty probability is then defined as
#' \deqn{
#' p(\mathrm{sign}_j \mid y) =
#' \begin{cases}
#'   S_j,       & \text{if } \mathrm{sign}(E[\beta_j \mid y]) > 0, \\
#'   1 - S_j,   & \text{if } \mathrm{sign}(E[\beta_j \mid y]) < 0, \\
#'   0.5,       & \text{if } E[\beta_j \mid y] = 0.
#' \end{cases}
#' }
#'
#' The intercept is included in all models. For slope coefficients that are
#' excluded from a given model, the contribution to \eqn{S_j} is set to
#' \eqn{F_t(0)=0.5}, reflecting a symmetric distribution centered at zero.
#'
#' @param pmp_uniform Numeric vector of length \code{MS} containing posterior
#'   model probabilities under a uniform model prior.
#' @param pmp_random Numeric vector of length \code{MS} containing posterior
#'   model probabilities under a random model prior.
#' @param betas Numeric matrix of dimension \code{MS x (K+1)} containing
#'   estimated coefficients for each model. Column 1 corresponds to the intercept.
#' @param VAR Numeric matrix of dimension \code{MS x (K+1)} containing variances
#'   of the coefficient estimates. Must have the same dimensions as \code{betas}.
#' @param DF Numeric vector of length \code{MS} giving the degrees of freedom
#'   associated with each model.
#' @param Reg_ID Numeric or integer matrix of dimension \code{MS x K} indicating
#'   regressor inclusion. Entry \code{Reg_ID[i, k] = 1} if regressor \code{k}
#'   is included in model \code{i}, and 0 otherwise.
#'
#' @return A list with four elements:
#' \describe{
#'   \item{PSC_uniform}{A \code{(K+1) x 1} numeric matrix containing posterior
#'     sign certainty probabilities under the uniform model prior.}
#'   \item{PSC_random}{A \code{(K+1) x 1} numeric matrix containing posterior
#'     sign certainty probabilities under the random model prior.}
#'   \item{PostMean_uniform}{A \code{(K+1) x 1} numeric matrix of posterior means
#'     under the uniform model prior.}
#'   \item{PostMean_random}{A \code{(K+1) x 1} numeric matrix of posterior means
#'     under the random model prior.}
#' }
#'
#' @references
#' Sala-i-Martin, X., Doppelhofer, G., and Miller, R. I. (2004).
#' Determinants of long-term growth: A Bayesian averaging of classical estimates.
#' \emph{American Economic Review}, 94(4), 813--835.
#'
#' Doppelhofer, G. and Weeks, M. (2009).
#' Jointness of growth determinants.
#' \emph{Journal of Applied Econometrics}, 24(2), 209--244.
#'
#' @export
posterior_sign_certainty <- function(
    pmp_uniform,
    pmp_random,
    betas,
    VAR,
    DF,
    Reg_ID
) {

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

  DF  <- as.numeric(DF)
  w_u <- as.numeric(pmp_uniform)
  w_r <- as.numeric(pmp_random)

  # ---------- Posterior means ----------
  mu_u <- numeric(K + 1)
  mu_r <- numeric(K + 1)

  mu_u[1] <- sum(w_u * betas[, 1])
  mu_r[1] <- sum(w_r * betas[, 1])

  B <- betas[, 2:(K + 1), drop = FALSE]
  mu_u[2:(K + 1)] <- colSums((B * Reg_ID) * w_u)
  mu_r[2:(K + 1)] <- colSums((B * Reg_ID) * w_r)

  # ---------- CDF contributions ----------
  # slopes: default 0.5 for excluded regressors
  C <- matrix(0.5, nrow = MS, ncol = K)

  V <- VAR[, 2:(K + 1), drop = FALSE]
  incl <- (Reg_ID == 1L)
  ok   <- incl & is.finite(B) & is.finite(V) & (V > 0)

  if (any(ok)) {
    t1 <- B[ok] / sqrt(V[ok])
    df_ok <- DF[row(B)[ok]]
    C[ok] <- stats::pt(t1, df = df_ok)
  }

  # intercept
  if (any(VAR[, 1] <= 0)) {
    stop("Intercept variance VAR[,1] must be strictly positive.")
  }
  t0 <- betas[, 1] / sqrt(VAR[, 1])
  c0 <- stats::pt(t0, df = DF)

  S_u <- numeric(K + 1)
  S_r <- numeric(K + 1)

  S_u[1] <- sum(w_u * c0)
  S_r[1] <- sum(w_r * c0)

  S_u[2:(K + 1)] <- colSums(C * w_u)
  S_r[2:(K + 1)] <- colSums(C * w_r)

  # ---------- Sign certainty ----------
  PSC_uniform <- ifelse(mu_u > 0, S_u,
                        ifelse(mu_u < 0, 1 - S_u, 0.5))
  PSC_random  <- ifelse(mu_r > 0, S_r,
                        ifelse(mu_r < 0, 1 - S_r, 0.5))

  list(
    PSC_uniform = matrix(PSC_uniform, ncol = 1),
    PSC_random  = matrix(PSC_random,  ncol = 1),
    PostMean_uniform = matrix(mu_u, ncol = 1),
    PostMean_random  = matrix(mu_r, ncol = 1)
  )
}
