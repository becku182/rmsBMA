#' Extreme bounds analysis summaries from a Bayesian model space
#'
#' Computes Extreme Bounds Analysis (EBA) summaries for the intercept and each
#' regressor across a model space. For each coefficient, the function reports:
#' the minimum coefficient ("Low"), maximum coefficient ("High"), the mean
#' coefficient ("Mean_coef"), and corresponding "extreme bounds" defined as
#' \eqn{\mathrm{Low} - 2\cdot \mathrm{SE}} and \eqn{\mathrm{High} + 2\cdot \mathrm{SE}},
#' where \eqn{\mathrm{SE}=\sqrt{\mathrm{VAR}}} is the standard error associated
#' with the coefficient estimate in the model attaining the minimum/maximum.
#'
#' The intercept (constant) is assumed to be included in all models. Each
#' regressor is summarized only over models in which it is included, as indicated
#' by the model-inclusion matrix \code{Reg_ID}.
#'
#' @param betas Numeric matrix of dimension \code{MS x (K+1)} containing
#'   estimated coefficients across models. Column 1 corresponds to the intercept,
#'   columns 2 to \code{K+1} correspond to regressors.
#' @param VAR Numeric matrix of dimension \code{MS x (K+1)} containing variances
#'   of the coefficient estimates. Must have the same dimensions as \code{betas}.
#' @param Reg_ID Numeric or integer matrix of dimension \code{MS x K} indicating
#'   regressor inclusion. Entry \code{Reg_ID[i,k]=1} if regressor \code{k} is
#'   included in model \code{i}, and 0 otherwise.
#' @param var_tol Nonnegative numeric scalar used as a tolerance when checking
#'   variance positivity. Entries with \code{VAR <= var_tol} are treated as invalid
#'   for bound calculations. Default is \code{0}.
#'
#' @return A numeric matrix of dimension \code{(K+1) x 5} with columns:
#' \describe{
#'   \item{Lower_bound}{\eqn{\min(\beta) - 2\cdot \mathrm{SE}} evaluated at the model
#'   where \eqn{\beta} is minimal.}
#'   \item{Low}{Minimum coefficient value across relevant models.}
#'   \item{Mean_coef}{Mean coefficient across relevant models (intercept: all models;
#'   regressor: included models only).}
#'   \item{High}{Maximum coefficient value across relevant models.}
#'   \item{Upper_bound}{\eqn{\max(\beta) + 2\cdot \mathrm{SE}} evaluated at the model
#'   where \eqn{\beta} is maximal.}
#' }
#' Rows correspond to the intercept (row 1) and regressors (rows 2..K+1).
#' If a regressor is never included (no 1s in its column of \code{Reg_ID}), its row
#' will contain \code{NA}.
#'
#' @examples
#' \dontrun{
#' eba_tbl <- eba_bounds(betas = betas, VAR = VAR, Reg_ID = Reg_ID)
#' rownames(eba_tbl) <- c("Const", x_names)
#' eba_tbl
#' }
#'
#' @export
eba <- function(betas, VAR, Reg_ID, var_tol = 0) {
  betas <- as.matrix(betas)
  VAR   <- as.matrix(VAR)
  Reg_ID <- as.matrix(Reg_ID)

  if (!is.numeric(var_tol) || length(var_tol) != 1L || is.na(var_tol) || var_tol < 0) {
    stop("var_tol must be a single nonnegative number.")
  }

  MS <- nrow(Reg_ID)
  K  <- ncol(Reg_ID)

  if (nrow(betas) != MS) stop("nrow(betas) must equal nrow(Reg_ID) = MS.")
  if (ncol(betas) != K + 1L) stop("ncol(betas) must equal K+1 (intercept + K regressors).")
  if (!all(dim(VAR) == dim(betas))) stop("VAR must have the same dimensions as betas.")
  if (any(VAR < 0, na.rm = TRUE)) stop("VAR contains negative values.")

  out <- matrix(NA_real_, nrow = K + 1L, ncol = 5)
  colnames(out) <- c("Lower_bound", "Low", "Mean_coef", "High", "Upper_bound")

  compute_stats <- function(b, se) {
    lo_idx <- which.min(b)
    hi_idx <- which.max(b)
    c(
      Lower_bound = b[lo_idx] - 2 * se[lo_idx],
      Low         = b[lo_idx],
      Mean_coef   = mean(b),
      High        = b[hi_idx],
      Upper_bound = b[hi_idx] + 2 * se[hi_idx]
    )
  }

  # Intercept (always included)
  b0  <- betas[, 1]
  v0  <- VAR[, 1]
  se0 <- sqrt(v0)

  keep0 <- is.finite(b0) & is.finite(se0) & (v0 > var_tol)
  if (!any(keep0)) stop("No valid intercept entries found (check VAR[,1]).")
  out[1, ] <- compute_stats(b0[keep0], se0[keep0])

  # Regressors (only where included)
  for (k in 1:K) {
    rows <- which(Reg_ID[, k] == 1L)
    if (length(rows) == 0L) next

    b  <- betas[rows, k + 1L]
    v  <- VAR[rows,   k + 1L]
    se <- sqrt(v)

    keep <- is.finite(b) & is.finite(se) & (v > var_tol)
    if (!any(keep)) next

    out[k + 1L, ] <- compute_stats(b[keep], se[keep])
  }

  out
}
