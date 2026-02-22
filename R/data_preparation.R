#' Fixed-effects demeaning and data standardization
#'
#' Prepares a dataset for econometric analysis by applying fixed-effects
#' demeaning (within transformation) and/or standardization to numeric
#' variables. The behavior of the function depends on whether panel identifiers
#' are supplied and whether fixed effects are explicitly requested.
#'
#' If both \code{id} and \code{time} are provided and \code{fixed_effects = TRUE},
#' the function applies section, time, or two-way fixed-effects demeaning and may
#' optionally standardize the transformed variables. If \code{fixed_effects = FALSE},
#' fixed-effects demeaning is skipped even when identifiers are present, and only
#' standardization (if requested) is applied.
#'
#' If either \code{id} or \code{time} is missing, fixed-effects demeaning is not
#' available and the function requires \code{standardize = TRUE}.
#'
#' For two-way fixed effects, the transformation is:
#' \deqn{x_{it}^{*} = x_{it} - \bar{x}_{i\cdot} - \bar{x}_{\cdot t} + \bar{x}_{\cdot\cdot}}
#'
#' Standardization consists of subtracting the mean and dividing by the standard
#' deviation of each variable and is applied after fixed-effects demeaning
#' (if any).
#'
#' @param data A data.frame containing the data.
#' @param id An optional character string specifying the cross-sectional
#'   (section) identifier. Must be supplied together with \code{time} to
#'   enable fixed-effects demeaning.
#' @param time An optional character string specifying the time identifier.
#'   Must be supplied together with \code{id} to enable fixed-effects demeaning.
#' @param fixed_effects Logical. If \code{TRUE}, fixed-effects demeaning is
#'   applied when both \code{id} and \code{time} are provided. If \code{FALSE},
#'   fixed-effects demeaning is skipped even when identifiers are present.
#' @param effect A character string indicating the fixed-effects structure
#'   when \code{fixed_effects = TRUE}. One of \code{"twoway"},
#'   \code{"section"}, or \code{"time"}.
#' @param standardize Logical. If \code{TRUE}, numeric variables are
#'   standardized by subtracting their mean and dividing by their standard
#'   deviation. When fixed effects are applied, standardization occurs
#'   after demeaning.
#'
#' @return
#' A data.frame containing only numeric variables used in estimation.
#' Panel identifiers (\code{id}, \code{time}) are removed from the output.
#' Transformed variables preserve their original column names.
#'
#' @details
#' The function operates in three modes:
#' \itemize{
#'   \item \strong{Fixed effects only}: \code{fixed_effects = TRUE},
#'   \code{standardize = FALSE}.
#'   \item \strong{Fixed effects + standardization}: \code{fixed_effects = TRUE},
#'   \code{standardize = TRUE}.
#'   \item \strong{Standardization only}: \code{fixed_effects = FALSE},
#'   \code{standardize = TRUE}.
#' }
#'
#' When \code{id} and \code{time} are not provided, only the standardization-only
#' mode is available.
#'
#' Missing values are ignored when computing means and standard deviations.
#' After fixed-effects demeaning, an intercept term is redundant in subsequent
#' linear regressions.
#'
#' @examples
#' \dontrun{
#' # Standardization only (panel identifiers present but FE skipped)
#' X <- data_preparation(
#'   df,
#'   id = "Section",
#'   time = "Time",
#'   fixed_effects = FALSE,
#'   standardize = TRUE
#' )
#'
#' # Two-way fixed effects with standardization
#' X <- data_preparation(
#'   df,
#'   id = "Section",
#'   time = "Time",
#'   fixed_effects = TRUE,
#'   effect = "twoway",
#'   standardize = TRUE
#' )
#'
#' # Section fixed effects only
#' X <- data_preparation(
#'   df,
#'   id = "Section",
#'   time = "Time",
#'   fixed_effects = TRUE,
#'   effect = "section"
#' )
#'
#' # Standardization only (no panel identifiers)
#' X <- data_preparation(df, standardize = TRUE)
#' }
#'
#' @export

data_preparation <- function(data,
                             id = NULL,
                             time = NULL,
                             fixed_effects = FALSE,
                             effect = c("twoway", "section", "time"),
                             standardize = FALSE) {

  stopifnot(is.data.frame(data))
  out <- data

  num_vars <- names(out)[sapply(out, is.numeric)]
  if (length(num_vars) == 0)
    stop("No numeric variables found in `data`.")

  has_ids <- !is.null(id) && !is.null(time)

  # If IDs are provided, we standardize/demean only numeric vars excluding id/time
  if (has_ids) {
    stopifnot(id %in% names(out), time %in% names(out))
    vars <- setdiff(num_vars, c(id, time))
    if (length(vars) == 0)
      stop("No numeric variables to transform (excluding `id` and `time`).")

    # Apply FE demeaning only if requested
    if (fixed_effects) {
      effect <- match.arg(effect)

      gmean <- function(x, g) stats::ave(x, g, FUN = function(z) mean(z, na.rm = TRUE))

      for (v in vars) {
        x <- out[[v]]

        mi <- if (effect %in% c("twoway", "section")) gmean(x, out[[id]])   else 0
        mt <- if (effect %in% c("twoway", "time"))    gmean(x, out[[time]]) else 0
        mg <- if (effect == "twoway") mean(x, na.rm = TRUE) else 0

        out[[v]] <- x - mi - mt + mg
      }
    }

    # Standardization (after demeaning if fixed_effects==TRUE)
    if (standardize) {
      for (v in vars) {
        x <- out[[v]]
        s <- stats::sd(x, na.rm = TRUE)
        if (is.finite(s) && s > 0) {
          out[[v]] <- (x - mean(x, na.rm = TRUE)) / s
        } else {
          out[[v]] <- NA_real_
        }
      }
    } else if (!fixed_effects) {
      # If neither FE nor standardization requested, do nothing but drop ids
      # (kept silent; remove this branch if you prefer an error)
    }

    # Drop id/time from final output
    out[[id]] <- NULL
    out[[time]] <- NULL

    return(out)
  }

  # No IDs provided: FE not available; only standardization allowed
  if (!standardize) {
    stop("When `id` and `time` are not provided, set `standardize = TRUE` (only standardization is available).")
  }

  for (v in num_vars) {
    x <- out[[v]]
    s <- stats::sd(x, na.rm = TRUE)
    if (is.finite(s) && s > 0) {
      out[[v]] <- (x - mean(x, na.rm = TRUE)) / s
    } else {
      out[[v]] <- NA_real_
    }
  }

  out
}
