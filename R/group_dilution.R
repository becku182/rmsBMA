#' Compute group-based dilution factors for a model space
#'
#' Computes a group-based dilution factor for each model (row) in a model-inclusion
#' matrix. Regressors are assigned to "relatedness groups" via `Nar_vec`. For each
#' model and each group, the dilution exponent increases by 1 for every additional
#' regressor from that group beyond the first. The model's dilution factor is the
#' product of group-specific penalties.
#'
#' Formally, for model i and group h >= 1, let
#' \deqn{c_{ih} = \sum_{j=1}^K \gamma_{ij} I\{g(j)=h\}}
#' and
#' \deqn{D_i = \sum_{h \ge 1} \max(0, c_{ih} - 1).}
#' If `p` is a scalar, the dilution factor is \deqn{p^{D_i}}.
#' If `p` is group-specific with values p_h, then
#' \deqn{\mathrm{NarDilution}_i = \prod_{h \ge 1} p_h^{\max(0, c_{ih}-1)}.}
#'
#' @param Reg_ID An `MS x K` numeric/integer matrix of model indicators. Each row
#'   corresponds to a model; each column corresponds to a regressor. Entries should
#'   be 0/1 (0 = excluded, 1 = included).
#' @param Nar_vec An integer vector of length `K` giving group assignments for each
#'   regressor. Use 0 for regressors not belonging to any group. Positive integers
#'   (1,2,...) denote group IDs.
#' @param p Either:
#'   \itemize{
#'     \item a single numeric scalar in \eqn{[0,1]} applied to all groups, or
#'     \item a numeric vector in \eqn{[0,1]} with one entry per group.
#'   }
#'   If `p` is a vector, it is matched to groups as follows:
#'   \itemize{
#'     \item If `p` has names, they must match group IDs (e.g., `c("1"=0.7,"2"=0.5)`),
#'     \item otherwise it is assumed to be in the order of `sort(unique(Nar_vec[Nar_vec>0]))`.
#'   }
#'
#' @return A numeric vector of length `MS` containing the dilution factor for each model.
#'
#' @examples
#' # Example model space: MS models, K regressors
#' Reg_ID <- rbind(
#'   c(0,0,0,0,0),  # null
#'   c(1,1,0,0,0),  # two from group 1 -> penalty p_1^(1)
#'   c(1,1,1,0,0),  # three from group 1 -> penalty p_1^(2)
#'   c(1,1,0,1,1)   # two from group 1 and two from group 2 -> p_1^1 * p_2^1
#' )
#' Nar_vec <- c(1,1,1,2,2)
#'
#' # Scalar p (same for all groups)
#' group_dilution(Reg_ID, Nar_vec, p = 0.7)
#'
#' # Group-specific p (vector in group order: group 1 then group 2)
#' group_dilution(Reg_ID, Nar_vec, p = c(0.7, 0.5))
#'
#' # Group-specific p with explicit mapping via names
#' group_dilution(Reg_ID, Nar_vec, p = c("1"=0.7, "2"=0.5))
#'
#' @export
group_dilution <- function(Reg_ID, Nar_vec, p) {

  Reg_ID <- as.matrix(Reg_ID)
  Nar_vec <- as.integer(Nar_vec)

  MS <- nrow(Reg_ID)
  K  <- ncol(Reg_ID)

  if (length(Nar_vec) != K) {
    stop("Nar_vec must have length equal to ncol(Reg_ID) = K.")
  }
  if (!is.numeric(p) || anyNA(p)) {
    stop("p must be numeric and not NA.")
  }
  if (any(p < 0 | p > 1)) {
    stop("All p values must lie in [0, 1].")
  }

  groups <- sort(unique(Nar_vec[Nar_vec > 0L]))
  G <- length(groups)

  if (G == 0L) {
    return(rep(1, MS))  # no groups => no dilution
  }

  # Build a vector p_by_group of length G aligned with `groups`
  if (length(p) == 1L) {
    p_by_group <- rep(p, G)
  } else if (length(p) == G) {
    if (!is.null(names(p))) {
      # Match by group IDs via names
      gp_names <- as.character(groups)
      if (!all(gp_names %in% names(p))) {
        stop("Named p provided, but not all group IDs are present in names(p).")
      }
      p_by_group <- as.numeric(p[gp_names])
    } else {
      # Assume p is in the same order as sorted group IDs
      p_by_group <- as.numeric(p)
    }
  } else {
    stop(sprintf(
      "p must be length 1 or length equal to the number of groups (%d).",
      G
    ))
  }

  # Total log-dilution per model (use logs for numerical stability)
  log_dil <- numeric(MS)

  for (gi in seq_along(groups)) {
    g <- groups[gi]
    cols <- which(Nar_vec == g)

    # counts of included regressors from group g
    cnt <- as.integer(Reg_ID[, cols, drop = FALSE] %*% rep(1L, length(cols)))

    expo <- pmax(cnt - 1L, 0L)

    # add expo * log(p_g), with careful handling of p_g = 0
    pg <- p_by_group[gi]
    if (pg == 0) {
      # if expo > 0 => log dilution is -Inf => dilution 0
      log_dil[expo > 0L] <- -Inf
      # expo == 0 contributes 0
    } else {
      # if already -Inf, keep -Inf
      ok <- is.finite(log_dil)
      log_dil[ok] <- log_dil[ok] + expo[ok] * log(pg)
    }
  }

  # convert back
  out <- exp(log_dil)
  out[!is.finite(log_dil)] <- 0
  out
}
