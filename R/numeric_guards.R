#' Snap a rounding-level sum of squares to exactly zero
#'
#' A perfectly fitting model has a residual sum of squares of zero in exact
#' arithmetic, but whether it comes out as exactly zero or as a rounding-level
#' residual depends on the BLAS: an error of 1e-15 in an estimated coefficient
#' is enough to turn a log likelihood of Inf into a large finite number. The
#' quantities built as \code{yty - m * mean(y)^2} can even come out slightly
#' negative through cancellation, which sends \code{log()} to NaN.
#'
#' Snapping anything at or below rounding level to zero makes the degenerate
#' case behave identically on every platform. The tolerance scales with the
#' data, so it cannot affect a genuine fit: for a sum of squares to be caught
#' it must be smaller than the rounding error of the data it was computed from.
#'
#' @param x A sum of squares, as a scalar or 1x1 matrix.
#' @param scale The scale of the data it came from, typically \code{sum(y^2)}.
#' @param m Number of observations.
#' @return \code{x} as a scalar, or 0 when it is at or below rounding level.
#' @keywords internal
#' @noRd
snap_zero <- function(x, scale, m) {
  x <- as.numeric(x)
  if (is.finite(x) && x <= .Machine$double.eps * m * max(1, abs(scale))) 0 else x
}
