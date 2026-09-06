# S3 methods for the model_space and bma objects.

fmt_g <- function(info) {
  if (is.character(info$g)) return(info$g)
  lab <- info$g_label
  if (is.null(lab) || identical(lab, "user-specified"))
    format(info$g, digits = 4)
  else paste0(lab, " (g = ", format(info$g, digits = 4), ")")
}

# pad labels so the colons line up, as badpse's print method does
kv <- function(label, ..., width = 20) {
  cat("  ", formatC(label, width = -width), ": ", ..., "\n", sep = "")
}

#' Compactly describe a model space
#'
#' @param x A \code{model_space} object.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @export
print.model_space <- function(x, ...) {
  info <- x$info
  mc3  <- identical(info$method, "mc3")
  n <- function(v) format(v, big.mark = ",", scientific = FALSE)
  cat("model space\n")
  kv("regressors", x$K)
  kv("maximum size", x$M, if (x$M == x$K) " (full model space)" else "")
  if (mc3) kv("models", n(x$MS), " visited by MC3, of ", n(info$space_size))
  else     kv("models", n(x$MS), " enumerated")
  kv("g prior", fmt_g(info))
  kv("covariance", if (isTRUE(info$HC)) "heteroscedasticity-consistent"
                   else "conventional")
  if (mc3) {
    kv("chain", info$draws, " draws after ", info$burn, " burn-in")
    kv("acceptance rate", format(info$acceptance, digits = 3))
    kv("convergence", "cor(analytic PMP, visit frequency) = ",
       if (is.na(info$cor_pmp)) "NA" else format(info$cor_pmp, digits = 4))
  }
  invisible(x)
}

#' Summarize a model space
#'
#' Reports how the model space was built and how the models in it are
#' distributed by size. For an MC^3 space the distribution is over the models
#' actually visited, weighted by how often the chain visited them.
#'
#' @param object A \code{model_space} object.
#' @param ... Ignored.
#' @return An object of class \code{summary.model_space}, invisibly printed.
#' @export
summary.model_space <- function(object, ...) {
  info <- object$info
  K    <- object$K
  r    <- rowSums(object$ols_results[, seq_len(K), drop = FALSE])
  w    <- if (identical(info$method, "mc3")) info$visits / sum(info$visits) else NULL

  sizes <- data.frame(size = 0:object$M)
  sizes$models <- as.integer(table(factor(r, levels = 0:object$M)))
  if (!is.null(w)) {
    sizes$visit_share <- as.numeric(tapply(w, factor(r, levels = 0:object$M), sum))
    sizes$visit_share[is.na(sizes$visit_share)] <- 0
  }

  ll <- object$ols_results[, 3 * K + 3]
  structure(list(K = K, M = object$M, MS = object$MS, info = info,
                 sizes = sizes,
                 loglik_range = range(ll[is.finite(ll)]),
                 best = object$x_names[object$ols_results[which.max(ll),
                                                          seq_len(K)] == 1]),
            class = "summary.model_space")
}

#' @param x A \code{summary.model_space} object.
#' @param ... Ignored.
#' @rdname summary.model_space
#' @export
print.summary.model_space <- function(x, ...) {
  print.model_space(list(K = x$K, M = x$M, MS = x$MS, info = x$info))
  cat("\n  log marginal likelihood: from",
      format(x$loglik_range[1], digits = 6), "to",
      format(x$loglik_range[2], digits = 6), "\n")
  cat("  highest-likelihood model:",
      if (length(x$best)) paste(x$best, collapse = ", ") else "(intercept only)", "\n")
  cat("\nmodels by size\n")
  print(x$sizes, row.names = FALSE)
  invisible(x)
}

#' Compactly describe a BMA object
#'
#' @param x A \code{bma} object.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @export
print.bma <- function(x, ...) {
  mc3 <- !is.null(x[[16]])
  pms <- x[[4]][, 2]
  cat("Bayesian model averaging\n")
  kv("regressors", x[[6]])
  kv("models", format(x[[7]], big.mark = ",", scientific = FALSE),
     if (mc3) " visited by MC3" else " enumerated")
  kv("prior model size", format(x[[8]], digits = 4))
  kv("posterior model size", format(pms[1], digits = 4), " (binomial), ",
     format(pms[2], digits = 4), " (binomial-beta)")
  kv("dilution prior", if (identical(x[[9]], 1)) "applied" else "none")
  kv("extreme bounds", if (is.null(x[[3]])) "not available (MC3 model space)"
                       else "in element 3")
  cat("\nUse summary() for the coefficient table, coef() for posterior means.\n")
  invisible(x)
}

#' Summarize a BMA object
#'
#' @param object A \code{bma} object.
#' @param prior Which model prior to report, \code{"uniform"} for the binomial
#'   prior or \code{"random"} for the binomial-beta prior.
#' @param ... Ignored.
#' @return An object of class \code{summary.bma}, invisibly printed.
#' @export
summary.bma <- function(object, prior = c("uniform", "random"), ...) {
  prior <- match.arg(prior)
  tab   <- if (prior == "uniform") object[[1]] else object[[2]]
  ord   <- order(tab[, "PIP"], decreasing = TRUE, na.last = TRUE)
  structure(list(table = tab[ord, , drop = FALSE], prior = prior,
                 K = object[[6]], MS = object[[7]], EMS = object[[8]],
                 pms = object[[4]], dilution = object[[9]],
                 eba = !is.null(object[[3]]), mc3 = object[[16]]),
            class = "summary.bma")
}

#' @param x A \code{summary.bma} object.
#' @param ... Ignored.
#' @rdname summary.bma
#' @export
print.summary.bma <- function(x, ...) {
  cat("Bayesian model averaging --",
      if (x$prior == "uniform") "binomial model prior"
      else "binomial-beta model prior", "\n")
  cat("  regressors:", x$K, "  models:",
      format(x$MS, big.mark = ",", scientific = FALSE),
      if (!is.null(x$mc3)) "(MC3)" else "(enumerated)", "\n")
  cat("  prior model size:", format(x$EMS, digits = 4),
      "  posterior model size:",
      format(x$pms[if (x$prior == "uniform") 1 else 2, 2], digits = 4), "\n")
  if (!is.null(x$mc3))
    cat("  chain: acceptance", format(x$mc3$acceptance, digits = 3),
        " cor(analytic PMP, visit frequency)",
        format(x$mc3$cor_pmp, digits = 4), "\n")
  if (!x$eba)
    cat("  extreme bounds analysis is not available for an MC3 model space\n")
  cat("\nregressors ordered by posterior inclusion probability\n")
  print(x$table)
  invisible(x)
}

#' Posterior means from a BMA object
#'
#' @param object A \code{bma} object.
#' @param prior Which model prior to use, \code{"uniform"} for the binomial
#'   prior or \code{"random"} for the binomial-beta prior.
#' @param conditional If \code{TRUE}, return posterior means conditional on
#'   inclusion (\code{PMcon}) rather than unconditional ones (\code{PM}).
#' @param ... Ignored.
#' @return A named numeric vector, the intercept first.
#' @export
coef.bma <- function(object, prior = c("uniform", "random"),
                     conditional = FALSE, ...) {
  prior <- match.arg(prior)
  tab   <- if (prior == "uniform") object[[1]] else object[[2]]
  cl    <- if (isTRUE(conditional)) "PMcon" else "PM"
  stats::setNames(as.numeric(tab[, cl]), rownames(tab))
}
