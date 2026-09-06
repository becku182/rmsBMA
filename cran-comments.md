## R CMD check results

0 errors | 0 warnings | 0 notes

## Previous check issues addressed

This version resolves both issues reported for 0.1.2.

* **Rd files without \usage** (NOTE on r-devel-linux-x86_64-debian-clang and
  r-devel-linux-x86_64-debian-gcc). Four files carried a
  `utils::globalVariables()` call between the roxygen block and the function it
  documented, so the block was attached to that call rather than to the
  function and no \usage could be derived. The declarations are now
  consolidated in `R/globals.R`, and `coef_hist.Rd`, `model_pmp.Rd`,
  `model_sizes.Rd` and `posterior_dens.Rd` all carry \usage sections.

* **Test failure on the BLIS BLAS machine** (ERROR in
  `test-fast_ols_const.R`). With a constant dependent variable the fit is
  exact, so the residual sum of squares is zero in exact arithmetic and the log
  marginal likelihood is Inf. Whether the BLAS returns exactly zero or a
  rounding-level residual decided whether the result was Inf or a large finite
  number. Every estimator now snaps a rounding-level sum of squares to zero, so
  the degenerate case behaves identically on any BLAS. The tolerance scales
  with the data and cannot affect a genuine fit. The same guard also removes a
  cancellation in the g-prior estimators that could make a sum of squares
  slightly negative and send log() to NaN.

## Changes in this version

* MC^3 model-space sampling (`model_space(mc3 = TRUE)`), for full and reduced
  model spaces, with chain diagnostics and convergence warnings. Extreme Bounds
  Analysis is withheld for a sampled model space, since a sampler does not
  visit the extremes.
* `print`, `summary` and `coef` methods for the `model_space` and `bma`
  objects.
* A `type` argument on `model_sizes()` and `model_pmp()` selecting between the
  existing line display and a histogram.

Results from the enumerated model space are unchanged from 0.1.2.

## Downstream dependencies

There are no downstream dependencies.
