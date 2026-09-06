## Submission

This is a feature and maintenance release. It resolves both issues reported by
the CRAN check machines for version 0.1.2.

## Test environments

* local: macOS 15.7.4, R 4.4.1, aarch64-apple-darwin20
* win-builder: R-devel                    [PENDING -- paste result]
* win-builder: R release                  [PENDING -- paste result]
* macOS builder: R release                [PENDING -- paste result]

## R CMD check results

Local check: 0 errors | 0 warnings | 2 notes.

* `checking for future file timestamps ... NOTE -- unable to verify current
  time`. The check machine could not reach the time server. Not reproducible
  with network access.

* `checking HTML version of manual ... NOTE`. The `tidy` binary shipped with
  macOS predates HTML5 and does not recognise the `<main>` element that R's own
  `Rd2HTML` emits. The same four complaints (`<main>`, `<link>` type attribute,
  `<script>` onload attribute, `<table>` summary attribute) are produced for
  every Rd file in the package and originate in R's HTML template rather than in
  the package. Two data help files also contain en dashes in year and page
  ranges, which that version of `tidy` misreads as Latin-1. The file encoding is
  declared as UTF-8 and `checking Rd files` passes. Not reproduced on a machine
  with a current `tidy`.

## Previous check issues addressed

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
  marginal likelihood is Inf. Whether the BLAS returned exactly zero or a
  rounding-level residual decided whether the result was Inf or a large finite
  number. Every estimator now snaps a rounding-level sum of squares to zero, so
  the degenerate case behaves identically on any BLAS. The tolerance scales with
  the data and cannot affect a genuine fit. The same guard also removes a
  cancellation in the g-prior estimators that could make a sum of squares
  slightly negative and send log() to NaN. Verified by running the test suite
  against BLIS as well as the reference BLAS.

## Changes in this version

* MC^3 model-space sampling (`model_space(mc3 = TRUE)`), for full and reduced
  model spaces, with chain diagnostics and convergence warnings. Extreme Bounds
  Analysis is withheld for a sampled model space, since a sampler does not
  visit the extremes.
* `print`, `summary` and `coef` methods for the `model_space` and `bma`
  objects.
* A `type` argument on `model_sizes()` and `model_pmp()` selecting between the
  existing line display and a histogram.

Results obtained from an enumerated model space are unchanged from 0.1.2.

## Downstream dependencies

There are none.
