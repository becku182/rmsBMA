# rmsBMA (development version)

## Bug fixes

* The list returned by `bma()` named its first two elements identically, both
  "Table with the binomial model prior results". The second element holds the
  binomial-beta results and is now named accordingly. Only the name changes;
  the contents were always correct.

# rmsBMA 0.2.0

## New features

* `model_space()` gains `mc3`, `draws` and `burn`. With `mc3 = TRUE` the model
  space is explored by MC^3 sampling (Madigan and York, 1995) instead of
  exhaustive enumeration, which makes a large number of regressors feasible:
  enumeration means fitting `2^K` models, already 131,072 at `K = 17`. Reduced
  model spaces are supported, with the proposal correction the size constraint
  requires. The chain reports its acceptance rate, the number of distinct
  models visited and a convergence diagnostic, and warns when the diagnostic
  falls below 0.99 or when the chain has barely moved.

* Extreme Bounds Analysis is withheld for a sampled model space. Extreme bounds
  are the smallest and largest coefficient estimates across all models, and a
  sampler does not visit the extremes of a distribution it is exploring by
  posterior mass, so bounds from a visited subset are systematically too
  narrow. Element 3 of the `bma` object is `NULL` and a message explains why.

* `print`, `summary` and `coef` methods for the `model_space` and `bma`
  objects. `summary()` on a `bma` object gives the coefficient table ordered by
  posterior inclusion probability, for either model prior; `coef()` returns the
  posterior means as a named vector, optionally conditional on inclusion.

* `model_sizes()` and `model_pmp()` gain `type`, selecting between the existing
  line display and a histogram.

## Bug fixes

* The estimators no longer depend on the BLAS for their treatment of a
  perfectly fitting model. With a constant dependent variable the residual sum
  of squares is zero in exact arithmetic, but whether it came out as exactly
  zero or as a rounding-level residual decided whether the log marginal
  likelihood was `Inf` or a large finite number. Rounding-level sums of squares
  are now snapped to zero, so the degenerate case behaves identically on any
  BLAS. Reported by CRAN's BLIS check machine.

* The same guard removes a cancellation in the g-prior estimators. The centered
  sum of squares was computed as a difference of two nearly equal numbers and
  could come out slightly negative, which would have produced a `NaN` log
  marginal likelihood and propagated silently into the posterior model
  probabilities.

* `coef_hist()`, `model_pmp()`, `model_sizes()` and `posterior_dens()` now have
  `\usage` sections in their help files. A `utils::globalVariables()` call
  between the roxygen block and the function had attached the block to that
  call, so no usage could be derived. Reported as a NOTE by CRAN's r-devel
  check machines.

* `model_pmp()` misspelled "ranking" on the x axis of two of its three
  graphs, and reported the wrong quantity when `top`
  exceeded the size of the model space.

## Other

* An enumerated `model_space` object now always carries a sixth element,
  `info`, recording how the space was built. Elements 1 to 5 are unchanged.
  The `bma` object gains a sixteenth element holding the chain diagnostics,
  `NULL` for an enumerated model space.

* Results obtained from an enumerated model space are unchanged from 0.1.2.

# rmsBMA 0.1.2

* First CRAN release.
