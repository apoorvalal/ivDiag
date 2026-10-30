# ivDiag 1.0.8
* `eff_F()` (and so the effective F in `ivDiag()`) now builds the instrument cross-product
  matrix from the instruments with the intercept, controls and fixed effects partialled
  out, weighted by `weights`. Earlier versions used the raw instruments (demeaned only
  within fixed effects), so with two or more instruments, controls, no fixed effects or
  weights gave a wrong value, e.g. 24.5452 instead of 85.7537 for two instruments with
  nonzero means and no controls. Results with one instrument are unchanged.
* `eff_F()` (and so `ivDiag()`) now stops with a clear error when the instruments are
  collinear with each other or with the controls or fixed effects, instead of returning
  `NaN`, `NA` or a wrong number. Collinear controls are still handled. `ivDiag()` checks
  this first, so it no longer fails with a matrix-inversion error.
* `eff_F()` (and so `ivDiag()`) no longer fails or returns a wrong value when instruments
  are on very different scales (e.g. 1e6 and 1e-6) or have a very large mean (1e5 or
  more). The effective F does not depend on the scale or mean of the instruments, and
  they are now standardized internally.

# ivDiag 1.0.7
* `AR_test()` (and `ivDiag()`) now find the Anderson-Rubin confidence set on the
  whole real line. The reduced-form coefficients are linear and their robust or
  clustered variance quadratic in the hypothesised effect, so with one instrument
  the set solves a quadratic inequality (an interval, two rays, the whole line or
  empty); with several instruments the whole line is scanned and each boundary is
  refined by root finding. Earlier versions inverted the test only on the 2SLS
  estimate +/- 10 SEs and could miss pieces of the set outside that range, e.g.
  reporting `[0.0642, Inf)` instead of `(-Inf, -1.1605] Union [0.0640, Inf)`.
  The AR test of zero effect (`Fstat`) is unchanged. `parallel` and `cores` in
  `AR_test()` are no longer used.
* `ivDiag()` gains `running`, `cutoff`, `rd.poly`, `bandwidth` and `kernel` for
  fuzzy regression discontinuity designs: the instrument is the indicator for
  crossing the cutoff, and a polynomial in the running variable, separate on each
  side, is partialled out, so only the jump at the cutoff identifies the effect.
* `ltz()` with `weights`: the adjustment term A now uses weighted cross-products,
  consistent with the weighted 2SLS fit (with one instrument, A = 1 / the weighted
  first-stage coefficient). Earlier versions used unweighted cross-products of the
  weighted residuals, e.g. LTZ estimate 3.6088 (SE 0.8112) instead of 3.7441 (SE
  0.7568) in the GSZ example. Results without weights are unchanged (#18).
* The parallel bootstrap starts its workers with `parallelly::makeClusterPSOCK()`;
  `future::makeClusterPSOCK()` is a re-export that will be removed (#15). **future**
  is no longer imported; **parallelly** is.
* New dataset `mcclendon` (McClendon 2014, a randomized encouragement design; CC0),
  now the main example in the tutorial and the help pages. The tutorial adds a weak-
  instrument AR example and a section on fuzzy RD designs. `rueda` is still shipped.
* The bootstrap is now reproducible: each draw gets its own random-number stream
  derived from `seed` (**doRNG**, a new import), so results are identical with or
  without parallel computing and for any number of cores, and the caller's
  random-number state is left as it was. Up to 1.0.6, `seed` did not reach the
  parallel workers, so the default (parallel) bootstrap changed from run to run.
  For a given `seed`, bootstrap results differ from those of 1.0.6.
* `tF()` below F = 4, where the Lee et al. (2022) table starts: the critical value is
  now `Inf` and the interval the whole real line (p-value 1). Earlier versions used
  18.66, the value at F = 4, for any smaller F, so the interval was too narrow; the
  valid critical value is infinite for F up to 3.84. `plot_coef()` labels such an
  interval "Unbounded CI".
* Silenced harmless "NaNs produced" warnings from felm's summary statistics in
  bootstrap resamples (e.g. the first stage under one-sided noncompliance).
* Fixed the parallel bootstrap in `ivDiag()`: since the GitHub change of July 2024
  (#12) the worker cluster was stopped twice and the call failed with "invalid
  connection". The cluster is now stopped once on exit, and the sequential
  `foreach` backend is restored. (CRAN 1.0.6 was not affected.)
* Added unit tests (`tests/testthat`).
* DESCRIPTION gives the GitHub address for bug reports (#8, #16).
* The tutorial is now a Quarto book in `tutorial/` (not part of the package build).
* CRAN housekeeping: `ivDiag()` no longer calls `set.seed()`, so the caller's
  random-number state is unchanged; unused imports ('fixest', 'glue') removed;
  'testthat' moved to Suggests and dropped from the examples; base packages declared;
  DESCRIPTION cites Lal et al. (2024, Political Analysis) by DOI; data help pages give
  their Dataverse sources and licenses.

# ivDiag 1.0.6
Fix clustering issues with both the main function and AR test (thanks to Michael
Allen @mmooddaa for reporting this bug).

# ivDiag 1.0.5
Fix bugs with `plot_coef()`.

# ivDiag 1.0.4
Allow users not to run bootstraps in `ivDiag`. First CRAN version.

# ivDiag 1.0.2
First complete version ready for CRAN submission.

