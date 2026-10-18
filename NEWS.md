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
* New dataset `mcclendon` (McClendon 2014, a randomized encouragement design; CC0),
  now the main example in the tutorial and the help pages. The tutorial adds a weak-
  instrument AR example and a section on fuzzy RD designs. `rueda` is still shipped.
* The bootstrap is now reproducible: each draw gets its own random-number stream
  derived from `seed` (**doRNG**, a new import), so results are identical with or
  without parallel computing and for any number of cores, and the caller's
  random-number state is left as it was. Up to 1.0.6, `seed` did not reach the
  parallel workers, so the default (parallel) bootstrap changed from run to run.
  For a given `seed`, bootstrap results differ from those of 1.0.6.
* Silenced harmless "NaNs produced" warnings from felm's summary statistics in
  bootstrap resamples (e.g. the first stage under one-sided noncompliance).
* Fixed the parallel bootstrap in `ivDiag()`: since the GitHub change of July 2024
  (#12) the worker cluster was stopped twice and the call failed with "invalid
  connection". The cluster is now stopped once on exit, and the sequential
  `foreach` backend is restored. (CRAN 1.0.6 was not affected.)
* Added unit tests (`tests/testthat`).
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

