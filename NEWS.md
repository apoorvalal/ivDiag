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
* Fixed the parallel bootstrap in `ivDiag()`: since the GitHub change of July 2024
  (#12) the worker cluster was stopped twice and the call failed with "invalid
  connection". The cluster is now stopped once on exit, and the sequential
  `foreach` backend is restored. (CRAN 1.0.6 was not affected.)
* Added unit tests (`tests/testthat`).

# ivDiag 1.0.6
Fix clustering issues with both the main function and AR test (thanks to Michael
Allen @mmooddaa for reporting this bug).

# ivDiag 1.0.5
Fix bugs with `plot_coef()`.

# ivDiag 1.0.4
Allow users not to run bootstraps in `ivDiag`. First CRAN version.

# ivDiag 1.0.2
First complete version ready for CRAN submission.

