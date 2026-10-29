## Submission

ivDiag 1.0.7 (update of 1.0.6).

* The Anderson-Rubin confidence set is now computed on the whole real line (earlier
  versions searched only near the 2SLS estimate and could miss unbounded pieces).
* The bootstrap is reproducible: random-number streams per draw via 'doRNG' (new
  import); the caller's random-number state is no longer changed (no `set.seed()`
  inside functions).
* Fixed the parallel bootstrap (cluster stopped twice), added a fuzzy regression
  discontinuity option, a new example dataset (`mcclendon`), and unit tests.
* `tF()` returns an unbounded interval below F = 4, where the Lee et al. (2022) table starts.
* Fixed `ltz()` with weights: the adjustment now uses the weighted first stage.
* The parallel bootstrap now uses 'parallelly' (replacing 'future', whose
  makeClusterPSOCK() is a re-export).
* Removed unused imports ('fixest', 'glue'); 'testthat' moved from Imports to Suggests.

## Test environments

* local: macOS (aarch64), R 4.5.3
* TODO before submitting: win-builder (R-devel), macbuilder

## R CMD check results

0 errors | 0 warnings | 0 notes

## Notes for the reviewer

* Examples and tests run serially or on at most 2 cores; the parallel bootstrap test is
  skipped on CRAN.
* Bundled data: `mcclendon` and `rueda` come from Harvard Dataverse deposits released
  under CC0 1.0 (doi:10.7910/DVN/LVOSAF, doi:10.7910/DVN/K6ZOOW); sources are in the
  help pages. `gsz` and `gsz_south` (Guiso, Sapienza and Zingales 2016) were shipped
  in earlier versions.

## Reverse dependencies

None on CRAN.
