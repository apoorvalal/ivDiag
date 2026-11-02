## Submission

ivDiag 1.0.8 (update of 1.0.7).

* Fixes the test failure in CRAN's BLIS check of 1.0.7 (deadline 2026-11-22). With
  one-sided noncompliance, the first-stage intercept's variance is zero in theory and
  can round to slightly below zero with some BLAS/LAPACK libraries; 'lfe' then warned
  "NaNs produced" when taking its square root. That one warning is now muffled in every
  `felm()` fit (it concerns a standard error the package never uses); other warnings
  pass through. Reported results are unchanged.
* The Description cites the published article, Lal, Lockhart, Xu and Zu (2024)
  <doi:10.1017/pan.2024.2>.
* The effective F statistic is corrected for two or more instruments with controls,
  without fixed effects, or with weights; collinear instruments now give a clear error;
  instruments are standardized internally (not with multiway clustering).
* Factor and character instruments are expanded into dummy variables for every statistic.

## Test environments

* local: macOS (aarch64), R 4.5.3; tests also run with Apple Accelerate and OpenBLAS
  0.3.31 in place of the reference BLAS (all pass, no warnings)
* win-builder: Windows Server 2022 (x86_64), R-devel (2026-10-31 r90820 ucrt)

## R CMD check results

0 errors | 0 warnings | 1 note (on both)

* Days since last update: 1.0.7 was published on 2026-10-30. This release fixes the
  BLIS issue CRAN asked us to correct.

## Notes for the reviewer

* Examples and tests run serially or on at most 2 cores; the parallel bootstrap test is
  skipped on CRAN.
* Bundled data: `mcclendon` and `rueda` come from Harvard Dataverse deposits released
  under CC0 1.0 (doi:10.7910/DVN/LVOSAF, doi:10.7910/DVN/K6ZOOW); sources are in the
  help pages. `gsz` and `gsz_south` (Guiso, Sapienza and Zingales 2016) were shipped
  in earlier versions.

## Reverse dependencies

None on CRAN.
