## Submission

ivDiag 1.0.9 (update of 1.0.8).

* Fixes a bug that changes results: with a factor clustering variable that has unused
  levels, the cluster bootstrap looked up the wrong rows for the clusters it drew, so
  bootstrap standard errors, intervals, p-values and the bootstrapped F statistic were
  wrong. Reported by a user (GitHub issue #24). Each drawn cluster now gets its own rows;
  results without unused factor levels are unchanged. New tests cover the case.

## Test environments

* local: macOS (aarch64), R 4.5.3
* win-builder: Windows Server 2022 (x86_64), R-devel --- TODO before submitting
* macbuilder: macOS, R-release --- TODO before submitting

## R CMD check results

Local `R CMD check --as-cran`: 0 errors | 0 warnings | 0 notes.
win-builder: TODO; expected 1 note:

* Days since last update: 1.0.8 was published on 2026-11-02. We are sorry for the quick
  update: 1.0.9 fixes a bug that gives users wrong bootstrap results (see above).

## Notes for the reviewer

* Examples and tests run serially or on at most 2 cores; the parallel bootstrap test is
  skipped on CRAN.
* Bundled data: `mcclendon` and `rueda` come from Harvard Dataverse deposits released
  under CC0 1.0 (doi:10.7910/DVN/LVOSAF, doi:10.7910/DVN/K6ZOOW); sources are in the
  help pages. `gsz` and `gsz_south` (Guiso, Sapienza and Zingales 2016) were shipped
  in earlier versions.

## Reverse dependencies

None on CRAN.
