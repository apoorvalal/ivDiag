# felm() can warn "NaNs produced" for a coefficient whose variance is zero in
# theory (the intercept under one-sided noncompliance) but rounds to slightly
# below zero with some BLAS/LAPACK libraries (CRAN's BLIS check of 1.0.7).
# quiet_nan() muffles that warning and nothing else.
test_that("quiet_nan muffles only the NaN warning", {
  expect_no_warning(ivDiag:::quiet_nan(sqrt(-1)))
  expect_true(is.nan(ivDiag:::quiet_nan(sqrt(-1))))
  expect_warning(ivDiag:::quiet_nan(warning("something else")), "something else")
})

test_that("quiet_nan matches the translated message", {
  skip_on_cran()
  old <- Sys.setLanguage("fr")
  on.exit(Sys.setLanguage(old), add = TRUE)
  expect_no_warning(ivDiag:::quiet_nan(sqrt(-1)))
})
