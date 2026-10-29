# tF(): critical values interpolated from the Lee et al. (2022) table by sqrt(F);
# below F = 4 the interval is the whole line.

test_that("tF matches the table at its ends and grid points", {
  expect_equal(tF(coef = 1, se = 0.5, Fstat = 4)[["cF"]], 18.66)
  expect_equal(tF(coef = 1, se = 0.5, Fstat = 2.1^2)[["cF"]], 9.74)
  expect_equal(tF(coef = 1, se = 0.5, Fstat = 10.3^2)[["cF"]], 1.96)
  expect_equal(tF(coef = 1, se = 0.5, Fstat = 500)[["cF"]], 1.96)
  # between grid points: linear in sqrt(F)
  expect_equal(tF(coef = 1, se = 0.5, Fstat = 2.05^2)[["cF"]], (18.66 + 9.74) / 2, tolerance = 1e-4)
})

test_that("cF falls as F rises and the interval is coef +/- cF * se", {
  Fs <- seq(4, 110, length.out = 200)
  cF <- sapply(Fs, function(f) tF(coef = 1, se = 0.5, Fstat = f)[["cF"]])
  expect_true(all(diff(cF) <= 1e-8))
  out <- tF(coef = 0.3149, se = 0.0890, Fstat = 50, prec = 10)
  expect_equal(unname(out[c("CI2.5%", "CI97.5%")]),
               0.3149 + c(-1, 1) * out[["cF"]] * 0.0890, tolerance = 1e-8)
})

test_that("below F = 4 the tF interval is unbounded", {
  for (f in c(0.5, 3.84, 3.9)) {
    out <- tF(coef = 1, se = 0.5, Fstat = f)
    expect_equal(out[["cF"]], Inf)
    expect_equal(unname(out[c("CI2.5%", "CI97.5%")]), c(-Inf, Inf))
    expect_equal(out[["p-value"]], 1)
  }
})

test_that("plot_coef draws a weak-instrument fit with an unbounded tF interval", {
  set.seed(2)
  n <- 300
  z <- rnorm(n); u <- rnorm(n)
  d <- 0.05 * z + 0.8 * u + rnorm(n)
  weak <- data.frame(y = d + u, d = d, z = z)
  g <- suppressMessages(ivDiag(data = weak, Y = "y", D = "d", Z = "z", bootstrap = FALSE))
  expect_lt(g$F_stat[["F.effective"]], 4)
  expect_equal(g$tF[["cF"]], Inf)
  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(plot_coef(g))
})
