# Bootstrap results depend only on `seed`: not on parallel computing, the
# number of cores, or the caller's random-number state.
boot_ci <- function(...) {
  data(ivDiag, package = "ivDiag", envir = environment())
  g <- suppressMessages(ivDiag(mcclendon, "intended", "openedesteem", "esteem",
    nboots = 60, run.AR = FALSE, ...))
  g$est_2sls["Boot.c", ]
}

test_that("same seed, same bootstrap, serial; caller's RNG untouched", {
  kind <- RNGkind()
  set.seed(123); a <- boot_ci(parallel = FALSE, seed = 1)
  set.seed(456); b <- boot_ci(parallel = FALSE, seed = 1)
  expect_equal(a, b)
  expect_false(isTRUE(all.equal(a, boot_ci(parallel = FALSE, seed = 2))))
  expect_equal(RNGkind(), kind)
})

test_that("parallel matches serial and does not depend on the number of cores", {
  skip_on_cran()
  s <- boot_ci(parallel = FALSE, seed = 1)
  expect_equal(boot_ci(cores = 2, seed = 1), s)
  expect_equal(boot_ci(cores = 2, seed = 1), boot_ci(cores = 2, seed = 1))
})

test_that("no warnings from the bootstrap with one-sided noncompliance", {
  expect_no_warning(boot_ci(parallel = FALSE, seed = 1))
})
