# The parallel bootstrap must finish, return finite bootstrap results, and
# leave a sequential foreach backend (on master after PR #12 the cluster was
# stopped twice and ivDiag() failed with "invalid connection").
test_that("parallel bootstrap runs twice in a row", {
  skip_on_cran()
  data(ivDiag, package = "ivDiag", envir = environment())
  for (k in 1:2) {
    g <- suppressMessages(ivDiag(rueda, "e_vote_buying", "lm_pob_mesa", "lz_pob_mesa_f",
      c("lpopulation", "lpotencial"), cl = "muni_code", nboots = 40, cores = 2,
      run.AR = FALSE))
    expect_true(all(is.finite(g$est_2sls["Boot.c", c("SE", "CI 2.5%", "CI 97.5%")])))
  }
  expect_equal(foreach::getDoParName(), "doSEQ")
})
