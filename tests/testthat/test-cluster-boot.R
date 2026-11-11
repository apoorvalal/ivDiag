# The cluster bootstrap must draw whole clusters that appear in the data. Unused
# levels of a factor clustering variable used to shift the lookup from a sampled
# cluster to its rows, so draws fetched empty groups or another cluster's rows
# (issue #24).

sim_cl <- function() {
  set.seed(127)
  d <- data.frame(group = rep(c("b", "d", "f", "h", "j"), each = 16),
    z = rnorm(80), x = rnorm(80))
  d$d <- .8 * d$z + .4 * d$x + rnorm(80)
  d$y <- 1.3 * d$d + .5 * d$x + rnorm(80)
  d
}

run_cl <- function(dat, ...) {
  suppressMessages(ivDiag(data = dat, Y = "y", D = "d", Z = "z", controls = "x",
    cl = "group", run.AR = FALSE, nboots = 100, seed = 94305, ...))
}

boot_parts <- function(g) list(g$est_ols, g$est_2sls, g$est_rf, g$est_fs, g$F_stat)

test_that("each cluster maps to its own rows, with or without unused levels", {
  g <- factor(c("b", "b", "d", "d", "d"), levels = c("a", "b", "c", "d", "e"))
  cg <- cluster_groups(g)
  expect_equal(as.character(cg$clusters), c("b", "d"))
  expect_true(all(lengths(cg$id.list) > 0))
  for (k in seq_along(cg$clusters)) {
    expect_true(all(g[cg$id.list[[k]]] == cg$clusters[k]))
  }
  expect_equal(sort(unlist(cg$id.list[match(cg$clusters, cg$clusters)],
    use.names = FALSE)), 1:5)
  # numeric identifiers that print alike stay separate clusters
  num <- cluster_groups(c(0.1 + 0.2, 0.3, 0.3))
  expect_equal(length(num$clusters), 2)
  expect_equal(unname(lengths(num$id.list)), c(1, 2))
})

test_that("unused factor levels do not change the bootstrap", {
  d <- sim_cl()
  base <- run_cl(transform(d, group = factor(group)), parallel = FALSE)
  # unused levels before, between and after the observed ones
  wide <- transform(d, group = factor(group, levels = letters[1:10]))
  expect_equal(boot_parts(run_cl(wide, parallel = FALSE)), boot_parts(base))
  # a character clustering variable gives the same
  expect_equal(boot_parts(run_cl(d, parallel = FALSE)), boot_parts(base))
  # the case of issue #24, before the fix: SE 0.0600 with unused levels
  expect_equal(unname(base$est_2sls["Boot.c", "SE"]), 0.1327, tolerance = 1e-3)
})

test_that("a cluster removed by the complete-case filter is not drawn", {
  d <- sim_cl()
  d$group <- factor(d$group)
  d$y[d$group == "f"] <- NA # every row of cluster f is incomplete
  kept <- d[d$group != "f", ]
  kept$group <- droplevels(kept$group)
  expect_equal(boot_parts(run_cl(d, parallel = FALSE)),
    boot_parts(run_cl(kept, parallel = FALSE)))
})

test_that("parallel bootstrap: unused levels do not change it either", {
  skip_on_cran()
  d <- sim_cl()
  wide <- transform(d, group = factor(group, levels = letters[1:10]))
  par <- run_cl(wide, cores = 2)
  expect_equal(boot_parts(par), boot_parts(run_cl(d, parallel = FALSE)))
})
