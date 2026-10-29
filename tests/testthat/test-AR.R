# AR confidence sets: the set must match a direct check with the regression
# itself, on the whole real line, not only near the 2SLS estimate.

# accept/reject of beta by the regression itself (no shortcut)
direct_accept <- function(data, Y, D, Z, cl = NULL, alpha = 0.05) {
  d <- data[, c(Y, D, Z, cl)]
  for (v in c(Y, D, Z)) d[, v] <- d[, v] - mean(d[, v])
  function(b) unname(ivDiag:::AR_one(b, d, Y, D, Z, cl, NULL)[4] >= alpha)
}
in_set <- function(ci, b) {
  if (length(ci) == 1 && is.na(ci)) return(FALSE)
  any(b >= ci[seq(1, length(ci), 2)] & b <= ci[seq(2, length(ci), 2)])
}

sim <- function(n, pi, seed, p = 1) {
  set.seed(seed)
  z <- matrix(rnorm(n * p), n, p)
  u <- rnorm(n)
  d <- c(z %*% rep(pi, p)) + 0.8 * u + rnorm(n)
  y <- 1 * d + u
  out <- data.frame(y = y, d = d, z, g = rep(1:(n / 5), each = 5))
  names(out)[3:(2 + p)] <- paste0("z", 1:p)
  out
}

test_that("strong instrument: bounded interval, ends where p = alpha", {
  df <- sim(500, 1, 1)
  ar <- suppressMessages(AR_test(df, "y", "d", "z1", controls = NULL))
  expect_true(ar$bounded)
  expect_length(ar$ci, 2)
  acc <- direct_accept(df, "y", "d", "z1")
  # just inside and just outside each end
  expect_true(acc(ar$ci[1] + 1e-3) && acc(ar$ci[2] - 1e-3))
  expect_false(acc(ar$ci[1] - 1e-3) || acc(ar$ci[2] + 1e-3))
  expect_match(ar$ci.print, "^\\[-?[0-9.]+, -?[0-9.]+\\]$")
})

test_that("clustered strong instrument: bounded, ends where p = alpha", {
  df <- sim(500, 1, 2)
  ar <- suppressMessages(AR_test(df, "y", "d", "z1", controls = NULL, cl = "g"))
  acc <- direct_accept(df, "y", "d", "z1", cl = "g")
  expect_true(ar$bounded)
  expect_true(acc(ar$ci[1] + 1e-3) && acc(ar$ci[2] - 1e-3))
  expect_false(acc(ar$ci[1] - 1e-3) || acc(ar$ci[2] + 1e-3))
})

test_that("no identification: the whole real line", {
  df <- sim(500, 1, 3)
  # make z exactly orthogonal to both d and y in the sample
  df$z1 <- resid(lm(z1 ~ d + y, data = df))
  # lfe warns that the 2SLS fit is rank deficient, as it should here
  ar <- suppressWarnings(suppressMessages(AR_test(df, "y", "d", "z1", controls = NULL)))
  expect_equal(ar$ci, c(-Inf, Inf))
  expect_equal(ar$ci.print, "(-Inf, Inf)")
  expect_false(ar$bounded)
})

test_that("weak instrument: every piece of the set, including far rays", {
  for (seed in 1:20) {
    df <- sim(200, 0.08, seed)
    ar <- suppressMessages(AR_test(df, "y", "d", "z1", controls = NULL))
    acc <- direct_accept(df, "y", "d", "z1")
    fin <- ar$ci[is.finite(ar$ci)]
    probes <- c(-1e6, -1e3, -10, 0, 10, 1e3, 1e6, fin - 1e-4, fin + 1e-4)
    for (b in probes) expect_equal(acc(b), in_set(ar$ci, b), info = paste("seed", seed, "b", b))
  }
})

test_that("two instruments: set matches the direct check", {
  for (seed in 1:5) {
    df <- sim(300, c(0.1), seed, p = 2)
    ar <- suppressMessages(AR_test(df, "y", "d", c("z1", "z2"), controls = NULL))
    acc <- direct_accept(df, "y", "d", c("z1", "z2"))
    fin <- ar$ci[is.finite(ar$ci)]
    probes <- c(-1e6, -1e2, -3, 0, 1, 3, 1e2, 1e6, fin - 1e-4, fin + 1e-4)
    for (b in probes) expect_equal(acc(b), in_set(ar$ci, b), info = paste("seed", seed, "b", b))
  }
})

# fixtures/hager_krakowski_2022.rds: the three variables of Hager and Krakowski (2022, APSR),
# column 1, as standardized in Lal et al. (2024); source data CC0 1.0, doi:10.7910/DVN/JWTRSC.
# ivDiag 1.0.6 searched only near the 2SLS estimate and missed the left ray.
test_that("Hager and Krakowski (2022), column 1: two rays", {
  df <- readRDS(test_path("fixtures", "hager_krakowski_2022.rds"))
  ar <- suppressMessages(AR_test(df, "y", "commanders", "priests_continuous", controls = NULL))
  expect_equal(ar$ci[c(1, 4)], c(-Inf, Inf))
  expect_equal(ar$ci[2], -1.16, tolerance = 0.005 / 1.16)
  expect_equal(ar$ci[3], 0.064, tolerance = 0.001 / 0.064)
  expect_false(ar$bounded)
  expect_match(ar$ci.print, "^\\(-Inf, -1\\.16[0-9]*\\] Union \\[0\\.06[0-9]*, Inf\\)$")
})
