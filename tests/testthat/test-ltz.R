# ltz(): with one instrument, A = 1 / (first-stage coefficient), so the LTZ
# estimate is the 2SLS estimate minus mu / pi and its variance adds (sd / pi)^2.
# With weights, pi must be the weighted first-stage coefficient.

sim_ltz <- function(n = 1000, seed = 1) {
  set.seed(seed)
  x <- rnorm(n)
  z <- rbinom(n, 1, 0.4)
  u <- rnorm(n)
  # the first stage differs across units, so weighting changes pi
  w <- exp(rnorm(n))
  d <- (0.2 + 0.6 * (w > 1)) * z + 0.5 * x + 0.5 * u + rnorm(n)
  y <- 1 * d + 0.3 * x + u
  data.frame(y = y, d = d, z = z, x = x, w = w)
}

ltz_by_hand <- function(df, weights = NULL, prior) {
  ww <- if (is.null(weights)) NULL else df[, weights]
  fs <- lfe::felm(d ~ z + x, data = df, weights = ww)
  pi_hat <- unname(coef(fs)["z"])
  iv <- lfe::felm(y ~ x | 0 | (d ~ z), data = df, weights = ww)
  b <- unname(tail(iv$coefficients, 1))
  se <- unname(tail(iv$rse, 1))
  c(Coef = b - prior[1] / pi_hat, SE = sqrt(se^2 + (prior[2] / pi_hat)^2))
}

test_that("ltz matches the one-instrument closed form without weights", {
  df <- sim_ltz()
  out <- ltz(df, Y = "y", D = "d", Z = "z", controls = "x", prior = c(0.1, 0.05), prec = 10)
  expect_equal(unname(out$ltz[c("Coef", "SE")]),
               unname(ltz_by_hand(df, prior = c(0.1, 0.05))), tolerance = 1e-8)
})

test_that("ltz uses the weighted first stage with weights", {
  df <- sim_ltz()
  out <- ltz(df, Y = "y", D = "d", Z = "z", controls = "x", weights = "w",
             prior = c(0.1, 0.05), prec = 10)
  expect_equal(unname(out$ltz[c("Coef", "SE")]),
               unname(ltz_by_hand(df, weights = "w", prior = c(0.1, 0.05))), tolerance = 1e-8)
  # the weighted and unweighted adjustments differ in this design
  out0 <- ltz(df, Y = "y", D = "d", Z = "z", controls = "x", prior = c(0.1, 0.05), prec = 10)
  expect_false(isTRUE(all.equal(out$ltz[["Coef"]] - out$iv[["Coef"]],
                                out0$ltz[["Coef"]] - out0$iv[["Coef"]])))
})

test_that("constant weights give the unweighted result", {
  df <- sim_ltz()
  df$one <- 3
  a <- ltz(df, Y = "y", D = "d", Z = "z", controls = "x", weights = "one", prior = c(0.1, 0.05), prec = 10)
  b <- ltz(df, Y = "y", D = "d", Z = "z", controls = "x", prior = c(0.1, 0.05), prec = 10)
  expect_equal(a$ltz, b$ltz, tolerance = 1e-8)
})

test_that("tutorial example (GSZ 2016, population weights)", {
  controls <- c("altitudine", "escursione", "costal", "nearsea",
                "population", "pop2", "gini_land", "gini_income")
  out <- ltz(gsz, Y = "totassoc_p", D = "libero_comune_allnord", Z = "bishopcity",
             controls = controls, weights = "population", prior = c(0.178, 0.137))
  expect_equal(unname(out$iv[c("Coef", "SE")]), c(4.2087, 0.6670), tolerance = 1e-4)
  expect_equal(unname(out$ltz[c("Coef", "SE")]), c(3.7441, 0.7568), tolerance = 1e-4)
})
