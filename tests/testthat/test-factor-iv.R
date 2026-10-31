# A factor (or character) instrument enters as dummies for all levels but the
# base level, so every statistic must equal what the same dummies give when
# they are passed as numeric instruments. Up to 1.0.8 the 2SLS fit and eff_F()
# (through felm) used the dummies, but the standard, robust and cluster F, the
# first stage, the AR test and tF used the factor as one numeric column.

sim_fac <- function(n = 400, seed = 1, K = 4) {
  set.seed(seed)
  x1 <- rnorm(n)
  g <- sample(1:6, n, replace = TRUE)
  ge <- rnorm(6)[g]
  cl <- sample(1:30, n, replace = TRUE)
  f <- sample(letters[1:K], n, replace = TRUE)
  eff <- c(0, 0.5, -0.4, 0.3, 0.6, -0.2)[seq_len(K)]
  u <- rnorm(n)
  d <- eff[match(f, letters)] + 0.5 * x1 + ge + 0.6 * u + rnorm(n)
  y <- d + 0.5 * x1 + ge + u
  df <- data.frame(y = y, d = d, x1 = x1, g = factor(g), cl = cl,
                   w = exp(0.5 * rnorm(n)), zch = f, stringsAsFactors = FALSE)
  # levels in a non-alphabetical order: the base level is the first one, "c"
  df$zf <- factor(f, levels = c("c", setdiff(letters[1:K], "c")))
  # the same instrument as numeric dummies, base level "c"
  for (l in setdiff(letters[1:K], "c")) df[[paste0("D", l)]] <- as.numeric(f == l)
  df
}
dums <- function(K = 4) paste0("D", setdiff(letters[1:K], "c"))

run_iv <- function(df, Z, ...) {
  suppressWarnings(suppressMessages(ivDiag(df, "y", "d", Z, ..., parallel = FALSE)))
}

expect_same_iv <- function(a, b) {
  expect_equal(a$p_iv, b$p_iv)
  expect_equal(a$F_stat, b$F_stat, tolerance = 1e-8)
  expect_equal(a$AR, b$AR, tolerance = 1e-8)
  expect_equal(a$est_2sls, b$est_2sls, tolerance = 1e-8)
  expect_equal(a$est_ols, b$est_ols, tolerance = 1e-8)
  expect_equal(unname(a$est_fs), unname(b$est_fs), tolerance = 1e-8)
  expect_equal(unname(a$est_rf), unname(b$est_rf), tolerance = 1e-8)
  expect_equal(a$rho, b$rho, tolerance = 1e-8)
  expect_equal(a$tF, b$tF, tolerance = 1e-8)
  expect_equal(a$N, b$N)
}

test_that("a factor instrument gives the same results as its dummies", {
  df <- sim_fac()
  fac <- run_iv(df, "zf", bootstrap = FALSE)
  num <- run_iv(df, dums(), bootstrap = FALSE)
  expect_same_iv(fac, num)
  expect_equal(fac$p_iv, 3)
  expect_null(fac$tF)  # tF needs one instrument
  expect_equal(rownames(fac$est_fs), c("zf_a", "zf_b", "zf_d"))
  # F statistics are those of the three dummies, not of one numeric column
  fs <- lfe::felm(d ~ Da + Db + Dd, data = df)
  ftest <- summary(lm(d ~ Da + Db + Dd, data = df))$fstatistic
  expect_equal(unname(fac$F_stat["F.standard"]), unname(ftest["value"]), tolerance = 1e-4)
  # number-like labels, as in Alt (2015): 1.0.8 ran and treated the factor as
  # the numeric codes 1, ..., 4 in the F statistics and the AR test
  df$zn <- factor(as.integer(df$zf))   # "1" is level "c", the base level
  fac2 <- run_iv(df, "zn", bootstrap = FALSE)
  expect_same_iv(fac2, num)
  expect_equal(fac2$p_iv, 3)
})

test_that("factor instrument with controls, fixed effects, clusters and weights", {
  df <- sim_fac(seed = 2)
  fac <- run_iv(df, "zf", controls = "x1", FE = "g", cl = "cl", weights = "w",
                bootstrap = FALSE)
  num <- run_iv(df, dums(), controls = "x1", FE = "g", cl = "cl", weights = "w",
                bootstrap = FALSE)
  expect_same_iv(fac, num)
  expect_true(is.finite(fac$F_stat["F.cluster"]))
})

test_that("a character instrument is expanded like a factor (sorted levels)", {
  df <- sim_fac(seed = 3)
  # base level "a" (alphabetical) for a character column
  for (l in c("b", "c", "d")) df[[paste0("E", l)]] <- as.numeric(df$zch == l)
  chr <- run_iv(df, "zch", controls = "x1", bootstrap = FALSE)
  num <- run_iv(df, c("Eb", "Ec", "Ed"), controls = "x1", bootstrap = FALSE)
  expect_same_iv(chr, num)
  expect_equal(rownames(chr$est_fs), c("zch_b", "zch_c", "zch_d"))
})

test_that("a factor next to a numeric instrument, and the bootstrap", {
  df <- sim_fac(seed = 4)
  df$z2 <- 0.5 * df$d + rnorm(nrow(df))
  fac <- run_iv(df, c("z2", "zf"), controls = "x1", bootstrap = TRUE, nboots = 50)
  num <- run_iv(df, c("z2", dums()), controls = "x1", bootstrap = TRUE, nboots = 50)
  expect_same_iv(fac, num)
  expect_equal(fac$p_iv, 4)
})

test_that("a two-level factor is one instrument, and tF is reported", {
  # values a and b; the declared base level "c" never occurs, so it is dropped
  # and "a" becomes the base level
  df <- sim_fac(seed = 5, K = 2)
  fac <- run_iv(df, "zf", controls = "x1", cl = "cl", bootstrap = FALSE)
  num <- run_iv(df, "Db", controls = "x1", cl = "cl", bootstrap = FALSE)
  expect_same_iv(fac, num)
  expect_equal(fac$p_iv, 1)
  expect_equal(rownames(fac$est_fs), "zf_b")
  expect_false(is.null(fac$tF))
  expect_equal(unname(fac$F_stat["F.effective"]), unname(fac$F_stat["F.cluster"]))
})

test_that("missing values: levels come from the estimation sample", {
  df <- sim_fac(seed = 6)
  df$zf[c(5, 17)] <- NA
  df[c(5, 17), dums()] <- NA
  df$y[c(9, 33)] <- NA
  # a level that occurs only in a row with a missing outcome is not an instrument
  df$zf <- factor(as.character(df$zf), levels = c("c", "a", "b", "d", "e"))
  df$zf[9] <- "e"
  fac <- run_iv(df, "zf", controls = "x1", bootstrap = FALSE)
  num <- run_iv(df, dums(), controls = "x1", bootstrap = FALSE)
  expect_same_iv(fac, num)
  expect_equal(fac$p_iv, 3)
})

test_that("AR_test(), eff_F() and ltz() expand a factor instrument too", {
  df <- sim_fac(seed = 7)
  a1 <- suppressMessages(AR_test(df, "y", "d", "zf", controls = "x1", cl = "cl"))
  a2 <- suppressMessages(AR_test(df, "y", "d", dums(), controls = "x1", cl = "cl"))
  expect_equal(a1, a2, tolerance = 1e-8)
  expect_equal(unname(a1$Fstat["df1"]), 3)
  expect_equal(suppressMessages(eff_F(df, "y", "d", "zf", controls = "x1", FE = "g",
                                      weights = "w", prec = 12)),
               eff_F(df, "y", "d", dums(), controls = "x1", FE = "g", weights = "w",
                     prec = 12), tolerance = 1e-8)
  df2 <- sim_fac(seed = 8, K = 2)
  l1 <- suppressMessages(ltz(df2, "y", "d", "zf", controls = "x1", prior = c(0.1, 0.05)))
  l2 <- ltz(df2, "y", "d", "Db", controls = "x1", prior = c(0.1, 0.05))
  expect_equal(l1, l2)
})

test_that("the expansion is reported once, and a one-level factor stops", {
  df <- sim_fac(seed = 9)
  expect_message(eff_F(df, "y", "d", "zf"), "3 dummy instruments \\(base level \"c\"\\)")
  df$zone <- factor("a")
  expect_error(suppressMessages(ivDiag(df, "y", "d", "zone", bootstrap = FALSE)),
               "fewer than two levels")
  # new names do not overwrite existing columns
  df$zf_a <- 0
  out <- run_iv(df, "zf", bootstrap = FALSE)
  expect_equal(rownames(out$est_fs), c("zf_a.1", "zf_b", "zf_d"))
  expect_same_iv(out, run_iv(df, dums(), bootstrap = FALSE))
})

test_that("numeric instruments are untouched (no copy, no message)", {
  df <- sim_fac()
  expect_identical(ivDiag:::expand_factor_iv(df, dums()), list(data = df, Z = dums()))
  expect_silent(eff_F(df, "y", "d", dums()))
})
