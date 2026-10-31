# eff_F(): effective F = pi' Q pi / tr(Sigma Q), Q = Z~' W Z~, where Z~ are the
# instruments with the intercept, controls and fixed effects partialled out.
# The reference below is computed from scratch with base R matrix algebra.

sim_effF <- function(n = 300, seed = 1, zmean = c(0, 0, 0)) {
  set.seed(seed)
  x1 <- rnorm(n)
  x2 <- rnorm(n) + 0.5 * x1
  g <- sample(1:8, n, replace = TRUE)
  ge <- rnorm(8)[g]
  cl <- sample(1:25, n, replace = TRUE)
  # instruments correlated with the controls and the fixed effects
  z1 <- zmean[1] + 0.9 * x1 + 0.2 * x2 + ge + rnorm(n)
  z2 <- zmean[2] - 0.7 * x1 + 0.8 * x2 - 0.8 * ge + rnorm(n)
  z3 <- zmean[3] + 0.5 * x1 - 0.6 * x2 + 0.6 * ge + rnorm(n)
  u <- rnorm(n) * (1 + 0.5 * abs(x1))
  d <- 0.5 * z1 + 0.3 * z2 + 0.2 * z3 + 0.6 * x1 + 0.3 * x2 + ge + u
  y <- d + 0.5 * x1 + ge + rnorm(n) + 0.5 * u
  w <- exp(0.7 * rnorm(n))
  data.frame(y = y, d = d, z1 = z1, z2 = z2, z3 = z3, x1 = x1, x2 = x2,
             g = factor(g), cl = cl, w = w / mean(w))
}

# from-scratch effective F; HC1 / CR1 variance as in the package
effF_ref <- function(df, Z, controls = NULL, FE = NULL, cl = NULL, weights = NULL) {
  used <- c("y", "d", Z, controls, FE, cl, weights)
  df <- df[complete.cases(df[, used]), ]
  w <- if (is.null(weights)) rep(1, nrow(df)) else df[[weights]]
  rhs <- c(controls, if (!is.null(FE)) paste0("factor(", FE, ")"))
  X0 <- model.matrix(reformulate(if (length(rhs)) rhs else "1"), df)  # with intercept
  Zm <- as.matrix(df[, Z, drop = FALSE])
  # weighted residuals of the instruments on [1, controls, FE]
  Zt <- Zm - X0 %*% solve(crossprod(X0, X0 * w), crossprod(X0, Zm * w))
  Q <- crossprod(Zt, Zt * w)
  X <- cbind(Zm, X0)
  bread <- solve(crossprod(X, X * w))
  b <- bread %*% crossprod(X, df$d * w)
  e <- c(df$d - X %*% b)
  n <- nrow(X); k <- ncol(X)
  S <- X * (w * e)                       # scores
  if (!is.null(cl)) {
    S <- rowsum(S, df[[cl]])
    G <- nrow(S)
    adj <- G / (G - 1) * (n - 1) / (n - k)
  } else {
    adj <- n / (n - k)
  }
  V <- adj * bread %*% crossprod(S) %*% bread
  iv <- seq_along(Z)
  pi <- b[iv, , drop = FALSE]
  c(t(pi) %*% Q %*% pi / sum(diag(V[iv, iv, drop = FALSE] %*% Q)))
}

check_effF <- function(df, Z, controls = NULL, FE = NULL, cl = NULL, weights = NULL) {
  got <- eff_F(df, Y = "y", D = "d", Z = Z, controls = controls, FE = FE,
               cl = cl, weights = weights, prec = 12)
  expect_equal(got, effF_ref(df, Z, controls, FE, cl, weights), tolerance = 1e-8)
}

test_that("one instrument with controls is unchanged (Q cancels)", {
  df <- sim_effF()
  # old value: Q = ivx'ivx, which is any positive scalar for one instrument
  fs <- lfe::felm(y ~ x1 + x2 | 0 | (d ~ z1), data = df)$stage1
  pi1 <- fs$coefficients["z1", ]
  old <- unname(pi1^2 * c(crossprod(fs$ivx)) / (fs$robustvcv["z1", "z1"] * c(crossprod(fs$ivx))))
  expect_equal(eff_F(df, "y", "d", "z1", controls = c("x1", "x2"), prec = 12), old,
               tolerance = 1e-10)
  check_effF(df, "z1", controls = c("x1", "x2"))
})

test_that("two instruments with correlated controls, no fixed effects", {
  check_effF(sim_effF(), c("z1", "z2"), controls = c("x1", "x2"))
})

test_that("instruments with nonzero means and no controls or fixed effects", {
  df <- sim_effF(zmean = c(3, -2, 0))
  check_effF(df, c("z1", "z2"))
  # the intercept is partialled out, so shifting an instrument changes nothing
  df2 <- df
  df2$z1 <- df2$z1 + 10
  expect_equal(eff_F(df2, "y", "d", c("z1", "z2"), prec = 10),
               eff_F(df, "y", "d", c("z1", "z2"), prec = 10), tolerance = 1e-8)
})

test_that("fixed effects, with and without controls", {
  df <- sim_effF()
  check_effF(df, c("z1", "z2"), controls = c("x1", "x2"), FE = "g")
  check_effF(df, c("z1", "z2"), FE = "g")
})

test_that("weights, with and without fixed effects", {
  df <- sim_effF()
  check_effF(df, c("z1", "z2"), controls = c("x1", "x2"), weights = "w")
  check_effF(df, c("z1", "z2"), controls = c("x1", "x2"), FE = "g", weights = "w")
})

test_that("clustered standard errors, three instruments", {
  df <- sim_effF()
  check_effF(df, c("z1", "z2", "z3"), controls = c("x1", "x2"), cl = "cl")
  check_effF(df, c("z1", "z2", "z3"), controls = c("x1", "x2"), FE = "g", cl = "cl",
             weights = "w")
})

test_that("rows with missing values are dropped as in the estimation sample", {
  df <- sim_effF()
  df$x1[c(3, 40)] <- NA
  df$y[c(7, 90)] <- NA
  df$d[15] <- NA
  df$z2[22] <- NA
  df$cl[31] <- NA
  check_effF(df, c("z1", "z2"), controls = c("x1", "x2"), cl = "cl")
  check_effF(df, c("z1", "z2"), controls = c("x1", "x2"), FE = "g", weights = "w")
})

# ---- unidentified instruments: a clear error instead of NaN / NA / a wrong number
# (lfe only warns about rank deficiency, and may alias a control instead of the
# instrument, so these cases need an explicit check)

test_that("collinear instruments stop with a clear error", {
  df <- sim_effF()
  df$z1b <- 2 * df$z1
  expect_error(suppressWarnings(eff_F(df, "y", "d", c("z1", "z1b"), controls = c("x1", "x2"))),
               "instruments are collinear")
  expect_error(suppressWarnings(eff_F(df, "y", "d", c("z1", "z1b"))),
               "instruments are collinear")
  # numerically almost exact
  df$z1c <- df$z1 + 1e-10 * rnorm(nrow(df))
  expect_error(suppressWarnings(eff_F(df, "y", "d", c("z1", "z1c"), FE = "g")),
               "instruments are collinear")
})

test_that("a constant instrument stops with a clear error", {
  df <- sim_effF()
  df$zc <- 5
  expect_error(suppressWarnings(eff_F(df, "y", "d", c("z1", "zc"), controls = c("x1", "x2"))),
               "instruments are collinear")
  expect_error(suppressWarnings(eff_F(df, "y", "d", "zc")), "instruments are collinear")
  expect_error(suppressWarnings(eff_F(df, "y", "d", c("z1", "zc"), FE = "g")),
               "instruments are collinear")
})

test_that("an instrument collinear with the fixed effects stops with a clear error", {
  df <- sim_effF()
  df$zfe <- ave(df$z1, df$g)   # constant within each fixed-effect group
  expect_error(suppressWarnings(eff_F(df, "y", "d", c("z1", "zfe"), FE = "g")),
               "instruments are collinear")
  expect_error(suppressWarnings(eff_F(df, "y", "d", "zfe", controls = "x1", FE = "g")),
               "instruments are collinear")
})

test_that("an instrument in the span of the controls stops with a clear error", {
  df <- sim_effF()
  df$zx <- df$x1             # identical to a control
  df$z4 <- df$x1 + df$x2     # linear combination of the controls
  expect_error(suppressWarnings(eff_F(df, "y", "d", c("z1", "zx"), controls = c("x1", "x2"))),
               "instruments are collinear")
  expect_error(suppressWarnings(eff_F(df, "y", "d", "zx", controls = c("x1", "x2"))),
               "instruments are collinear")
  expect_error(suppressWarnings(eff_F(df, "y", "d", c("z1", "z4"), controls = c("x1", "x2"))),
               "instruments are collinear")
  expect_error(suppressWarnings(eff_F(df, "y", "d", c("z1", "z4"), controls = c("x1", "x2"),
                                      weights = "w", cl = "cl")),
               "instruments are collinear")
})

test_that("collinear controls are harmless", {
  df <- sim_effF()
  df$x3 <- 2 * df$x1
  df$x4 <- df$x1 - df$x2
  want <- eff_F(df, "y", "d", c("z1", "z2"), controls = c("x1", "x2"), prec = 12)
  got <- suppressWarnings(eff_F(df, "y", "d", c("z1", "z2"),
                                controls = c("x1", "x2", "x3"), prec = 12))
  expect_equal(got, want, tolerance = 1e-8)
  got2 <- suppressWarnings(eff_F(df, "y", "d", c("z1", "z2"),
                                 controls = c("x1", "x2", "x3", "x4"), prec = 12))
  expect_equal(got2, want, tolerance = 1e-8)
  expect_equal(got, effF_ref(df, c("z1", "z2"), controls = c("x1", "x2")), tolerance = 1e-8)
  # with fixed effects, weights and clusters; a constant control is all zero within FE
  df$kc <- 1
  want <- eff_F(df, "y", "d", c("z1", "z2"), controls = c("x1", "x2"), FE = "g",
                weights = "w", cl = "cl", prec = 12)
  got <- suppressWarnings(eff_F(df, "y", "d", c("z1", "z2"),
                                controls = c("x1", "x2", "x3", "kc"), FE = "g",
                                weights = "w", cl = "cl", prec = 12))
  expect_equal(got, want, tolerance = 1e-8)
})

test_that("weak or highly correlated but non-degenerate instruments still work", {
  df <- sim_effF()
  # correlated with z1 at about 0.99999, but not collinear
  df$z1n <- df$z1 + 0.005 * rnorm(nrow(df))
  check_effF(df, c("z1", "z1n"), controls = c("x1", "x2"))
  # an instrument with no first-stage effect gives a number, not an error
  df$zw <- rnorm(nrow(df))
  f <- eff_F(df, "y", "d", c("z1", "zw"), controls = c("x1", "x2"), prec = 12)
  expect_true(is.finite(f) && f > 0)
  check_effF(df, c("z1", "zw"), controls = c("x1", "x2"))
  # a large mean with a small spread is not degenerate
  df$zbig <- 1000 + df$z2
  check_effF(df, c("z1", "zbig"), controls = c("x1", "x2"))
})

# ---- instruments on very different scales or with a huge mean
# The effective F is invariant to nonsingular linear transformations of the
# instruments, so the reference is computed on standardised instruments (a
# base-R computation on the raw columns itself loses precision at these scales).

effF_ref_std <- function(df, Z, ...) {
  for (z in Z) df[[z]] <- (df[[z]] - mean(df[[z]])) / sd(df[[z]])
  effF_ref(df, Z, ...)
}

check_scaled <- function(df0, df, Z, ..., label = "") {
  got <- suppressWarnings(eff_F(df, "y", "d", Z, ..., prec = 12))
  expect_equal(got, effF_ref_std(df, Z, ...), tolerance = 1e-8, label = label)
  # and the same as with the original scales
  expect_equal(got, eff_F(df0, "y", "d", Z, ..., prec = 12), tolerance = 1e-8, label = label)
}

test_that("instruments on very different scales give the same value", {
  df0 <- sim_effF()
  df <- df0
  df$z1 <- df$z1 * 1e6
  df$z2 <- df$z2 * 1e-6
  check_scaled(df0, df, c("z1", "z2"), controls = c("x1", "x2"))
  check_scaled(df0, df, c("z1", "z2"))
  check_scaled(df0, df, c("z1", "z2"), controls = c("x1", "x2"), FE = "g", cl = "cl", weights = "w")
  df$z3 <- df$z3 * 1e-6
  check_scaled(df0, df, c("z1", "z2", "z3"), controls = c("x1", "x2"))
})

test_that("both instruments on a tiny scale, and on a huge scale, work", {
  df0 <- sim_effF()
  df <- df0
  df$z1 <- df$z1 * 1e-6
  df$z2 <- df$z2 * 1e-6
  check_scaled(df0, df, c("z1", "z2"), controls = c("x1", "x2"))
  df$z1 <- df0$z1 * 1e6
  df$z2 <- df0$z2 * 1e6
  check_scaled(df0, df, c("z1", "z2"), controls = c("x1", "x2"))
})

test_that("instruments with a huge mean work", {
  df0 <- sim_effF()
  df <- df0
  df$z1 <- df$z1 + 1e6
  df$z2 <- df$z2 + 1e6
  check_scaled(df0, df, c("z1", "z2"), controls = c("x1", "x2"))
  check_scaled(df0, df, c("z1", "z2"), FE = "g", weights = "w")
})

test_that("a missing value in an instrument does not break the standardisation", {
  df0 <- sim_effF()
  df0$z1[c(4, 50)] <- NA
  df <- df0
  df$z1 <- df$z1 * 1e6
  got <- eff_F(df, "y", "d", c("z1", "z2"), controls = c("x1", "x2"), prec = 12)
  expect_equal(got, eff_F(df0, "y", "d", c("z1", "z2"), controls = c("x1", "x2"), prec = 12),
               tolerance = 1e-8)
  expect_true(is.finite(got))
})

test_that("a constant instrument gives the collinearity error, never Inf or NaN", {
  df <- sim_effF()
  for (k in c(0, 5, 1e9, 0.1 + 0.2)) {
    df$zk <- k
    expect_error(suppressWarnings(eff_F(df, "y", "d", "zk", controls = "x1")),
                 "instruments are collinear")
    expect_error(suppressWarnings(eff_F(df, "y", "d", c("z1", "zk"), controls = "x1")),
                 "instruments are collinear")
  }
})

test_that("the identification tolerance behaves on standardised instruments", {
  df <- sim_effF()
  set.seed(5)
  e <- rnorm(nrow(df))
  # small but real independent variation is not collinear ...
  df$z1n <- df$z1 + 1e-3 * e
  expect_true(is.finite(eff_F(df, "y", "d", c("z1", "z1n"), controls = c("x1", "x2"))))
  df$zxn <- df$x1 + 1e-3 * e
  expect_true(is.finite(eff_F(df, "y", "d", c("z1", "zxn"), controls = c("x1", "x2"))))
  # ... and variation at rounding level is
  df$z1m <- df$z1 + 1e-10 * e
  expect_error(suppressWarnings(eff_F(df, "y", "d", c("z1", "z1m"), controls = c("x1", "x2"))),
               "instruments are collinear")
  df$zxm <- df$x1 + 1e-10 * e
  expect_error(suppressWarnings(eff_F(df, "y", "d", c("z1", "zxm"), controls = c("x1", "x2"))),
               "instruments are collinear")
})

# ---- ivDiag() gives the same clear error before any fit or inversion fails

test_that("ivDiag() stops with the clear error for unidentified instruments", {
  df <- sim_effF()
  df$zdup <- df$z1
  df$zmul <- -3.7 * df$z1
  df$zctl <- df$x1
  df$zspan <- 2 * df$x1 - 0.5 * df$x2
  df$zfe <- ave(df$z1, df$g)
  df$zk <- 5
  cases <- list(
    list(c("z1", "zdup"), NULL, NULL),
    list(c("z1", "zmul"), NULL, NULL),
    list(c("z1", "zctl"), c("x1", "x2"), NULL),
    list("zctl", c("x1", "x2"), NULL),
    list(c("z1", "zspan"), c("x1", "x2"), NULL),
    list("zspan", c("x1", "x2"), NULL),
    list(c("z1", "zfe"), "x1", "g"),
    list("zfe", "x1", "g"),
    list(c("z1", "zk"), "x1", NULL),
    list("zk", NULL, NULL)
  )
  for (cs in cases) {
    for (run_ar in c(FALSE, TRUE)) {
      expect_error(
        suppressWarnings(suppressMessages(ivDiag(df, "y", "d", cs[[1]], controls = cs[[2]],
          FE = cs[[3]], bootstrap = FALSE, run.AR = run_ar, parallel = FALSE))),
        "instruments are collinear",
        info = paste(cs[[1]], collapse = "+"))
    }
  }
  # with weights and clusters, and with the bootstrap switched on
  expect_error(
    suppressWarnings(suppressMessages(ivDiag(df, "y", "d", c("z1", "zdup"), weights = "w",
      cl = "cl", bootstrap = FALSE))),
    "instruments are collinear")
  expect_error(
    suppressWarnings(suppressMessages(ivDiag(df, "y", "d", c("z1", "zdup"), bootstrap = TRUE,
      nboots = 20, parallel = FALSE))),
    "instruments are collinear")
})

test_that("ivDiag() still works with collinear controls, and reports eff_F()", {
  df <- sim_effF()
  df$x3 <- 2 * df$x1
  out <- suppressWarnings(suppressMessages(
    ivDiag(df, "y", "d", c("z1", "z2"), controls = c("x1", "x2", "x3"),
           bootstrap = FALSE, run.AR = FALSE)))
  expect_equal(unname(out$F_stat["F.effective"]),
               eff_F(df, "y", "d", c("z1", "z2"), controls = c("x1", "x2")))
  out2 <- suppressMessages(ivDiag(df, "y", "d", c("z1", "z2"), controls = c("x1", "x2"),
    FE = "g", cl = "cl", weights = "w", bootstrap = FALSE, run.AR = FALSE))
  expect_equal(unname(out2$F_stat["F.effective"]),
               eff_F(df, "y", "d", c("z1", "z2"), controls = c("x1", "x2"), FE = "g",
                     cl = "cl", weights = "w"))
})

# ---- factor, character and logical instruments are passed to felm as before
# (a factor or character instrument is expanded into dummies; only numeric and
# logical instruments are standardised). Expected values were computed with
# eff_F() from commit ea04722, before the standardisation was added.

test_that("factor and character instruments give the values of ea04722", {
  set.seed(1)
  n <- 300
  x1 <- rnorm(n); x2 <- rnorm(n) + 0.5 * x1
  g <- sample(1:8, n, replace = TRUE); ge <- rnorm(8)[g]
  f <- sample(c("a", "b", "c"), n, replace = TRUE)
  z2 <- 0.7 * x1 + ge + rnorm(n)
  d <- 0.5 * z2 + 0.8 * (f == "b") - 0.6 * (f == "c") + 0.6 * x1 + ge + rnorm(n)
  y <- d + 0.5 * x1 + ge + rnorm(n)
  df <- data.frame(y = y, d = d, z2 = z2, x1 = x1, x2 = x2, g = factor(g),
                   zf = factor(f), zch = f, zl = z2 > 0, w = exp(0.5 * rnorm(n)),
                   cl = sample(1:25, n, replace = TRUE), stringsAsFactors = FALSE)
  e <- function(...) suppressWarnings(eff_F(df, "y", "d", ..., prec = 12))
  expect_equal(e(c("zf", "z2"), controls = c("x1", "x2")), 239.112982831, tolerance = 1e-8)
  expect_equal(e("zf", controls = c("x1", "x2")), 11.4875501166, tolerance = 1e-8)
  expect_equal(e(c("zch", "z2"), controls = c("x1", "x2")), 239.112982831, tolerance = 1e-8)
  expect_equal(e("zch"), 9.71963995531, tolerance = 1e-8)
  expect_equal(e(c("zf", "z2"), controls = "x1", FE = "g", weights = "w", cl = "cl"),
               47.5693144438, tolerance = 1e-8)
  expect_equal(e(c("zl", "z2"), controls = "x1"), 317.270741014, tolerance = 1e-8)
})

# ---- one instrument and two or more clustering variables: eff_F = F.cluster
# With multiway clusters felm repairs a non-positive-definite variance by
# zeroing its negative eigenvalues, which depends on the scale of the
# regressors, so eff_F() must not rescale the instrument there (1.0.8 did, and
# gave 55.1894 instead of 50.8149 for seed 29 below).

sim_2way <- function(seed, n = 300, G = 12, H = 8) {
  set.seed(seed)
  g <- sample(1:G, n, TRUE); h <- sample(1:H, n, TRUE)
  x1 <- rnorm(n) + rnorm(G)[g]; x2 <- rnorm(n) + rnorm(H)[h]
  z <- 3 + 0.5 * x1 + rnorm(G)[g] + rnorm(H)[h] + rnorm(n)
  u <- rnorm(n) + rnorm(G)[g] + rnorm(H)[h]
  d <- 0.4 * z + 0.5 * x1 + 0.3 * x2 + u + rnorm(n)
  y <- d + x1 + u + rnorm(n)
  data.frame(y, d, z, x1, x2, g, h)
}

test_that("one instrument, two-way clusters: eff_F equals the cluster-robust F", {
  for (s in c(13, 29, 35, 37)) {
    df <- sim_2way(s)
    out <- suppressWarnings(suppressMessages(ivDiag(df, "y", "d", "z", controls = c("x1", "x2"),
      cl = c("g", "h"), bootstrap = FALSE, run.AR = FALSE)))
    expect_equal(unname(out$F_stat["F.effective"]), unname(out$F_stat["F.cluster"]),
                 info = paste("seed", s))
    # and felm's own pi^2 / V on the raw instrument
    fs <- suppressWarnings(lfe::felm(y ~ x1 + x2 | 0 | (d ~ z) | g + h, data = df))$stage1
    ref <- unname(fs$coefficients["z", ]^2 / fs$clustervcv["z", "z"])
    expect_equal(suppressWarnings(eff_F(df, "y", "d", "z", controls = c("x1", "x2"),
                                        cl = c("g", "h"), prec = 12)),
                 ref, tolerance = 1e-8, info = paste("seed", s))
  }
  expect_equal(suppressWarnings(eff_F(sim_2way(29), "y", "d", "z", controls = c("x1", "x2"),
                                      cl = c("g", "h"))), 50.8149)
  # with fixed effects and weights
  df <- sim_2way(29)
  df$w <- exp(0.5 * rnorm(nrow(df)))
  out <- suppressWarnings(suppressMessages(ivDiag(df, "y", "d", "z", controls = "x2", FE = "g",
    cl = c("g", "h"), weights = "w", bootstrap = FALSE, run.AR = FALSE)))
  expect_equal(unname(out$F_stat["F.effective"]), unname(out$F_stat["F.cluster"]))
})
