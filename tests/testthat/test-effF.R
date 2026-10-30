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
