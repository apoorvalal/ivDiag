#' Weak-IV Robust Anderson Rubin Test
# fork of `ivmodel::AR.test` that accepts FELM objects and fits models with FEs
#' @param data dataframe
#' @param Y outcome (string)
#' @param D treatment (string)
#' @param Z instrument (string)
#' @param controls control variables (character vector)
#' @param FE fixed effects (character vector)
#' @param cl Cluster name
#' @param weights weighting vector
#' @param prec precision of CI in string
#' @param CI whether to conduct inversion to get CI
#' @param alpha level of statitical significance
#' @param parallel no longer used; kept for backward compatibility
#' @param cores no longer used; kept for backward compatibility
#' @details The confidence set collects every value \eqn{\beta_0} that the AR
#'   test does not reject. The reduced-form coefficients on the instruments are
#'   linear in \eqn{\beta_0}, and their robust or clustered variance is
#'   quadratic in \eqn{\beta_0}. The function recovers both exactly from three
#'   regressions. With one instrument the set then solves a quadratic
#'   inequality: a bounded interval, two unbounded rays, the whole real line,
#'   or empty. With several instruments the AR statistic is evaluated on the
#'   whole real line (after the change of variable \eqn{\beta_0 = \tan\theta})
#'   and each boundary is refined by root finding.
#' @return A list with \code{Fstat} (the AR test of \eqn{\beta = 0}) and, when
#'   \code{CI = TRUE}, \code{ci} (the end points of the set, in increasing
#'   order, with \code{-Inf}/\code{Inf} for unbounded pieces), \code{ci.print},
#'   and \code{bounded}.
#' @importFrom lfe felm
#' @export

AR_test = function(
    data, Y, D, Z, controls, FE = NULL, cl = NULL,
    weights = NULL, prec = 4, CI = TRUE, alpha = 0.05,
    parallel = NULL, cores = NULL) {
  # keep rows with complete data
  X <- controls
  data <- data[, c(Y, D, Z, X, FE, cl, weights)]
  data <- data[complete.cases(data), ]
  p_iv <- length(Z)

  # residualise
  Ytil = suppressWarnings(partialer(Y, X = X, FE = FE, data = data, weights = weights))
  Dtil = suppressWarnings(partialer(D, X = X, FE = FE, data = data, weights = weights))
  Ztil = matrix(NA, nrow = nrow(data), ncol = p_iv)
  for (i in 1:length(Z)) {
    Ztil[, i] = suppressWarnings(partialer(Z[i], X = X, FE = FE, data = data, weights = weights))
  }
  d <- cbind.data.frame(Ytil, Dtil, Ztil, data[, c(cl, weights)])
  colnames(d) <- c(Y, D, Z, cl, weights)

  # AR test of beta = 0
  Fstat <- AR_one(0, d = d, Y = Y, D = D, Z = Z, cl = cl, weights = weights)

  # Confidence set
  if (CI == TRUE) {
    message("AR Test Inversion...\n")
    # centre and scale of the search: the 2SLS estimate and its SE
    ivfit <- IV(data = data, Y = Y, D = D, Z = Z, X = X, FE = FE, cl = cl, weights = weights)
    b0 <- ivfit$coef
    s <- ivfit$se
    if (!is.finite(b0)) b0 <- 0
    if (!is.finite(s) || s <= 0) s <- 1
    df2 <- Fstat[3]
    crit <- qf(1 - alpha, p_iv, df2)
    q <- AR_pieces(d = d, Y = Y, D = D, Z = Z, cl = cl, weights = weights, b0 = b0, s = s)
    if (q$exact == TRUE && p_iv == 1) {
      sets <- AR_set_quadratic(q, crit)
    } else if (q$exact == TRUE) {
      Ffun <- function(theta) AR_F_theta(theta, q)
      sets <- AR_set_theta(Ffun, crit, ngrid = 20001)
    } else {
      # the variance is not quadratic in beta (e.g. a multiway clustered
      # variance repaired to be positive definite): call the regression at
      # every point instead
      message("AR set from direct evaluation (variance is not quadratic in beta)\n")
      Ffun <- function(theta) {
        t <- tan(pmax(pmin(theta, pi / 2 - 1e-6), -pi / 2 + 1e-6))
        sapply(t, function(tt) {
          AR_one(b0 + s * tt, d = d, Y = Y, D = D, Z = Z, cl = cl, weights = weights)[1]
        })
      }
      sets <- AR_set_theta(Ffun, crit, ngrid = 801)
    }
    # from the scaled parameter t back to beta = b0 + s * t
    sets$ends <- b0 + s * sets$ends
    out.set <- AR_format(sets, prec)
    ci <- out.set$ci
    ci.print <- out.set$ci.print
    bounded <- out.set$bounded
  }

  # output
  if (CI == TRUE) {
    out <- list(Fstat = round(Fstat, prec), ci.print = ci.print, ci = ci, bounded = bounded)
  } else {
    out <- list(Fstat = round(Fstat, prec))
  }
  return(out)
}

# AR test of a single value beta on the residualised data
AR_one <- function(beta, d, Y, D, Z, cl, weights) {
  d[, Y] <- d[, Y] - beta * d[, D]
  p_iv <- length(Z)
  reg <- OLS(data = d, Y = Y, D = Z, X = NULL, FE = NULL, cl = cl, weights = weights)
  coef <- reg$coef
  vcov <- reg$vcv
  df2 <- reg$df
  df1 <- p_iv
  Fstat <- c((t(coef) %*% solve(vcov) %*% coef) / p_iv)
  pval <- pf(Fstat, df1, df2, lower.tail = FALSE)
  output <- c(Fstat, df1, df2, pval)
  names(output) <- c("F", "df1", "df2", "p")
  return(output)
}

# Write beta = b0 + s * t. The reduced-form coefficients are linear in t,
# coef(t) = a0 + a1 t, and their variance (robust, clustered, weighted) is
# quadratic, V(t) = v0 + v1 t + v2 t^2. Recover both from the regressions at
# t = -1, 0, 1, then check them at t = 2.5.
AR_pieces <- function(d, Y, D, Z, cl, weights, b0, s) {
  fit <- function(t) {
    dd <- d
    dd[, Y] <- dd[, Y] - (b0 + s * t) * dd[, D]
    reg <- OLS(data = dd, Y = Y, D = Z, X = NULL, FE = NULL, cl = cl, weights = weights)
    list(coef = c(reg$coef), vcv = as.matrix(reg$vcv))
  }
  m <- fit(-1); z <- fit(0); p <- fit(1)
  a0 <- z$coef
  a1 <- (p$coef - m$coef) / 2
  v0 <- z$vcv
  v1 <- (p$vcv - m$vcv) / 2
  v2 <- (p$vcv + m$vcv) / 2 - z$vcv
  # check at a fourth point
  tc <- 2.5
  chk <- fit(tc)
  coef.hat <- a0 + a1 * tc
  vcv.hat <- v0 + v1 * tc + v2 * tc^2
  err.coef <- max(abs(coef.hat - chk$coef)) / max(abs(chk$coef), .Machine$double.eps)
  err.vcv <- max(abs(vcv.hat - chk$vcv)) / max(abs(chk$vcv), .Machine$double.eps)
  exact <- is.finite(err.coef) && is.finite(err.vcv) && err.coef < 1e-6 && err.vcv < 1e-6
  list(a0 = a0, a1 = a1, v0 = v0, v1 = v1, v2 = v2, exact = exact)
}

# AR statistic at beta = b0 + s * tan(theta), theta in [-pi/2, pi/2].
# Multiplying coef by cos(theta) and V by cos(theta)^2 leaves the statistic
# unchanged and keeps it finite at theta = +-pi/2 (beta = +-Inf).
AR_F_theta <- function(theta, q) {
  p_iv <- length(q$a0)
  sapply(theta, function(th) {
    cs <- cos(th); sn <- sin(th)
    a <- q$a0 * cs + q$a1 * sn
    V <- q$v0 * cs^2 + q$v1 * cs * sn + q$v2 * sn^2
    c(t(a) %*% solve(V, a)) / p_iv
  })
}

# One instrument: accept t when (a0 + a1 t)^2 <= crit (v0 + v1 t + v2 t^2),
# i.e. A t^2 + B t + C <= 0.
AR_set_quadratic <- function(q, crit) {
  A <- c(q$a1^2 - crit * q$v2)
  B <- c(2 * q$a0 * q$a1 - crit * q$v1)
  C <- c(q$a0^2 - crit * q$v0)
  scale <- max(abs(c(A, B, C)))
  if (abs(A) <= 1e-12 * scale) { # linear: B t + C <= 0
    if (abs(B) <= 1e-12 * scale) {
      if (C <= 0) return(list(ends = c(-Inf, Inf))) else return(list(ends = numeric(0)))
    }
    r <- -C / B
    if (B > 0) return(list(ends = c(-Inf, r))) else return(list(ends = c(r, Inf)))
  }
  disc <- B^2 - 4 * A * C
  if (disc < 0) { # no real root
    if (A > 0) return(list(ends = numeric(0))) else return(list(ends = c(-Inf, Inf)))
  }
  r <- sort((-B + c(-1, 1) * sqrt(disc)) / (2 * A))
  if (A > 0) {
    list(ends = r) # bounded interval
  } else {
    list(ends = c(-Inf, r[1], r[2], Inf)) # two rays
  }
}

# Several instruments (or the direct fallback): scan theta in [-pi/2, pi/2],
# where theta = +-pi/2 is beta = +-Inf, and refine each change of status.
AR_set_theta <- function(Ffun, crit, ngrid) {
  theta <- seq(-pi / 2, pi / 2, length.out = ngrid)
  g <- Ffun(theta) - crit
  acc <- g <= 0
  # the statistic has the same limit at +-Inf: use one value at both ends
  acc[ngrid] <- acc[1]
  if (all(acc)) return(list(ends = c(-Inf, Inf)))
  if (!any(acc)) return(list(ends = numeric(0)))
  gf <- function(th) Ffun(th) - crit
  roots <- c()
  for (i in which(acc[-1] != acc[-ngrid])) {
    r <- tryCatch(uniroot(gf, c(theta[i], theta[i + 1]), tol = 1e-12)$root,
      error = function(e) (theta[i] + theta[i + 1]) / 2)
    roots <- c(roots, r)
  }
  ends <- tan(roots)
  # rays when the status at +-Inf is "accept"
  if (acc[1] == TRUE) ends <- c(-Inf, ends, Inf)
  list(ends = ends)
}

# Turn end points into the output (ci, ci.print, bounded) of earlier versions
AR_format <- function(sets, prec) {
  ends <- sets$ends
  fmt <- function(x) {
    if (x == -Inf) return("-Inf")
    if (x == Inf) return("Inf")
    sprintf(paste0("%.", prec, "f"), x)
  }
  bounded <- FALSE
  if (length(ends) == 0) {
    return(list(ci = NA, ci.print = "empty", bounded = FALSE))
  }
  if (length(ends) == 2 && ends[1] == -Inf && ends[2] == Inf) {
    return(list(ci = c(-Inf, Inf), ci.print = "(-Inf, Inf)", bounded = FALSE))
  }
  pieces <- c()
  for (k in seq(1, length(ends), by = 2)) {
    lo <- ends[k]; hi <- ends[k + 1]
    pieces <- c(pieces, paste0(ifelse(lo == -Inf, "(", "["), fmt(lo), ", ",
      fmt(hi), ifelse(hi == Inf, ")", "]")))
  }
  ci.print <- paste(pieces, collapse = " Union ")
  ci <- ends
  ci[is.finite(ci)] <- round(ci[is.finite(ci)], prec)
  bounded <- all(is.finite(ends))
  list(ci = ci, ci.print = ci.print, bounded = bounded)
}
