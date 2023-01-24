#' Stitches together formula for use in felm
#' @param y The dependent variable
#' @param X vector of controls
#' @param W treatment variable
#' @param D vector of factor variables to be partialed out
#' @param Z vector of instruments
#' @param C vector of variables cluster standard errors (multi-way permitted by LFE)
#' @export
formula_lfe = function(Y, X, W = NULL, D = NULL, Z = NULL, C = NULL) {
  # 'second stage' step
  if (!is.null(W) & is.null(Z)) { # separate treatment dummy only
    felm_ss = paste(c(Y, paste(c(W, X), collapse = "+")), collapse = "~")
  } else { # no instrumented variable
    if (!is.null(X)) {
      felm_ss = paste(c(Y, paste(X, collapse = "+")), collapse = "~")
    } else {
      felm_ss = paste(c(Y, 1), collapse = "~")
    }
  }
  # first stage
  if (!is.null(Z)) {
    felm_fs = paste(c("(", paste(c(W, paste(Z, collapse = "+")),
      collapse = "~"
    ), ")"), collapse = "")
  } else { # no instrument
    felm_fs = "0"
  }
  # FEs
  if (!is.null(D)) {
    facs = paste(D, collapse = "+")
  } else {
    facs = "0"
  }
  # clusters
  if (!is.null(C)) {
    clusts = paste(C, collapse = "+")
  } else {
    clusts = "0"
  }
  # return formula
  as.formula(paste(c(felm_ss, facs, felm_fs, clusts), collapse = "|"))
}

# %%
#' Stitches together formula for use in fixest
#' @param y The dependent variable
#' @param X vector of controls
#' @param W treatment variable
#' @param D vector of factor variables to be partialed out
#' @param Z vector of instruments
formula_fixest = function(y, X, W = NULL, D = NULL, Z = NULL) {
  if (!is.null(W) & is.null(Z)) {
    fixest_ss = paste(c(y, paste(c(W, X), collapse = "+")),
      collapse = "~"
    )
  } else {
    fixest_ss = paste(c(y, paste(X, collapse = "+")), collapse = "~")
  }
  if (!is.null(D))
    facs = paste(D, collapse = "+")
  else facs = "0"
  if (!is.null(Z)) {
    fixest_fs = paste(c(paste(c(W, paste(Z, collapse = "+")),
      collapse = "~"
    )), collapse = "")
    as.formula(paste(c(fixest_ss, facs, fixest_fs), collapse = "|"))
  } else {
    as.formula(paste(c(fixest_ss, facs), collapse = "|"))
  }
}


# %% partialer
#' partial out controls and covariates from y
#' @param y The dependent variable
#' @param X vector of controls
#' @param FEs fixed effects
#' @param data Dataframe
#' @param weights name of the weighting variable
#' @export
partialer = function(Y, X, FE = NULL, data, weights = NULL) {
  # regress variable on controls and FEs and return residuals
  f = formula_lfe(Y = Y, X = X, D = FE)
  if (is.null(weights) == TRUE) {
    m <- lfe::felm(f, data)
  } else {
    m <- lfe::felm(f, data, weights = data[, weights])
  }
  return(m$residuals[, 1])
}

# %% run IV through FELM
IV <- function(data, D, Y, Z, X = NULL, FE = NULL, cl = NULL, weights = NULL # weights is a string
) {
  fmla = formula_lfe(Y = Y, W = D, Z = Z, X = X, D = FE, C = cl)
  if (is.null(weights)) {
    m2 = robustify(lfe::felm(fmla, data = data))
  } else {
    m2 = robustify(lfe::felm(fmla, data = data, weights = data[, weights]))
  }
  ## Output
  coef <- c(tail(m2$coefficients, n = 1)) # felm IV fits have endog coef at the tail
  se <- c(tail(sqrt(diag(m2$vcv)), n = 1))
  df <- m2$df.residual
  names(coef) <- names(se) <- names(df) <- NULL
  return(list(coef = coef, se = se, df = df))
}

# %% get robust and clustered F
get.vcov <- function(data, D, Y, Z, X = NULL, FE = NULL, cl = NULL, weights = NULL # weights is a string
) {
  p_iv <- length(Z)
  fmla = formula_lfe(Y = Y, W = D, Z = Z, X = X, D = FE, C = cl)
  if (is.null(weights)) {
    m2 = robustify(lfe::felm(fmla, data = data))
  } else {
    m2 = robustify(lfe::felm(fmla, data = data, weights = data[, weights]))
  }
  stage1 <- m2$stage1
  p <- nrow(stage1$vcv)
  iv.pos <- (p - p_iv + 1):p
  vcov.standard <- stage1$vcv[iv.pos, iv.pos]
  vcov.robust <- stage1$robustvcv[iv.pos, iv.pos]
  if (is.null(cl) == FALSE) {
    vcov.cluster <- stage1$clustervcv[iv.pos, iv.pos]
  } else {
    vcov.cluster <- NA
  }
  ## Output
  out <- c(list(
    vcov.standard = vcov.standard,
    vcov.robust = vcov.robust,
    vcov.cluster = vcov.cluster
  ))
  return(out)
}


# %% wrapper around felm to run OLS with and without controls
OLS <- function(data, Y, D, X = NULL, FE = NULL, cl = NULL, weights = NULL # weights is a string
) {
  p_D <- length(D)
  fmla = formula_lfe(Y = Y, W = D, X = X, D = FE, C = cl)
  if (is.null(weights) == TRUE) {
    m1 = robustify(lfe::felm(fmla, data = data))
  } else {
    m1 = robustify(lfe::felm(fmla, data = data, weights = data[, weights]))
  }
  if (is.null(FE) == TRUE) {
    coef <- m1$coefficients[2:(1 + p_D)]
    se <- sqrt(diag(m1$vcv)[2:(1 + p_D)])
  } else {
    coef = m1$coefficients[1:p_D]
    se <- sqrt(diag(m1$vcv)[1:p_D])
  }
  df <- m1$df.residual
  out <- list(coef = coef, se = se, df = df)
  return(out)
}

# %%
first_stage_coefs <- function(data, D, Z, X, FE = NULL, weights = NULL # weights is a string
) {
  p_iv <- length(Z)
  formula <- formula_lfe(Y = D, W = Z, X = X, D = FE)
  if (is.null(weights)) {
    reg = robustify(lfe::felm(formula, data = data))
  } else {
    reg = robustify(lfe::felm(formula, data = data, weights = data[, weights]))
  }
  # slice model fit
  if (is.null(FE) == FALSE) {
    coefs <- reg$coefficients[1:p_iv]
  } else {
    coefs <- reg$coefficients[2:(1 + p_iv)]
  }
  return(coefs)
}

# %% first stage correlation coefficient for each IV
first_stage_rho = function(data, D, Z, X, FE = NULL, weights = NULL # weights is a string
) {
  p_iv <- length(Z)
  rho <- rep(NA, p_iv)
  if (p_iv > 1) {
    names(rho) <- Z
  }
  res.d <- partialer(Y = D, X = X, FE = FE, data = data, weights = weights)
  for (i in 1:p_iv) {
    res.z <- partialer(Y = Z[i], X = X, FE = FE, data = data, weights = weights)
    if (is.null(weights) == TRUE) {
      rho[i] <- cor(res.d, res.z)
    } else {
      rho[i] <- wCorr::weightedCorr(x = res.d, y = res.z, weights = data[, weights], method = "Pearson")
    }
  }
  return(rho)
}

# %% # function to replace SE and t-stats of FELM models inplace
#' @param model FELM fit object
robustify = function(model) {
  model$se = model$rse
  model$tval = model$rtval
  model$pval = model$rpval
  return(model)
}

# LMMP's tf procedure
#' @param coef 2SLS coefficient
#' @param se SE for 2SLS estimate
#' @param Fstat first stage F
tF <- function(coef, se, Fstat) {
  tstat <- coef / se
  F.sqrt <- sqrt(Fstat)
  F0.sqrt <- seq(2, 10.3, 0.1)
  cF0 <- c(
    18.66, 9.74, 7.37, 6.18, 5.43, 4.92, 4.54, 4.25, 4.01, 3.82, 3.65, 3.51, 3.39, 3.29, 3.19, 3.11, 3.03,
    2.97, 2.91, 2.85, 2.80, 2.75, 2.71, 2.67, 2.63, 2.60, 2.57, 2.54, 2.51, 2.48, 2.46, 2.43, 2.41, 2.39, 2.37,
    2.35, 2.33, 2.32, 2.30, 2.29, 2.27, 2.26, 2.24, 2.23, 2.22, 2.21, 2.20, 2.19, 2.17, 2.16, 2.16, 2.15, 2.14,
    2.13, 2.12, 2.11, 2.10, 2.10, 2.09, 2.08, 2.08, 2.07, 2.06, 2.06, 2.05, 2.04, 2.04, 2.03, 2.03, 2.02, 2.02,
    2.01, 2.01, 2.00, 2.00, 1.99, 1.99, 1.99, 1.98, 1.98, 1.97, 1.97, 1.97, 1.96
  )
  if (F.sqrt <= 2) {
    cF <- cF0[1]
  } else if (F.sqrt >= 10.3) {
    cF <- 1.96
  } else {
    pos.lower <- max(which(F0.sqrt < F.sqrt))
    pos.upper <- pos.lower + 1
    h1 <- abs(F.sqrt - F0.sqrt[pos.lower])
    h2 <- abs(F.sqrt - F0.sqrt[pos.upper])
    # critical value
    cF <- (cF0[pos.upper] * h1 + cF0[pos.lower] * h2) / (h1 + h2) # weighted average
  }
  ci <- c(coef - cF * se, coef + cF * se)
  p <- (1 - pnorm(abs(tstat) / (cF / 1.96))) * 2 # adjusted p value
  out <- c(Fstat, cF, coef, se, tstat, ci, p)
  names(out) <- c("F", "cF", "Coef", "SE", "t", "CI2.5%", "CI97.5%", "p-value")
  return(out)
}
