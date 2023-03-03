#' Stitches together formula for use in felm
#' @param y The dependent variable
#' @param X vector of controls
#' @param W treatment variable
#' @param FE vector of factor variables to be partialed out
#' @param Z vector of instruments
#' @param Cl vector of variables cluster standard errors (multi-way permitted by LFE)
#' @export
formula_lfe <- function(Y, X, W = NULL, FE = NULL, Z = NULL, Cl = NULL) {
  # 'second stage' step
  if (!is.null(W) & is.null(Z)) { # there is W, but no instrument, such as OLS
    felm_ss = paste(c(Y, paste(c(W, X), collapse = "+")), collapse = "~")
  } else { # either there is Z or there is no W 
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
  # FEs (factorial variables)
  if (!is.null(FE)) {
    facs = paste(FE, collapse = "+")
  } else {
    facs = "0"
  }
  # clusters
  if (!is.null(Cl)) {
    clusts = paste(Cl, collapse = "+")
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
partialer <- function(Y, X, FE = NULL, data, weights = NULL) {
  # regress variable on controls and FEs and return residuals
  f = formula_lfe(Y = Y, X = X, FE = FE)
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
  fmla = formula_lfe(Y = Y, W = D, Z = Z, X = X, FE = FE, Cl = cl)
  if (is.null(weights)) {
    m2 = robustify(lfe::felm(fmla, data = data))
  } else {
    m2 = robustify(lfe::felm(fmla, data = data, weights = data[, weights]))
  }
  ## Output
  coef <- c(tail(m2$coefficients, n = 1)) # felm IV fits have endog coef at the tail
  se <- c(tail(sqrt(diag(m2$robustvcv)), n = 1))
  df <- m2$df.residual
  names(coef) <- names(se) <- names(df) <- NULL
  return(list(coef = coef, se = se, df = df))
}

# %% get robust and clustered F
get.vcov <- function(data, D, Y, Z, X = NULL, FE = NULL, cl = NULL, weights = NULL # weights is a string
) {
  p_iv <- length(Z)
  fmla = formula_lfe(Y = Y, W = D, Z = Z, X = X, FE = FE, Cl = cl)
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
  fmla = formula_lfe(Y = Y, W = D, X = X, FE = FE, Cl = cl)
  if (is.null(weights) == TRUE) {
    m1 = robustify(lfe::felm(fmla, data = data))
  } else {
    m1 = robustify(lfe::felm(fmla, data = data, weights = data[, weights]))
  }
  if (is.null(FE) == TRUE) {
    coef <- m1$coefficients[2:(1 + p_D)]
    se <- sqrt(diag(m1$robustvcv)[2:(1 + p_D)])
  } else {
    coef = m1$coefficients[1:p_D]
    se <- sqrt(diag(m1$robustvcv)[1:p_D])
  }
  df <- m1$df.residual
  out <- list(coef = coef, se = se, df = df)
  return(out)
}

# %%
first_stage_coefs <- function(data, D, Z, X, FE = NULL, weights = NULL # weights is a string
) {
  p_iv <- length(Z)
  formula <- formula_lfe(Y = D, W = Z, X = X, FE = FE)
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

# %% first stage correlation coefficient between D and predicted D
first_stage_rho = function(data, D, Z, X, FE = NULL, weights = NULL # weights is a string
) {
  p_iv <- length(Z)
  # partial out covariates
  res.d <- partialer(Y = D, X = X, FE = FE, data = data, weights = weights)
  res.z <- matrix(NA, nrow(data), p_iv)
  for (i in 1:p_iv) {
    res.z[,i] <- partialer(Y = Z[i], X = X, FE = FE, data = data, weights = weights)
  }
  d0 <- cbind.data.frame(res.d, res.z); colnames(d0) <- c(D, Z)
  # first stage
  fmla = formula_lfe(Y = D, W = Z, X = NULL, FE = NULL, Cl = NULL)
  if (is.null(weights) == TRUE) {
    m = lfe::felm(fmla, data = d0)
  } else {
    m = lfe::felm(fmla, data = d0, weights = data[, weights])
  } 
  # rho
  if (is.null(weights) == TRUE) {
    rho <- cor(res.d, m$fitted.values, method = c("pearson")) 
  } else {
    rho <- wCorr::weightedCorr(x = res.d, y = m$fitted.values, weights = data[, weights], method = "Pearson")
  }
  rho <- c(rho)
  names(rho) <- NULL   
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

