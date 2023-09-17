#' Stitches together formula for use in felm
#' @param y The dependent variable
#' @param X vector of controls
#' @param D treatment variable
#' @param FE vector of factor variables to be partialed out
#' @param Z vector of instruments
#' @param cl vector of variables cluster standard errors (multi-way permitted by LFE)
#' @export
formula_lfe <- function(Y, X, D = NULL, FE = NULL, Z = NULL, cl = NULL) {
  # 'second stage' step
  if (!is.null(D) & is.null(Z)) { # there is W, but no instrument, such as OLS
    felm_ss = paste(c(Y, paste(c(D, X), collapse = "+")), collapse = "~")
  } else { # either there is Z or there is no W
    if (!is.null(X)) {
      felm_ss = paste(c(Y, paste(X, collapse = "+")), collapse = "~")
    } else {
      felm_ss = paste(c(Y, 1), collapse = "~")
    }
  }
  # first stage
  if (!is.null(Z)) {
    felm_fs = paste(c("(", paste(c(D, paste(Z, collapse = "+")),
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
  if (!is.null(cl)) {
    clusts = paste(cl, collapse = "+")
  } else {
    clusts = "0"
  }
  # return formula
  as.formula(paste(c(felm_ss, facs, felm_fs, clusts), collapse = "|"))
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


# $$$$$$\ $$\    $$\ 
# \_$$  _|$$ |   $$ |
#   $$ |  $$ |   $$ |
#   $$ |  \$$\  $$  |
#   $$ |   \$$\$$  / 
#   $$ |    \$$$  /  
# $$$$$$\    \$  /   
# \______|    \_/    
                   

# %% run IV through FELM
IV <- function(data, D, Y, Z, X = NULL, FE = NULL, cl = NULL, weights = NULL # weights is a string
) {
  fmla = formula_lfe(Y = Y, D = D, Z = Z, X = X, FE = FE, cl = cl)
  if (is.null(weights)) {
    m2 = lfe::felm(fmla, data = data)
  } else {
    m2 = lfe::felm(fmla, data = data, weights = data[, weights])
  }
  ## Output
  coef <- c(tail(m2$coefficients, n = 1)) # felm IV fits have endog coef at the tail
  if (is.null(cl) == TRUE) {
    se <- tail(m2$rse, n = 1)
  } else {
    se <- tail(m2$cse, n = 1)
  }
  df <- m2$df.residual
  names(coef) <- names(se) <- names(df) <- NULL
  return(list(coef = coef, se = se, df = df))
}

# %% get robust and clustered F
get.vcov <- function(data, D, Y, Z, X = NULL, FE = NULL, cl = NULL, weights = NULL # weights is a string
) {
  p_iv <- length(Z)
  fmla = formula_lfe(Y = Y, D = D, Z = Z, X = X, FE = FE, cl = cl)
  if (is.null(weights)) {
    m2 = lfe::felm(fmla, data = data)
  } else {
    m2 = lfe::felm(fmla, data = data, weights = data[, weights])
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

#  $$$$$$\  $$\       $$$$$$\  
# $$  __$$\ $$ |     $$  __$$\ 
# $$ /  $$ |$$ |     $$ /  \__|
# $$ |  $$ |$$ |     \$$$$$$\  
# $$ |  $$ |$$ |      \____$$\ 
# $$ |  $$ |$$ |     $$\   $$ |
#  $$$$$$  |$$$$$$$$\\$$$$$$  |
#  \______/ \________|\______/ 
                             
 # %% wrapper around felm to run OLS with and without controls
OLS <- function(data, Y, D, X = NULL, FE = NULL, cl = NULL, weights = NULL # weights is a string
) {
  p_D <- length(D)
  fmla = formula_lfe(Y = Y, D = D, X = X, FE = FE, cl = cl)
  if (is.null(weights) == TRUE) {
    m1 = lfe::felm(fmla, data = data)
  } else {
    m1 = lfe::felm(fmla, data = data, weights = data[, weights])
  }
  if (is.null(cl) == TRUE) {
    SE <- m1$rse # robust SE
    VCV <- m1$robustvcv # robust VCV
  } else {
    SE <- m1$cse # clustered SE
    VCV <- m1$clustervcv # clustered VCV
  }
  if (is.null(FE) == TRUE) {
    coef <- m1$coefficients[2:(1 + p_D)]
    se <- SE[2:(1 + p_D)]
    vcv <- VCV[2:(1 + p_D), 2:(1 + p_D)]
  } else {
    coef = m1$coefficients[1:p_D]
    se <- SE[1:p_D]
    vcv <- VCV[1:p_D, 1:p_D]
  }
  df <- m1$df.residual
  out <- list(coef = coef, se = se, vcv = vcv, df = df)
  return(out)
}

# %%
first_stage_coefs <- function(data, D, Z, X, FE = NULL, weights = NULL # weights is a string
) {
  p_iv <- length(Z)
  formula <- formula_lfe(Y = D, D = Z, X = X, FE = FE)
  if (is.null(weights)) {
    reg = lfe::felm(formula, data = data)
  } else {
    reg = lfe::felm(formula, data = data, weights = data[, weights])
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
  data <- data[, unique(c(D, Z, X, FE, weights))]
  data <- data[complete.cases(data), ]
  p_iv <- length(Z)
  # partial out covariates
  res.d <- partialer(Y = D, X = X, FE = FE, data = data, weights = weights)
  res.z <- matrix(NA, nrow(data), p_iv)
  for (i in 1:p_iv) {
    res.z[, i] <- partialer(Y = Z[i], X = X, FE = FE, data = data, weights = weights)
  }
  d0 <- cbind.data.frame(res.d, res.z); colnames(d0) <- c(D, Z)
  # first stage
  fmla = formula_lfe(Y = D, D = Z, X = NULL, FE = NULL, cl = NULL)
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


