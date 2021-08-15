#' Stitches together formula for use in felm
#' @param y The dependent variable
#' @param X vector of controls
#' @param W treatment variable
#' @param D vector of factor variables to be partialed out
#' @param Z vector of instruments
#' @param C vector of variables cluster standard errors (multi-way permitted by LFE)
#' @export
#' @examples
#' formula_lfe(Y='mpg', X = c('hp', 'drat'), D = c('wt', 'vs'))
#' formula_lfe(Y='mpg', X = c('hp', 'drat'), W = 'gear', Z = c('cyl', 'carb'), D = c('wt', 'vs'), C = c('cyl', 'wt'))

formula_lfe = function (Y, X, W = NULL, D = NULL, Z = NULL, C = NULL) {
  # 'second stage' step
  if (!is.null(W) & is.null(Z)) { # separate treatment dummy only
    felm_ss = paste(c(Y, paste(c(W, X), collapse = "+")), collapse = "~")

    } else { # no instrumented variable
    if (!is.null(X)){
    felm_ss = paste(c(Y, paste(X, collapse = "+")), collapse = "~")
    } else {
      felm_ss = paste(c(Y, 1), collapse="~")
    }
  }
  # first stage
  if (!is.null(Z)) {
    felm_fs = paste(c("(", paste(c(W, paste(Z, collapse = "+")),
                                 collapse = "~"), ")"), collapse = "")
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

# %% partialer
#' partial out controls and covariates from y
#' @param y The dependent variable
#' @param X vector of controls
#' @param FEs fixed effects
#' @param data Dataframe
#' @param wegihts name of the weighting variable
partialer = function(Y, X, FE = NULL, data, weights = NULL){
  # regress variable on controls and FEs and return residuals
  f = formula_lfe(Y = Y, X = X, D = FE)
  if(is.null(weights) == TRUE){
    m <- lfe::felm(f, data)
  } else {
    m <- lfe::felm(f, data, weights = data[,weights])
  }
  return(m$residuals[, 1])
}

# %% wrapper around felm to run OLS with and without controls
OLS <- function(data, D, Y, X=NULL, FE=NULL, cl=NULL, weights=NULL # weights is a string
  ){
  fmla = formula_lfe(Y = Y, W = D, X = X, D = FE, C = cl)
  if (is.null(weights)==TRUE) {
    m1 = robustify(lfe::felm(fmla, data = data))
  } else{
    m1 = robustify(lfe::felm(fmla, data = data, weights = data[,weights]))
  }
  if(is.null(FE)==TRUE){
    coef <- m1$coefficients[2]
    se <- sqrt(diag(m1$vcv)[2])
  } else{
    coef = m1$coefficients[1]
    se <- sqrt(diag(m1$vcv)[1])
  }
  out <- list(coef = coef, se = se)
  return(out)
}

# %% run IV through FELM
IV <- function(data, D, Y, Z, X=NULL, FE=NULL, cl=NULL, weights=NULL # weights is a string
  ){
  fmla = formula_lfe(Y = Y, W = D, Z = Z, X = X, D = FE, C = cl)
  if(is.null(weights)){
    m2 = robustify(lfe::felm(fmla, data = data))
  }
  else{
    m2 = robustify(lfe::felm(fmla, data = data, weights = data[,weights]))
  }
  ##Output
  coef <- c(tail(m2$coefficients, n = 1)) # felm IV fits have endog coef at the tail
  se <- c(tail(sqrt(diag(m2$vcv)), n = 1))
  names(coef) <- names(se) <- NULL
  return(list(coef = coef, se = se))
}

# %% get robust and clustered F
get.vcov <- function(data, D, Y, Z, X=NULL, FE=NULL, cl=NULL, weights=NULL # weights is a string
  ){
  p_iv <- length(Z)
  fmla = formula_lfe(Y = Y, W = D, Z = Z, X = X, D = FE, C = cl)
  if(is.null(weights)){
    m2 = robustify(lfe::felm(fmla, data = data))
  }
  else{
    m2 = robustify(lfe::felm(fmla, data = data, weights = data[,weights]))
  }
  stage1 <- m2$stage1
  p <- nrow(stage1$vcv)
  iv.pos <- (p-p_iv+1):p
  vcov.standard <- stage1$vcv[iv.pos, iv.pos]
  vcov.robust <- stage1$robustvcv[iv.pos, iv.pos]
  if (is.null(cl)==FALSE) {
    vcov.cluster <- stage1$clustervcv[iv.pos, iv.pos]
  } else {
    vcov.cluster <- NA
  }
  ##Output
  out <- c(list(vcov.standard = vcov.standard,
    vcov.robust = vcov.robust,
    vcov.cluster = vcov.cluster))
  return(out)
}


# %%
first_stage_coefs <- function(data, D, Z, X, FE = NULL, weights = NULL # weights is a string
  ) {
  p_iv <- length(Z)
  formula <- formula_lfe(Y = D, W = Z, X = X, D = FE)
  if(is.null(weights)){
    reg = robustify(lfe::felm(formula, data = data))
  } else{
    reg = robustify(lfe::felm(formula, data = data, weights = data[,weights]))
  }
  # slice model fit
  if (is.null(FE) == FALSE) {
    coefs <- reg$coefficients[1:p_iv]
  } else {
    coefs <- reg$coefficients[2:(1+p_iv)]
  }
  return(coefs)
}

# %% first stage correlation coefficient for each IV
first_stage_rho = function(data, D, Z, X, FE = NULL, weights = NULL # weights is a string
  ) {
  p_iv <- length(Z)
  rho <- rep(NA, p_iv)
  if (p_iv > 1) {names(rho) <- Z}
  res.d <- partialer(Y = D, X = X, FE = FE, data = data, weights = weights)
  for (i in 1:p_iv) {
    res.z <- partialer(Y = Z[i], X = X, FE = FE, data = data, weights = weights)
    if (is.null(weights) == TRUE) {
      rho[i] <- cor(res.d, res.z)
    } else {
      rho[i] <- weightedCorr(x = res.d, y = res.z, weights = data[,weights], method = "Pearson")
    }
  }
  return(rho)
}

# %% # function to replace SE and t-stats of FELM models inplace
#' @param model FELM fit object
robustify = function (model) {
    model$se = model$rse
    model$tval = model$rtval
    model$pval = model$rpval
    return(model)
}
