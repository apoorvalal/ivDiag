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
#' @param alpha level of statitical significance
#' @importFrom lfe felm
#' @export

AR_test = function(data, Y, D, Z, controls, FE = NULL, cl = NULL,
   weights = NULL, prec = 4, alpha = 0.05 
   ) {

  # keep rows with complete data
  X <- controls
  data <- data[, c(Y, D, Z, X, FE, weights)]
  data <- data[complete.cases(data), ]  
  p_iv <- length(Z)

  # Run IV regression once
  ivfit <- IV(data = data, Y = Y, D = D, Z = Z, X = X, FE = FE, cl = cl, weights = weights)
  beta_iv <- ivfit$coef
  se_iv <- ivfit$se
  beta_seq <- c(
    seq(beta_iv - 10*se_iv, beta_iv -3.1*se_iv, length.out = 100),
    seq(beta_iv - 3*se_iv, beta_iv + 3*se_iv, by = 0.02*se_iv),
    seq(beta_iv + 3.1*se_iv, beta_iv + 10*se_iv, length.out = 100)
  )
  ngrid <- length(beta_seq)
   
  # residualise
  Ytil = suppressWarnings(partialer(Y, X = X, FE = FE, data = data, weights = weights))
  Dtil = suppressWarnings(partialer(D, X = X, FE = FE, data = data, weights = weights))
  Ztil = matrix(NA, nrow = nrow(data), ncol = p_iv)
  for (i in 1:length(Z)) {
    Ztil[, i] = suppressWarnings(partialer(Z[i], X = X, FE = FE, data = data, weights = weights))
  }
  if (p_iv != 1) {Ztil = qr(Ztil)[[1]]}

  d <- cbind.data.frame(Ytil, Dtil, Ztil, data[, c(cl, weights)])
  colnames(d) <- c(Y, D, Z, cl, weights)

  # reduced form
  fmla = formula_lfe(Y = Y, W = Z, X = NULL, D = NULL, C = cl)
  if (is.null(weights) == TRUE) {
    m1 = lfe::felm(fmla, data = d)
  } else {
    m1 = lfe::felm(fmla, data = d, weights = d[, weights])
  }
  m1 <- robustify(m1)
  s <- summary(m1)
  Fstat <- s$F.fstat

  # Confidence intervals
  accept <- rep(NA, ngrid)
  for (i in 1:ngrid) {
    d[, Y] <- Ytil - beta_seq[i] * Dtil
    if (is.null(weights) == TRUE) {
      m2 = lfe::felm(fmla, data = d)
    } else {
      m2 = lfe::felm(fmla, data = d, weights = d[, weights])
    }
    m2 <- robustify(m2)
    s2 <- summary(m2)
    accept[i] <- ifelse(s2$pval >= alpha, 1, 0)    
  }

  # summarize
  bounded <- FALSE
  if (sum(accept)==ngrid) {
    ci <- "(-Inf, Inf)" # all accepted
  } else if (sum(accept)==0) {
    ci <- "empty"
  } else if (accept[1]==0 && accept[ngrid]==0) {
    betas <- range(beta_seq[accept == 1])
    ci <- paste0("[", sprintf(paste0("%.",prec,"f"), betas[1]), ", ", sprintf(paste0("%.",prec,"f"), betas[2]),"]") # bounded interval
    bounded <- TRUE
  } else if (accept[1]==1 && accept[ngrid]==1) {
    betas <- range(beta_seq[accept == 0]) # e.g. 1 1 1 1 0 0 0 1 1 1 
    ci <- paste0("(-Inf, ", sprintf(paste0("%.",prec,"f"), betas[1]), "] Union [", sprintf(paste0("%.",prec,"f"), betas[2]),", Inf)") 
  } else if (accept[1]==0 && accept[ngrid]==1) {
    betas <- range(beta_seq[accept == 1]) # e.g. 0 0 0 1 1 1 1 1
    ci <- paste0("[", sprintf(paste0("%.",prec,"f"), betas[1]), ", Inf)") 
  } else if (accept[1]==1 && accept[ngrid]==0) {
    betas <- range(beta_seq[accept == 1]) # e.g. 1 1 1 1 1 0 0 0 
    ci <- paste0("(-Inf, ", sprintf(paste0("%.",prec,"f"), betas[2]), "]") 
  }  

  return(list(Fstat = round(Fstat, prec),  ci = ci, bounded = bounded))
}
