# %% Local to zero IV estimate with direct effect of instrument set to μ
#' Estimates Local-to-Zero IV coefficients and SEs for a single instrument
#' @param data dataframe
#' @param Y outcome (string)
#' @param D treatment (string)
#' @param Z instrument (string)
#' @param controls control variables (character vector)
#' @param FE fixed effects (character vector)
#' @param cl clustering column for SE (character vector)
#' @param mu prior point estimate of direct effect of instrument on outcome
#' @param Sig variance of above prior
#' @export
ltz = function(data, Y, D, Z, controls, FE = NULL, cl = NULL, weights = NULL, mu, Sig) {
  # requirements: 2SLS coef and se, residualized Z and D, and prior on μ, Sigma
  
  # drop missingness
  data <- data[, unique(c(Y, D, Z, controls, FE, cl, weights))]
  data <- haven::zap_labels(data)
  data <- as.data.frame(data[complete.cases(data), ])
  n <- nrow(data); p_iv <- length(Z)

  # IV fit
  fmla = formula_lfe(Y = Y, W = D, Z = Z, X = controls, FE = FE, Cl = cl)
  if (is.null(weights)) {
    m2 = robustify(lfe::felm(fmla, data = data))
  } else {
    m2 = robustify(lfe::felm(fmla, data = data, weights = data[, weights]))
  }
  beta_2sls <- c(tail(m2$coefficients, n = 1)) # felm IV fits have endog coef at the tail
  Var_2sls <- c(tail(diag(m2$robustvcv), n = 1))

  # residualised D and Z
  resD <- matrix(partialer(Y = D, X = controls, FE = FE, data = data, weights = weights), n, 1)
  resZ <- matrix(NA, n, p_iv)
  for (i in 1:p_iv) {
    resZ[, i] <- partialer(Y = Z[i], X = controls, FE = FE, data = data, weights = weights)
  }

  # main calculation
  A <- solve(t(resD) %*% resZ %*% solve(t(resZ) %*% resZ) %*% t(resZ) %*% resD) %*% (t(resD) %*% resZ)
  # point est, var
  ltz_pt_est <- beta_2sls - A %*% mu
  ltz_var_est <- Var_2sls + A %*% Sig %*% t(A)
  ltz_se <- sqrt(ltz_var_est)
  ltz_ci <- qnorm(c(0.025, 0.975), ltz_pt_est, ltz_se)
  ltz_out <- matrix(c(ltz_pt_est, ltz_se, ltz_ci), 1, 4)
  colnames(ltz_out) = c("Coef", "SE", "CI 2.5%", "CI 97.5%")
  return(ltz_out)
}
