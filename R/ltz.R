# %% Local to zero IV estimate with direct effect of instrument set to μ
#' Estimates Local-to-Zero IV coefficients and SEs for a single instrument
#' @param data dataframe
#' @param Y outcome (string)
#' @param D treatment (string)
#' @param Z instrument (string)
#' @param controls control variables (character vector)
#' @param FE fixed effects (character vector)
#' @param cl clustering column for SE (character vector)
#' @param prior prior point and SD estimate of direct effect of instrument on outcome
#' @param prec precision of summary (4 by default)
#' @export
ltz = function(data, Y, D, Z, controls, FE = NULL, cl = NULL, weights = NULL, prior, prec = 4) {
  # requirements: iv coef and se, residualized Z and D, and prior on μ, Sigma
  
  # drop missingness
  data <- data[, unique(c(Y, D, Z, controls, FE, cl, weights))]
  data <- haven::zap_labels(data)
  data <- as.data.frame(data[complete.cases(data), ])
  n <- nrow(data); p_iv <- length(Z)
  names(prior) <- c("Mean", "SD")
  mu <- prior[1]
  Sig <- prior[2]^2

  # IV fit
  fmla = formula_lfe(Y = Y, D = D, Z = Z, X = controls, FE = FE, cl = cl)
  if (is.null(weights)) {
    m2 = lfe::felm(fmla, data = data)
  } else {
    m2 = lfe::felm(fmla, data = data, weights = data[, weights])
  }
  iv_beta <- c(tail(m2$coefficients, n = 1)) # felm IV fits have endog coef at the tail
  if (is.null(cl) == TRUE) {
    iv_se <- tail(m2$rse, n = 1)
  } else {
    iv_se <- tail(m2$cse, n = 1)
  }
  iv_Var <- iv_se^2
  iv_t <- iv_beta / iv_se
  iv_ci <- qnorm(c(0.025, 0.975), iv_beta, iv_se)
  iv_p <- (1 - pnorm(abs(iv_t))) * 2
  iv_out <- c(iv_beta, iv_se, iv_t, iv_ci, iv_p)

  # residualised D and Z
  resD <- matrix(partialer(Y = D, X = controls, FE = FE, data = data, weights = weights), n, 1)
  resZ <- matrix(NA, n, p_iv)
  for (i in 1:p_iv) {
    resZ[, i] <- partialer(Y = Z[i], X = controls, FE = FE, data = data, weights = weights)
  }

  # main calculation
  A <- solve(t(resD) %*% resZ %*% solve(t(resZ) %*% resZ) %*% t(resZ) %*% resD) %*% (t(resD) %*% resZ)
  # point est, var
  ltz_pt_est <- iv_beta - A %*% mu
  ltz_var_est <- iv_Var + A %*% Sig %*% t(A)
  ltz_se <- sqrt(ltz_var_est)
  ltz_t <- ltz_pt_est / ltz_se
  ltz_ci <- qnorm(c(0.025, 0.975), ltz_pt_est, ltz_se)
  ltz_p <- (1 - pnorm(abs(ltz_t))) * 2
  ltz_out <- c(ltz_pt_est, ltz_se, ltz_t, ltz_ci, ltz_p)
  names(iv_out) <- names(ltz_out) <- c("Coef", "SE", "t", "CI 2.5%", "CI 97.5%", "p-value")
  
  out <- list(iv = round(iv_out, prec),
    ltz = round(ltz_out, prec), 
    prior = round(prior, prec))
  class(out) <- "ltz"
  return(out)
}
