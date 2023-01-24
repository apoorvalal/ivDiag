# %% Local to zero IV estimate with direct effect of instrument set to μ
#' Estimates Local-to-Zero IV coefficients and SEs for a single instrument
#' @param df dataframe
#' @param mu prior point estimate of direct effect of instrument on outcome
#' @param W treatment (string)
#' @param Z instrument (string)
#' @param sig variance of above prior
#' @param ivmod felm iv model object
#' @export
ltz = function(df, W, Z, mu, sig, ivmod) {
  # requiremens: prior on μ, σ,  ivmod object, and 3 column data frame with
  # residualised treatment, Y, Z
  # felm stores IV coefficient as last value
  coefs = summary(ivmod)$coefficients; iv_ests = coefs[nrow(coefs), 1:2]
  β_tsls = iv_ests[1]; Σ_tsls = iv_ests[2]^2
  # prior
  Ω = sig
  # data
  Z = matrix(df[[Z]]); X = matrix(df[[W]])
  # main calculation
  A = solve(t(X) %*% Z %*% solve(t(Z) %*% Z) %*% t(Z) %*% X) %*% (t(X) %*% Z)
  # point est, var
  ltz_pt_est = β_tsls - A %*% mu
  ltz_var_est = Σ_tsls + A %*% Ω %*% t(A)
  ltz_se = sqrt(ltz_var_est)
  ltz_ci = qnorm(c(0.025, 0.975), ltz_pt_est, ltz_se)
  ltz_out = list(est = ltz_pt_est, se = ltz_se, ci = ltz_ci)
  return(ltz_out)
}
