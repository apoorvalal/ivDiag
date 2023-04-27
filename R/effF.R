#' Computes effective F statistic
#' @param data dataframe
#' @param Y outcome (string)
#' @param D treatment (string)
#' @param Z instrument (string)
#' @param controls control variables (character vector)
#' @param FE fixed effects (character vector)
#' @param cl clustering column for SE (character vector)
#' @param weights weights name (string)
#' @export
eff_F = function(data, Y, D, Z, controls = NULL, FE = NULL, cl = NULL, weights = NULL # weights is a string
) {
  fmla = formula_lfe(Y = Y, W = D, Z = Z, X = controls, FE = FE, Cl = cl)
  if (is.null(weights)) {
    ivmod = lfe::felm(fmla, data = data)
  } else {
    ivmod = lfe::felm(fmla, data = data, weights = data[, weights])
  }
  # first stage model object
  fsmod = ivmod$stage1
  # instrument
  Z = fsmod$ivx
  p_iv = length(fsmod$instruments)
  if (is.null(cl) == TRUE) {
    vcv = fsmod$robustvcv
  } else {
    vcv = fsmod$clustervcv
  }
  # variance covariance matrix
  p <- nrow(vcv)
  iv_pos <- (p - p_iv + 1):p
  # first stage coef and vcov
  pi = matrix(fsmod$coefficients[iv_pos], p_iv, 1)
  Sigma = vcv[iv_pos, iv_pos, drop = FALSE]
  # instrument matrix
  Q_zz = (t(Z) %*% Z)
  eff_F = c(t(pi) %*% Q_zz %*% pi / sum(diag(Sigma %*% Q_zz)))
  return(eff_F)
}
