#' Computes effective F statistic
#' @param ivmod fitted felm IV object with residualised Y, D, Z
#' @param noi noisy, defaults to TRUE
#' @export
eff_F_stat = function(ivmod) {
  # trace function
  trace = \(x) sum(diag(x))
  # first stage model object
  fsmod = ivmod$stage1
  scalingFactor = (fsmod$N - fsmod$p) / fsmod$N
  # instrument
  Z = fsmod$ivx
  nZ = length(fsmod$instruments)
  if (nZ == 1) { # just identified
    Eff_F = ivmod$stage1$rob.iv1fstat[[1]]['F'] # already computed in felm
  } else {
    # variance covariance matrix
    vcv = fsmod$robustvcv
    nrV = nrow(vcv)
    startid = (nrV - nZ + 1)
    Σ = vcv[startid:nrV, startid:nrV, drop = F]
    # first stage coef and vcov
    π = fsmod$coefficients %>% .[startid:nrV]
    # instrument matrix
    Q_zz = (t(Z) %*% Z)
    Eff_F = trace(Σ %*% Q_zz) * scalingFactor
  }
  return(Eff_F)
}
