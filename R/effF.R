#' Computes effective F statistic
#' @param ivmod fitted felm IV object with residualised Y, D, Z
#' @param noi noisy, defaults to TRUE
#' @export
eff_F_stat = function(ivmod, noi = T){
  fsmod = ivmod$stage1
  n_inst = length(fsmod$instruments)
  # overidentified
  if(n_inst > 1){
    if (noi) cat("Pflueger/Montiel-Olea Effective F statistic")
    # first stage coef and vcov
    π = fsmod$coefficients
    Σ = fsmod$robustvcv
    Z = cbind(1, fsmod$ivx)
    Q_zz = (t(Z) %*% Z)
    trace = function(x) sum(diag(x))
    Eff_F = (t(π)  %*% Q_zz  %*% π) / trace(Σ %*% Q_zz)
  }
  else{ # just identified
    if (noi) cat("Kleibergen-Paap F statistic")
    Eff_F = fsmod$rob.iv1fstat$x['F']
  }
  return(as.numeric(Eff_F))
}
