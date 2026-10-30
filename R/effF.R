# Relative tolerance for "collinear": a column is treated as linearly dependent
# on the columns before it when what is left of it, after projecting them out,
# has a norm below this fraction of its own norm (the default of lm()/qr(),
# about 1 - R^2 < 1e-14). Weak but non-degenerate instruments are far above it.
EFF_F_COLLINEAR_TOL = 1e-7
# An instrument whose sd is below this fraction of its largest absolute value
# is constant (up to rounding) and cannot be standardised.
EFF_F_CONSTANT_TOL = 1e-12
EFF_F_COLLINEAR_MSG = paste0(
  "eff_F(): the instruments are collinear with each other or with the ",
  "controls/fixed effects after partialling them out; drop the redundant ",
  "instrument(s).")

# TRUE if the last p_iv columns of felm's first stage (the instruments) are
# linearly independent of each other and of the other columns (intercept,
# controls, fixed effects). Aliased controls are ignored: they do not matter.
# lfe's own handling of such columns is not reliable: it sets the coefficient
# of one column of a dependent set to NaN, which can be a control rather than
# the instrument, so the NaN pattern cannot be used. Instead the instruments
# are residualised on the controls with a rank-revealing QR (weighted,
# fixed effects already demeaned by lfe) and compared with their raw size.
instruments_identified = function(fsmod, p_iv, tol = EFF_F_COLLINEAR_TOL) {
  cX = fsmod$cX  # demeaned design matrix: intercept, controls, instruments
  X = fsmod$X    # the same, not demeaned
  if (is.null(cX) || is.null(X)) return(TRUE)  # cannot check; later checks apply
  sw = if (is.null(fsmod$weights)) 1 else sqrt(fsmod$weights)
  iv = (ncol(cX) - p_iv + 1):ncol(cX)
  rel_size = function(a, b) sqrt(colSums((a * sw)^2)) / sqrt(colSums((b * sw)^2))
  # controls that are (numerically) all zero after demeaning, e.g. constants
  # with fixed effects, are dropped; they carry no information
  ctrl = cX[, -iv, drop = FALSE]
  keep = rel_size(ctrl, X[, -iv, drop = FALSE]) >= tol
  ctrl = ctrl[, which(keep), drop = FALSE]
  Zc = cX[, iv, drop = FALSE] * sw
  if (ncol(ctrl) > 0) {
    # pivoting QR sets aliased controls aside
    Zt = qr.resid(qr(ctrl * sw, tol = tol), Zc)
  } else {
    Zt = Zc
  }
  # size of each partialled-out instrument relative to the raw instrument
  ratio = sqrt(colSums(Zt^2)) / sqrt(colSums((X[, iv, drop = FALSE] * sw)^2))
  if (anyNA(ratio) || any(ratio < tol)) return(FALSE)
  # instruments independent of each other (scale-free, relative to own size)
  qr(Zt, tol = tol)$rank == p_iv
}

#' Computes effective F statistic
#' @param data dataframe
#' @param Y outcome (string)
#' @param D treatment (string)
#' @param Z instrument (string)
#' @param controls control variables (character vector)
#' @param FE fixed effects (character vector)
#' @param cl clustering column for SE (character vector)
#' @param weights weights name (string)
#' @param prec precision for output
#' @details The instruments must be identified: after the intercept, controls
#'   and fixed effects are partialled out, they must be linearly independent and
#'   not (numerically) zero. Otherwise \code{eff_F()} stops with an error.
#'   Collinear controls are harmless and are ignored. The instruments are centred
#'   and scaled internally; the effective F does not depend on their scale or mean.
#' @export
eff_F = function(data, Y, D, Z, controls = NULL, FE = NULL, 
  cl = NULL, weights = NULL, prec = 4
  ) {
  fmla = formula_lfe(Y = Y, D = D, Z = Z, X = controls, FE = FE, cl = cl)
  # The effective F is invariant to any nonsingular linear transformation of
  # the instruments (pi and Q change together), and a shift is absorbed by the
  # intercept or the fixed effects. So centre and scale each instrument (on its
  # non-missing rows) in a local copy of the data: instruments on very
  # different scales or with a huge mean otherwise make lfe's Cholesky step
  # alias them or lose precision. A constant instrument cannot be scaled.
  # Only numeric and logical instruments: a factor or character instrument is
  # left as it is, so that felm expands it into dummies.
  for (z in Z) {
    if (!(is.numeric(data[[z]]) || is.logical(data[[z]]))) next
    v = as.numeric(data[[z]])
    s = stats::sd(v, na.rm = TRUE)
    if (!is.finite(s) || s <= EFF_F_CONSTANT_TOL * max(abs(v), na.rm = TRUE)) {
      stop(EFF_F_COLLINEAR_MSG, call. = FALSE)
    }
    data[[z]] = (v - mean(v, na.rm = TRUE)) / s
  }
  # with a single instrument that is a linear combination of the controls,
  # felm itself fails inside chol2inv() ("... the inverse cannot be computed")
  ivmod = tryCatch({
    if (is.null(weights)) {
      lfe::felm(fmla, data = data, keepX = TRUE, keepCX = TRUE)
    } else {
      lfe::felm(fmla, data = data, weights = data[, weights],
                keepX = TRUE, keepCX = TRUE)
    }
  }, error = function(e) {
    if (grepl("inverse cannot be computed", conditionMessage(e), fixed = TRUE)) {
      stop(EFF_F_COLLINEAR_MSG, call. = FALSE)
    }
    stop(e)
  })
  # first stage model object
  fsmod = ivmod$stage1
  p_iv = length(fsmod$instruments)
  # stop early if the instruments are not identified (lfe only warns, and may
  # alias a control instead of the instrument, leaving a number that is wrong)
  if (!instruments_identified(fsmod, p_iv)) stop(EFF_F_COLLINEAR_MSG, call. = FALSE)
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
  # Q = Z~' W Z~, Z~ = instruments with intercept, controls and FE partialled
  # out (weighted), on felm's sample: inverse of the instrument block of
  # (X'WX)^{-1} (Frisch-Waugh-Lovell); instruments are the last columns
  Q_zz = tryCatch(solve(fsmod$inv[iv_pos, iv_pos, drop = FALSE]),
                  error = function(e) NULL)
  if (is.null(Q_zz) || anyNA(pi) || anyNA(Sigma) || anyNA(Q_zz)) {
    stop(EFF_F_COLLINEAR_MSG, call. = FALSE)
  }
  eff_F = c(t(pi) %*% Q_zz %*% pi / sum(diag(Sigma %*% Q_zz)))
  return(round(eff_F, prec))
}
