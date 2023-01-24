#' Weak-IV Robust Anderson Rubin Test
# fork of `ivmodel::AR.test` that accepts FELM objects and fits models with FEs
#' @param data dataframe
#' @param Y outcome (string)
#' @param D treatment (string)
#' @param Z instrument (string)
#' @param X control variables (character vector)
#' @param FE fixed effects (character vector)
#' @param Cl Cluster name
#' @param weights weighting vector
#' @param prec precision of CI in string
#' @param hetF Use Wald instead of F statistic to account for heteroskedasticity?
#' @importFrom lfe felm
#' @export

AR_test = function(
    data, Y, D, Z, X, FE,
    weights = NULL, Cl = NULL, prec = 3, hetF = FALSE) {
  # keep rows with complete data
  data <- data[, c(Y, D, Z, X, FE, weights)]
  data <- data[complete.cases(data), ]

  # fit IV model in FELM to get correct DoFs
  mod = suppressWarnings(
    felm(formula_lfe(Y = Y, W = D, X = X, Z = Z, D = FE, C = Cl), data = data) %>% robustify()
  )
  # residualise
  Ytil = suppressWarnings(partialer(Y, X = X, FE = FE, data = data, weights = weights))
  Dtil = suppressWarnings(partialer(D, X = X, FE = FE, data = data, weights = weights))

  if (length(Z) == 1) { # scalar partialling out
    Ztil = suppressWarnings(partialer(Z, X = X, FE = FE, data = data, weights = weights))
  } else { # partialling out sequentially
    Ztil = matrix(NA, nrow = nrow(data), ncol = length(Z))
    for (i in 1:length(Z)) {
      Ztil[, i] = suppressWarnings(partialer(Z[i],
        X = X, FE = FE,
        data = data, weights = weights
      ))
    }
  }

  alpha = 0.05; n = length(Ztil);  k = mod$p;  l = length(Z)
  ZtilQR = qr(Ztil)

  # compute F
  tmp = Ytil - 0 * Dtil # test null of \beta = 0
  # this is homoskedastic
  Fstat = c(sum(qr.fitted(ZtilQR, tmp)^2)) / c(sum(tmp^2) - sum(qr.fitted(ZtilQR, tmp)^2)) * (n - k - l) / l

  # Replace F with Wald test because of heteroskedasticity
  if (hetF) {
    wtype <- if (!(is.null(Cl))) "cluster" else "robust"
    Fstat = waldtest(mod, "endovars", type = wtype)
  }

  p.value = 1 - pf(Fstat, df1 = l, df2 = n - k - l)
  # CI construction ingredients
  cval = qf(1 - alpha, df1 = l, df2 = n - k - l) * l / (n - k - l)

  coef.beta0sq = cval * sum(Dtil^2) - (cval + 1) * sum(qr.fitted(ZtilQR, Dtil)^2)
  coef.beta0 = -2 * cval * sum(Dtil * Ytil) + 2 * (cval + 1) * sum(Dtil * qr.fitted(ZtilQR, Ytil))
  coef.constant = cval * sum(Ytil^2) - (cval + 1) * sum(qr.fitted(ZtilQR, Ytil)^2)

  Delta = coef.beta0^2 - 4 * coef.constant * coef.beta0sq
  ci = matrix(NA, ncol = 2); colnames(ci) <- c("lower", "upper")

  # CI construction gymnastics
  if (coef.beta0sq == 0) {
    if (coef.beta0 > 0) {
      info = c("[", round(-coef.constant / coef.beta0, prec), ",Infinity)")
      ci[1, ] = c(-coef.constant / coef.beta0, Inf)
    }
    if (coef.beta0 < 0) {
      info = c("(-Infinity,", round(-coef.constant / coef.beta0, prec), "]");
      ci[1, ] = c(-Inf, -coef.constant / coef.beta0)
    }
    if (coef.beta0 == 0) {
      if (coef.constant >= 0) {
        info = "Whole Real Line"
        ci[1, ] = c(-Inf, Inf)
      }
      if (coef.constant < 0) {
        info = "Empty Set"
      }
    }
  } else if (coef.beta0sq != 0) {
    if (Delta <= 0) {
      if (coef.beta0sq > 0) {
        info = "Whole Real Line"
        ci[1, ] = c(-Inf, Inf)
      }
      if (coef.beta0sq < 0) {
        info = "Empty Set"
      }
    }
    if (Delta > 0) {
      # Roots of quadratic equation
      root1 = (-coef.beta0 + sqrt(Delta)) / (2 * coef.beta0sq)
      root2 = (-coef.beta0 - sqrt(Delta)) / (2 * coef.beta0sq)
      upper.root = max(root1, root2)
      lower.root = min(root1, root2)
      if (coef.beta0sq < 0) {
        info = paste("[", round(lower.root, prec), ", ", round(upper.root, prec), "]", sep = "")
        ci[1, ] = c(lower.root, upper.root)
      }
      if (coef.beta0sq > 0) {
        info = paste("(-Infinity,", round(lower.root, prec), "] union [", round(upper.root, prec), ",Infinity)")
        ci[1, ] = c(-Inf, lower.root)
        ci <- rbind(ci, c(upper.root, Inf))
      }
    }
  }
  return(list(
    Fstat = Fstat,
    df = c(l, n - k - l),
    p.value = p.value,
    ci = ci,
    ci.info = info
  ))
}
