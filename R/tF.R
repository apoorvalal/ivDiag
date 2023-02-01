# LMMP's tf procedure
#' @param coef 2SLS coefficient
#' @param se SE for 2SLS estimate
#' @param Fstat first stage F
tF <- function(coef, se, Fstat) {
  tstat <- coef / se
  F.sqrt <- sqrt(Fstat)
  F0.sqrt <- seq(2, 10.3, 0.1)
  cF0 <- c(
    18.66, 9.74, 7.37, 6.18, 5.43, 4.92, 4.54, 4.25, 4.01, 3.82, 3.65, 3.51, 3.39, 3.29, 3.19, 3.11, 3.03,
    2.97, 2.91, 2.85, 2.80, 2.75, 2.71, 2.67, 2.63, 2.60, 2.57, 2.54, 2.51, 2.48, 2.46, 2.43, 2.41, 2.39, 2.37,
    2.35, 2.33, 2.32, 2.30, 2.29, 2.27, 2.26, 2.24, 2.23, 2.22, 2.21, 2.20, 2.19, 2.17, 2.16, 2.16, 2.15, 2.14,
    2.13, 2.12, 2.11, 2.10, 2.10, 2.09, 2.08, 2.08, 2.07, 2.06, 2.06, 2.05, 2.04, 2.04, 2.03, 2.03, 2.02, 2.02,
    2.01, 2.01, 2.00, 2.00, 1.99, 1.99, 1.99, 1.98, 1.98, 1.97, 1.97, 1.97, 1.96
  )
  if (F.sqrt <= 2) {
    cF <- cF0[1]
  } else if (F.sqrt >= 10.3) {
    cF <- 1.96
  } else {
    pos.lower <- max(which(F0.sqrt < F.sqrt))
    pos.upper <- pos.lower + 1
    h1 <- abs(F.sqrt - F0.sqrt[pos.lower])
    h2 <- abs(F.sqrt - F0.sqrt[pos.upper])
    # critical value
    cF <- (cF0[pos.upper] * h1 + cF0[pos.lower] * h2) / (h1 + h2) # weighted average
  }
  ci <- c(coef - cF * se, coef + cF * se)
  p <- (1 - pnorm(abs(tstat) / (cF / 1.96))) * 2 # adjusted p value
  out <- c(Fstat, cF, coef, se, tstat, ci, p)
  names(out) <- c("F", "cF", "Coef", "SE", "t", "CI2.5%", "CI97.5%", "p-value")
  return(out)
}
