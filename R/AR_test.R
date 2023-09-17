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
#' @param CI whether to conduct inversion to get CI
#' @param alpha level of statitical significance
#' @param parallel parallel computing
#' @param cores number of cores
#' @importFrom lfe felm
#' @export

AR_test = function(
    data, Y, D, Z, controls, FE = NULL, cl = NULL,
    weights = NULL, prec = 4, CI = TRUE, alpha = 0.05,
    parallel = NULL, cores = NULL) {
  # keep rows with complete data
  X <- controls
  data <- data[, c(Y, D, Z, X, FE, cl, weights)]
  data <- data[complete.cases(data), ]
  p_iv <- length(Z)

  # parallelising
  if (is.null(parallel) == TRUE) {
    if (nrow(data) > 5000 & CI == TRUE) {
      parallel <- TRUE
    } else {
      parallel <- FALSE
    }
  }

  # Run IV regression once
  ivfit <- IV(data = data, Y = Y, D = D, Z = Z, X = X, FE = FE, cl = cl, weights = weights)
  beta_iv <- ivfit$coef
  se_iv <- ivfit$se
  beta_seq <- c(
    seq(beta_iv - 10 * se_iv, beta_iv - 3.1 * se_iv, length.out = 100),
    seq(beta_iv - 3 * se_iv, beta_iv + 3 * se_iv, by = 0.02 * se_iv),
    seq(beta_iv + 3.1 * se_iv, beta_iv + 10 * se_iv, length.out = 100)
  )
  ngrid <- length(beta_seq)

  # residualise
  Ytil = suppressWarnings(partialer(Y, X = X, FE = FE, data = data, weights = weights))
  Dtil = suppressWarnings(partialer(D, X = X, FE = FE, data = data, weights = weights))
  Ztil = matrix(NA, nrow = nrow(data), ncol = p_iv)
  for (i in 1:length(Z)) {
    Ztil[, i] = suppressWarnings(partialer(Z[i], X = X, FE = FE, data = data, weights = weights))
  }
  d <- cbind.data.frame(Ytil, Dtil, Ztil, data[, c(cl, weights)])
  colnames(d) <- c(Y, D, Z, cl, weights)

  # reduced form
  fmla = formula_lfe(Y = Y, D = Z, X = NULL, FE = NULL, cl = cl)
  
  
  # function to test each value on the real line
  one.AR <- function(beta, d, Y, D, Z, cl, weights) {
    d[, Y] <- d[, Y] - beta * d[, D]
    p_iv <- length(Z)
    reg <- OLS(data = d, Y = Y, D = Z, X = NULL, FE = NULL, cl = cl, weights = weights)
    coef <- reg$coef
    vcov <- reg$vcv
    df2 <- reg$df
    df1 <- p_iv
    Fstat <- c((t(coef) %*% solve(vcov) %*% coef) / p_iv)
    pval <- pf(Fstat, df1, df2, lower.tail = FALSE)
    output <- c(Fstat, df1, df2, pval)
    names(output) <- c("F", "df1", "df2", "p")
    return(output)
  } 

  # AR test
  Fstat <- one.AR(0, d = d, Y = Y, D = D, Z = Z, cl = cl, weights = weights)

  
  # Confidence intervals
  if (CI == TRUE) {
    message("AR Test Inversion...\n")
    if (parallel == FALSE) {
      accept <- rep(NA, ngrid)
      for (i in 1:ngrid) {    
        test.out <-  one.AR(beta_seq[i], d = d, Y = Y, D = D, Z = Z, cl = cl, weights = weights) 
        accept[i] <- ifelse(test.out[4] >= alpha, 1, 0)
      }
    } else {
      # parallel computing
      if (is.null(cores)) {
        cores <- parallel::detectCores() - 1
      }
      message("Parallelising on ", cores, " cores \n\n", sep = "")
      # register
      cl.parallel <- future::makeClusterPSOCK(cores, verbose = FALSE)
      doParallel::registerDoParallel(cl.parallel)
      expfun <- c("OLS", "IV", "formula_lfe", "OLS")
      accept <- foreach(
        i = 1:ngrid, .combine = c, .inorder = FALSE,
        .export = expfun,
        .packages = c("lfe")
      ) %dopar% {
        test.out <-  one.AR(beta_seq[i], d = d, Y = Y, D = D, Z = Z, cl = cl, weights = weights) 
        return(ifelse(test.out[4] >= alpha, 1, 0))
      }
      doParallel::stopImplicitCluster()
    }
    # summarize ("accept" means reject the null)
    bounded <- FALSE
    if (sum(accept) == ngrid) {
      ci <- c(-Inf, Inf)
      ci.print <- "(-Inf, Inf)" # all accepted
    } else if (sum(accept) == 0) {
      ci <- NA
      ci.print <- "empty"
    } else if (accept[1] == 0 && accept[ngrid] == 0) { # e.g. 0 0 0 1 1 1 0 0 0
      betas <- range(beta_seq[accept == 1])
      ci <- round(betas, prec)
      ci.print <- paste0("[", sprintf(paste0("%.", prec, "f"), betas[1]), ", ", sprintf(paste0("%.", prec, "f"), betas[2]), "]") # bounded interval
      bounded <- TRUE
    } else if (accept[1] == 1 && accept[ngrid] == 1) {
      betas <- round(range(beta_seq[accept == 0]), prec) # e.g. 1 1 1 1 0 0 0 1 1 1
      ci <- c(-Inf, betas[1], betas[2], Inf)
      ci.print <- paste0("(-Inf, ", sprintf(paste0("%.", prec, "f"), betas[1]), "] Union [", sprintf(paste0("%.", prec, "f"), betas[2]), ", Inf)")
    } else if (accept[1] == 0 && accept[ngrid] == 1) {
      betas <- round(range(beta_seq[accept == 1]), prec) # e.g. 0 0 0 1 1 1 1 1
      ci <- c(betas[1], Inf)
      ci.print <- paste0("[", sprintf(paste0("%.", prec, "f"), betas[1]), ", Inf)")
    } else if (accept[1] == 1 && accept[ngrid] == 0) {
      betas <- round(range(beta_seq[accept == 1]), prec) # e.g. 1 1 1 1 1 0 0 0
      ci <- c(-Inf, betas[2])
      ci.print <- paste0("(-Inf, ", sprintf(paste0("%.", prec, "f"), betas[2]), "]")
    }
  }

  ## Print acceptance ##
  # accept <- cbind(beta_seq, accept)
  # rownames(accept) <- NULL
  # colnames(accept) <- c("beta", "accept")
  
  # output
  if (CI == TRUE) {
    out <- list(Fstat = round(Fstat, prec), ci.print = ci.print, ci = ci, bounded = bounded)
  } else {
    out <- list(Fstat = round(Fstat, prec))
  }
  return(out)
}
