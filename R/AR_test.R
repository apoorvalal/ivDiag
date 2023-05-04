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
  fmla = formula_lfe(Y = Y, W = Z, X = NULL, FE = NULL, Cl = cl)
  if (is.null(weights) == TRUE) {
    m1 = lfe::felm(fmla, data = d)
  } else {
    m1 = lfe::felm(fmla, data = d, weights = d[, weights])
  }
  m1 <- robustify(m1)
  s <- summary(m1)
  Fstat <- s$F.fstat

  # Confidence intervals
  if (CI == TRUE) {
    message("AR Test Inversion...\n")
    if (parallel == FALSE) {
      accept <- rep(NA, ngrid)
      for (i in 1:ngrid) {
        d[, Y] <- Ytil - beta_seq[i] * Dtil
        if (is.null(weights) == TRUE) {
          m2 = lfe::felm(fmla, data = d)
        } else {
          m2 = lfe::felm(fmla, data = d, weights = d[, weights])
        }
        m2 <- robustify(m2)
        s2 <- summary(m2)
        accept[i] <- ifelse(s2$pval >= alpha, 1, 0)
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
      expfun <- c("OLS", "IV", "formula_lfe", "robustify")
      accept <- foreach(
        i = 1:ngrid, .combine = c, .inorder = FALSE,
        .export = expfun,
        .packages = c("lfe")
      ) %dopar% {
        d[, Y] <- Ytil - beta_seq[i] * Dtil
        if (is.null(weights) == TRUE) {
          m2 = lfe::felm(fmla, data = d)
        } else {
          m2 = lfe::felm(fmla, data = d, weights = d[, weights])
        }
        m2 <- robustify(m2)
        s2 <- summary(m2)
        return(ifelse(s2$pval >= alpha, 1, 0))
      }
      doParallel::stopImplicitCluster()
    }
    # summarize
    bounded <- FALSE
    if (sum(accept) == ngrid) {
      ci <- c(-Inf, Inf)
      ci.print <- "(-Inf, Inf)" # all accepted
    } else if (sum(accept) == 0) {
      ci <- NA
      ci.print <- "empty"
    } else if (accept[1] == 0 && accept[ngrid] == 0) {
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
  
  # output
  if (CI == TRUE) {
    out <- list(Fstat = round(Fstat, prec), ci.print = ci.print, ci = ci, bounded = bounded)
  } else {
    out <- list(Fstat = round(Fstat, prec))
  }
  return(out)
}
