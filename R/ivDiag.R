#' All IV diagnostics in one
#' @param data dataframe
#' @param Y outcome (string)
#' @param D treatment (string)
#' @param Z instrument (string)
#' @param controls control variables (character vector)
#' @param FE fixed effects (character vector)
#' @param cl clustering column for SE (character vector)
#' @param weights weights name (string)
#' @param bootstrap whether to turn on bootstrap (TRUE by default)
#' @param run.AR whether to run AR test (TRUE by default)
#' @param nboots number of bootstrap reps (1000 by default)
#' @param parallel boolean for parallel bootstrap (on by default)
#' @param cores number of cores to parallelise across (defaults to all - 1)
#' @param prec precision of summary (4 by default)
#' @param seed seed
#' @return list of results
#' @importFrom lfe felm
#' @importFrom foreach foreach `%dopar%`
#' @export
ivDiag <- function(
    data, Y, D, Z, controls = NULL, FE = NULL, cl = NULL,
    weights = NULL, bootstrap = TRUE, run.AR = TRUE, 
    nboots = 1000, parallel = TRUE, cores = NULL, seed = 94305, 
    prec = 4, debug = FALSE) {
  
  if (bootstrap == FALSE & run.AR == FALSE) {
    parallel <- FALSE # no need to parallelise if not bootstrapping or running AR
  }    
  
  t0 <- Sys.time()
  set.seed(seed)
  ##############################
  # data prep
  ##############################
  for (var in unique(c(Y, D, Z, controls, FE, cl, weights))) {
    if (!var %in% names(data)) {
      stop(paste0("\"", var, "\" is not in the dataset; please double check."))
    }
  }

  # drop missingness
  data <- data[, unique(c(Y, D, Z, controls, FE, cl, weights))]
  data <- haven::zap_labels(data)
  d0 <- as.data.frame(data[complete.cases(data), ])
  n <- nrow(d0); p_iv <- length(Z)

  ### prep (using the first clustering variable in bootstrap)
  if (is.null(cl) == FALSE) { # find clusters
    d0 <- d0[order(d0[, cl[1]]), ]
    clusters <- unique(d0[, cl[1]])
    id.list <- split(1:n, d0[, cl[1]])
    ncl <- length(clusters)
  } else {
    ncl <- NULL
  }

  # IV levels, Treatment levels and Outcome levels
  nvalues <- matrix(rep(NA, 2 + p_iv), 1, 2 + p_iv)
  nvalues[, 1] <- length(unique(d0[, Y]))
  nvalues[, 2] <- length(unique(d0[, D]))
  for (i in 1:p_iv) {
    nvalues[, 2 + i] <- length(unique(d0[, Z[i]]))
  }
  colnames(nvalues) <- c(Y, D, Z)

  #####################################################
  ## point estimates - main data
  #####################################################

  # OLS
  olsfit <- OLS(data = d0, Y = Y, D = D, X = controls, FE = FE, cl = cl, weights = weights)
  OLS.Coef <- olsfit$coef
  OLS.SE <- olsfit$se
  OLS.t <- OLS.Coef / OLS.SE
  OLS.p <- (1 - pnorm(abs(OLS.t))) * 2
  OLS.df <- olsfit$df
  # 2SLS
  ivfit <- IV(data = d0, Y = Y, D = D, Z = Z, X = controls, FE = FE, cl = cl, weights = weights)
  IV.Coef <- ivfit$coef
  IV.SE <- ivfit$se
  IV.t <- IV.Coef / IV.SE
  IV.p <- (1 - pnorm(abs(IV.t))) * 2
  IV.df <- ivfit$df
  # reduced form
  rffit <- OLS(data = d0, Y = Y, D = Z, X = controls, FE = FE, cl = cl, weights = weights)
  RF.Coef <- matrix(rffit$coef, p_iv, 1)
  RF.SE <- matrix(rffit$se, p_iv, 1)
  RF.t <- RF.Coef / RF.SE
  RF.p <- (1 - pnorm(abs(RF.t))) * 2
  # first stage
  fsfit <- OLS(data = d0, Y = D, D = Z, X = controls, FE = FE, cl = cl, weights = weights)
  FS.Coef <- matrix(fsfit$coef, p_iv, 1)
  FS.SE <- matrix(fsfit$se, p_iv, 1)
  FS.t <- FS.Coef / FS.SE
  FS.p <- (1 - pnorm(abs(FS.t))) * 2
  rho <- first_stage_rho(data = d0, D = D, Z = Z, X = controls, FE = FE, weights = weights)


  #####################################################
  # Bootstrap function
  #####################################################

  if (bootstrap == TRUE) {
    message("Bootstrapping:\n")
    boot.core <- function() {
      if (is.null(cl) == TRUE) {
        # draw bootstrap sample
        smp <- sample(1:n, n, replace = TRUE)
      } else {
        # block bootstrap
        cluster.boot <- sample(clusters, ncl, replace = TRUE)
        smp <- unlist(id.list[match(cluster.boot, clusters)]) # match to locate the position of the clusterin the list
      }
      s <- d0[smp, ]
      ## init container vector for results (4 + 2 * p_iv)
      ncoef <- (2 + 2 * p_iv); res <- rep(NA, ncoef + 2)
      # first 2 cols are OLS and IV coefs -- 2
      # then Reduced Form & First Stage coefs -- p_iv + p_iv
      # then two t stats (OLS & IV) -- 2
      reg.ols <- OLS(data = s, Y = Y, D = D, X = controls, FE = FE, cl = cl, weights = weights)
      reg.iv <- IV(data = s, Y = Y, D = D, Z = Z, X = controls, FE = FE, cl = cl, weights = weights)
      res[1] <- reg.ols$coef
      res[2] <- reg.iv$coef
      res[1 + ncoef] <- (reg.ols$coef - OLS.Coef) / reg.ols$se # t stat for ols
      res[2 + ncoef] <- (reg.iv$coef - IV.Coef) / reg.iv$se # t stat for iv
      # reduced form
      reg.rf <- OLS(data = s, Y = Y, D = Z, X = controls, FE = FE, cl = cl, weights = weights)
      res[3:(2 + p_iv)] <- reg.rf$coef
      # first stage coefficients
      reg.fs <- OLS(data = s, Y = D, D = Z, X = controls, FE = FE, cl = cl, weights = weights)
      res[(3 + p_iv):(2 + p_iv * 2)] <- reg.fs$coef
      return(res)
    }
    # in case there's an error (need to comment out in debugging)
    one.boot <- function() {
      est <- try(boot.core(), silent = TRUE)
      if ('try-error' %in% class(est)) {
        est0 <- rep(NA, (4 + 2 * p_iv))
        return(est0)
      } else {
        return(est)
      }
    }
    if (debug) print(c("one boot", one.boot()))
    #####################################################
    # boot looper - series or parallel
    #####################################################
    if (parallel == FALSE) {
      boot.out = matrix(NA, nboots, 4 + 2 * p_iv)
      for (i in 1:nboots) {
        boot.out[i, ] <- one.boot()
      }
    } else { # parallel computing
      # parallelisation setup
      # let user specify ncores ; default to all - 1
      if (is.null(cores)) {
        cores <- parallel::detectCores() - 1
      }
      message("Parallelising ", nboots, " reps on ", cores, " cores \n", sep = "")
      # register
      cl.parallel <- future::makeClusterPSOCK(cores, verbose = FALSE)
      doParallel::registerDoParallel(cl.parallel)
      expfun <- c("OLS", "IV", "formula_lfe", "robustify")
      boot.out <- foreach::foreach(
        i = 1:nboots, .combine = rbind, .inorder = FALSE,
        .export = expfun,
        .packages = c("lfe")
      ) %dopar% {
        return(one.boot())
      }
      doParallel::stopImplicitCluster()
    }
    ##############################
    # post-bootstrap processing
    ##############################
    if (debug == TRUE) {
      print(boot.out)
    }
    # drop NAs
    boot.out <- boot.out[complete.cases(boot.out), ]
    # OLS and IV mean and standard errors
    boot.coefs <- boot.out[, 1:(2 + 2 * p_iv)]
    boot.tstat <- boot.out[, (3 + 2 * p_iv):(4 + 2 * p_iv)]
    bootSE <- apply(boot.coefs, 2, sd)
    OLS.boot.SE <- bootSE[1]
    IV.boot.SE <- bootSE[2]
    RF.boot.SE <- bootSE[3:(p_iv + 2)]
    FS.boot.SE <- bootSE[(3 + p_iv):(p_iv * 2 + 2)]
    # OLS and IV CIs
    CI.lvl <- c((1 - 0.95) / 2, (1 - (1 - 0.95) / 2))
    OLS.boot.ci <- quantile(boot.coefs[, 1], CI.lvl)
    IV.boot.ci <- quantile(boot.coefs[, 2], CI.lvl)
    # reduced form and 1st stage
    RF.boot.ci <- t(apply(boot.coefs[, 3:(p_iv + 2), drop = FALSE], 2, quantile, CI.lvl))
    FS.boot.ci <- t(apply(boot.coefs[, (3 + p_iv):(p_iv * 2 + 2), drop = FALSE], 2, quantile, CI.lvl))

    # Calculate first stage F
    FStat_bootout <- boot.out[, (3 + p_iv):(p_iv * 2 + 2)] # first stage
    F.boot = c((t(FS.Coef) %*% solve(var(FStat_bootout)) %*% FS.Coef) / p_iv)
    names(F.boot) <- "F.boot"

    # timing
    t1 <- Sys.time() - t0
    message("Bootstrap took", sprintf("%.3f", t1), "sec.\n")

    ## Bootstrap outputs
    # OLS
    OLS.boot.t <- OLS.Coef / OLS.boot.SE
    tmp.p <- sum(boot.coefs[, 1] > 0) / nrow(boot.coefs)
    OLS.boot.p <- ifelse(tmp.p > 0.5, 2 * (1 - tmp.p), 2 * tmp.p)
    # IV
    IV.boot.t <- IV.Coef / IV.boot.SE
    tmp.p <- sum(boot.coefs[, 2] > 0) / nrow(boot.coefs)
    IV.boot.p <- ifelse(tmp.p > 0.5, 2 * (1 - tmp.p), 2 * tmp.p)
    # Reduced form
    tmp.p <- apply(boot.coefs[, 3:(p_iv + 2), drop = FALSE] > 0, 2, sum) / nrow(boot.coefs)
    RF.boot.p <- ifelse(tmp.p > 0.5, 2 * (1 - tmp.p), 2 * tmp.p)
    # First stage
    tmp.p <- apply(boot.coefs[, (3 + p_iv):(p_iv * 2 + 2), drop = FALSE] > 0, 2, sum) / nrow(boot.coefs)
    FS.boot.p <- ifelse(tmp.p > 0.5, 2 * (1 - tmp.p), 2 * tmp.p)
    ## Bootstrap refinement
    # OLS
    ct.ols <- quantile(abs(boot.tstat[, 1]), 0.95) # critical value
    OLS.rf.ci <- c(OLS.Coef - ct.ols * OLS.SE, OLS.Coef + ct.ols * OLS.SE)
    OLS.rf.p <- sum(abs(OLS.t) < abs(boot.tstat[, 1])) / nrow(boot.tstat) # bootstrap refined p-value
    # IV
    ct.2sls <- quantile(abs(boot.tstat[, 2]), 0.95) # critical value
    IV.rf.ci <- c(IV.Coef - ct.2sls * IV.SE, IV.Coef + ct.2sls * IV.SE)
    IV.rf.p <- sum(abs(IV.t) < abs(boot.tstat[, 2])) / nrow(boot.tstat) # bootstrap refined p-value  

  } # end of bootstrap
  
  ##############################
  # prep output
  ##############################

  # F.stats
  vcov <- get.vcov(data, D, Y, Z, X = controls, FE = FE, cl = cl, weights = weights)
  F.standard <- c((t(FS.Coef) %*% solve(vcov[[1]]) %*% FS.Coef) / p_iv)
  F.robust <- c((t(FS.Coef) %*% solve(vcov[[2]]) %*% FS.Coef) / p_iv)
  if (is.null(cl) == FALSE) {
    F.cluster <- c((t(FS.Coef) %*% solve(vcov[[3]]) %*% FS.Coef) / p_iv)
  } else {
    F.cluster <- NA
  }
  F.effective <- eff_F(data, Y = Y, D = D, Z = Z, controls = controls, FE = FE, cl = cl, weights = weights)
  if (bootstrap == TRUE) {
    F_stat <- c(F.standard, F.robust, F.cluster, F.boot, F.effective)
    names(F_stat) <- c("F.standard", "F.robust", "F.cluster", "F.bootstrap", "F.effective")
  } else {
    F_stat <- c(F.standard, F.robust, F.cluster, F.effective)
    names(F_stat) <- c("F.standard", "F.robust", "F.cluster", "F.effective")  
  }
  
  # AR test
  if (run.AR == TRUE) {
    AR <- AR_test(data = data, Y = Y, D = D, Z = Z, controls = controls, FE = FE,
      cl = cl, weights = weights, prec = prec, CI = TRUE, alpha = 0.05, 
      parallel = parallel, cores = cores)
  } else {
    AR <- NULL
  }

  # calculate ratio
  if (p_iv == 1) {
    ratio <- rep(NA, 1)
    ratio[1] <- IV.Coef / OLS.Coef
  } else {
    ratio <- rep(NA, p_iv + 1)
    names(ratio) <- c("together", Z)
    ratio[1] <- IV.Coef / OLS.Coef
    # run individual IV regressions
    for (i in 1:p_iv) {
      iv.coef.tmp <- IV(data = d0, Y = Y, D = D, Z = Z[i], X = controls, FE = FE, cl = cl, weights = weights)$coef
      ratio[(i + 1)] <- iv.coef.tmp / OLS.Coef
    }
  }

  # put together (OLS and IV estimates)
  if (bootstrap == TRUE) {
    est_ols <- matrix(NA, 3, 6)
    est_ols[1, ] <- c(OLS.Coef, OLS.SE, OLS.t, OLS.Coef - 1.96 * OLS.SE, OLS.Coef + 1.96 * OLS.SE, OLS.p)
    est_ols[2, ] <- c(OLS.Coef, OLS.boot.SE, OLS.boot.t, OLS.boot.ci, OLS.boot.p)
    est_ols[3, ] <- c(OLS.Coef, OLS.SE, OLS.t, OLS.rf.ci, OLS.rf.p)
    est_2sls <- matrix(NA, 3, 6)
    est_2sls[1, ] <- c(IV.Coef, IV.SE, IV.t, IV.Coef - 1.96 * IV.SE, IV.Coef + 1.96 * IV.SE, IV.p)
    est_2sls[2, ] <- c(IV.Coef, IV.boot.SE, IV.boot.t, IV.boot.ci, IV.boot.p)
    est_2sls[3, ] <- c(IV.Coef, IV.SE, IV.t, IV.rf.ci, IV.rf.p)
    colnames(est_ols) <- colnames(est_2sls) <- c("Coef", "SE", "t", "CI 2.5%", "CI 97.5%", "p.value")
    rownames(est_ols) <- c("Analytic", "Boot.c", "Boot.t")
    rownames(est_2sls) <- c("Analytic", "Boot.c", "Boot.t")
  } else {
    est_ols <- matrix(NA, 1, 6)
    est_ols[1, ] <- c(OLS.Coef, OLS.SE, OLS.t, OLS.Coef - 1.96 * OLS.SE, OLS.Coef + 1.96 * OLS.SE, OLS.p)
    est_2sls <- matrix(NA, 1, 6)
    est_2sls[1, ] <- c(IV.Coef, IV.SE, IV.t, IV.Coef - 1.96 * IV.SE, IV.Coef + 1.96 * IV.SE, IV.p)
    colnames(est_ols) <- colnames(est_2sls) <- c("Coef", "SE", "t", "CI 2.5%", "CI 97.5%", "p.value")
    rownames(est_ols) <- c("Analytic")
    rownames(est_2sls) <- c("Analytic")
  }
  

  # tF procedure
  if (p_iv == 1) {
    # using analytic t and robust/cluster (effective) F
    tF.out <- tF(coef = IV.Coef, se = IV.SE, Fstat = F.effective)
  }

  # put together (reduced form and first stage)
  if (bootstrap == TRUE) {
    est_rf <- cbind(RF.Coef, RF.SE, RF.p, RF.boot.SE, RF.boot.ci, RF.boot.p)
    est_fs <- cbind(FS.Coef, FS.SE, FS.p, FS.boot.SE, FS.boot.ci, FS.boot.p)
    colnames(est_rf) <- colnames(est_fs) <- c("Coef", "SE", "p.value", "SE.b", "CI.b2.5%", "CI.b97.5%", "p.value.b")
    rownames(est_rf) <- rownames(est_fs) <- Z
  } else {
    est_rf <- cbind(RF.Coef, RF.SE, RF.p, RF.Coef - 1.96 * RF.SE, RF.Coef + 1.96 * RF.SE)
    est_fs <- cbind(FS.Coef, FS.SE, FS.p, FS.Coef - 1.96 * FS.SE, FS.Coef + 1.96 * FS.SE)
    colnames(est_rf) <- colnames(est_fs) <- c("Coef", "SE", "p.value", "CI.2.5%", "CI.97.5%")
    rownames(est_rf) <- rownames(est_fs) <- Z
  }  

  # save results
  if (p_iv == 1) {
    output <- list(
      # OLS and IV results
      est_ols = round(est_ols, prec),
      est_2sls = round(est_2sls, prec),
      # AR test
      AR = AR,
      # F stats
      F_stat = round(F_stat, prec),
      # first strage rho
      rho = round(rho, prec),
      # tF procedure
      tF = round(tF.out, prec),
      # reduced form and first stage
      est_rf = round(est_rf, prec),
      est_fs = round(est_fs, prec),
      # number of instruments
      p_iv = p_iv,
      # number of observations
      N = n,
      # number of clusters
      N_cl = ncl,
      # degrees of freedom
      df = IV.df,
      # number of unique values
      nvalues = nvalues
    )
  } else { # p_iv >1
    output <- list(
      # OLS and IV results
      est_ols = round(est_ols, prec),
      est_2sls = round(est_2sls, prec),
      # AR test
      AR = AR,
      # F stats
      F_stat = round(F_stat, prec),
      # first strage rho
      rho = round(rho, prec),
      # reduced form and first stage
      est_rf = round(est_rf, prec),
      est_fs = round(est_fs, prec),
      # number of instruments
      p_iv = p_iv,
      # number of observations
      N = n,
      # number of clusters
      N_cl = ncl,
      # degrees of freedom
      df = IV.df,
      # number of unique values
      nvalues = nvalues
    )
  }
  class(output) <- "ivDiag"
  return(output)
}


# %%
