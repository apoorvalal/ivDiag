#' bootstrap IV estimator
#' @param data dataframe
#' @param Y outcome (string)
#' @param D treatment (string)
#' @param Z instrument (string)
#' @param controls control variables (character vector)
#' @param FE fixed effects (character vector)
#' @param cl clustering column for SE (character vector)
#' @param weights weights name (string)
#' @param nboots number of bootstrap reps (1000 by default)
#' @param parallel boolean for parallel bootstrap (on by default)
#' @param cores number of cores to parallelise across (defaults to all - 1)
#' @param prec precision of summary (4 by default)
#' @param seed seed
#' @return list of results
#' @importFrom lfe felm
#' @importFrom foreach foreach `%dopar%`
#' @export
boot_IV <- function(data, Y, D, Z, controls=NULL, FE = NULL, cl = NULL,
    weights = NULL, nboots = 1000, parallel = TRUE, seed = 94305, cores = NULL,
    prec = 4, debug = FALSE, public = TRUE) {
    ## Bootstrap OLS and IV SE/CI + F Stat for a single-instrument, single-treatment setting
    t0 <- Sys.time()
    set.seed(seed)
    ##############################
    # data prep
    ##############################
    # drop missingness
    data <- data[, c(Y, D, Z, controls, FE, cl, weights)]
    d0 <- data[complete.cases(data), ]
    n <- nrow(d0); p_iv <- length(Z)
    if (debug) print(c("p_iv :", p_iv))
    # fit first stage and store
    out0 <- first_stage_coefs(data = d0, D = D, Z = Z,
        X = controls, FE = FE, weights = weights)
    ## Sensitivity analysis - bias threshold
    rho  <- first_stage_rho(data = d0, D = D, Z = Z, X = controls, FE = FE, weights = weights)
    x_tilde = partialer(Y = D, X = controls, FE = FE, data = d0, weights = weights)
    sig_x = sd(x_tilde)
    # error variance in naive OLS
    fmla = formula_lfe(Y = Y, W = D, X = controls, D = FE, C = cl)
    if (is.null(weights)==TRUE) {
      m1 = robustify(lfe::felm(fmla, d0))
    } else {
      m1 = robustify(lfe::felm(fmla, d0, weights = d0[,weights]))
    }
    sig_epsi = sd(m1$residuals)
    ### prep
    fs_coefs0 <- matrix(out0, p_iv, 1)
    if (debug) print(c("fs_coefs0", fs_coefs0))
    if (is.null(cl)==FALSE) { # find clusters
        d0 <- d0[order(d0[,cl]),]
        clusters <- unique(d0[,cl])
        id.list <- split(1:n,d0[,cl])
        ncl <- length(clusters)
    } else {
        ncl <- NULL
    }
    #####################################################
    # Bootstrap function
    #####################################################
    cat("Bootstrapping:\n")
    boot.core <- function(){
        if (is.null(cl)==TRUE) {
          # draw bootstrap sample
          smp <- sample(1:n, n, replace=TRUE)
        } else {
          # block bootstrap
          cluster.boot <- sample(clusters, ncl, replace=TRUE)
          smp <- unlist(id.list[match(cluster.boot, clusters)])   # match to locate the position of the clusterin the list
        }
        s <- d0[smp,]
        # init container vector for results
        res <- rep(NA, (2 + p_iv)) # first 2 cols are OLS and IV estimates, last p_iv are first stage IV coefs
        res[1] <- OLS(data=s, D=D, Y=Y, X=controls, FE=FE, cl=cl, weights=weights)$coef
        res[2] <- IV(data=s, D=D, Y=Y, Z=Z, X=controls, FE=FE, cl=cl, weights=weights)$coef
        # first stage coefficients
        res[3:(2+p_iv)] <- first_stage_coefs(data=s, D=D, Z=Z, X=controls, FE = FE, weights = weights)
        return(res)
    }
    # in case there's an error
    one.boot <- function() {
      est <- try(boot.core(), silent = TRUE)
      if ('try-error' %in% class(est)) {
        est0 <- rep(NA, (2 + p_iv))
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
        coefs.boot = matrix(NA, nboots, 2 + p_iv)
        for (i in 1:nboots) {
            coefs.boot[i,] <- one.boot()
        }
    } else { # parallel computing
        # parallelisation setup
        # let user specify ncores ; default to all - 1
        if (is.null(cores)) {cores <- parallel::detectCores() - 1}
        cat("Parallelising ", nboots, " reps on ", cores, " cores \n",sep="")
        # register
        cl.parallel <- future::makeClusterPSOCK(cores, verbose = FALSE)
        doParallel::registerDoParallel(cl.parallel)
        expfun <- c("OLS", "IV", "first_stage_coefs", "formula_lfe", "robustify")
        coefs.boot <- foreach (i=1:nboots,.combine=rbind, .inorder=FALSE,
            .export = expfun,
            .packages = c("lfe")) %dopar% {
            return(one.boot())
        }
        doParallel::stopImplicitCluster()
    }
    ##############################
    # post-bootstrap processing
    ##############################
    if (debug) {return(coefs.boot)}
    # drop NAs
    coefs.boot <- coefs.boot[complete.cases(coefs.boot),]
    # OLS and IV mean and standard errors
    OLS_IV_bootout <- coefs.boot[, 1:2]
    bootMeans <- apply(OLS_IV_bootout, 2, mean)
    bootSE <- apply(OLS_IV_bootout, 2, sd)
    # OLS and IV CIs
    CI.lvl <- c((1-0.95)/2, (1-(1-0.95)/2))
    ols_ci <- round(quantile(OLS_IV_bootout[, 1], CI.lvl), 3)
    iv_ci  <- round(quantile(OLS_IV_bootout[, 2], CI.lvl), 3)
    # Calculate F'
    FStat_bootout <- coefs.boot[, 3:ncol(coefs.boot)]
    F.boot = c((t(fs_coefs0) %*% solve(var(FStat_bootout)) %*% fs_coefs0)/p_iv)
    names(F.boot) <- "F.boot"
    # timing
    t1 <- Sys.time() - t0
    cat("Bootstrap took", sprintf("%.3f",t1), "sec.\n\n")

    ##############################
    # prep output
    ##############################
    # point estimates - main data
    olsfit <- OLS(data = d0, D=D, Y=Y, X=controls, FE=FE, cl=cl, weights=weights)
    OLS.Coef <- olsfit$coef
    OLS.SE <- olsfit$se
    ivfit <- IV(data = d0, D=D, Y=Y, Z=Z, X=controls, FE=FE, cl=cl, weights=weights)
    IV.Coef <- ivfit$coef
    IV.SE <- ivfit$se
    # put together
    est_ols <- c(OLS.Coef, OLS.SE, bootSE[1], ols_ci)
    est_2sls <- c(IV.Coef, IV.SE, bootSE[2], iv_ci)
    names(est_ols) <-  c("Coef", "SE.t", "SE.b", "CI.b 2.5%", "CI.b 97.5%")
    names(est_2sls) <- c("Coef", "SE.t", "SE.b", "CI.b 2.5%", "CI.b 97.5%")
    # F.stats
    vcov <- get.vcov(data, D, Y, Z, X=controls, FE=FE, cl=cl, weights=weights)
    F.standard <- c((t(fs_coefs0) %*% solve(vcov[[1]]) %*% fs_coefs0)/p_iv)
    F.robust <- c((t(fs_coefs0) %*% solve(vcov[[2]]) %*% fs_coefs0)/p_iv)
    if (is.null(cl)==FALSE) {
      F.cluster <- c((t(fs_coefs0) %*% solve(vcov[[3]]) %*% fs_coefs0)/p_iv)
    } else {
      F.cluster <- NA
    }
    F_stat <- c(F.standard, F.robust, F.cluster, F.boot)
    names(F_stat) <- c("F.standard", "F.robust", "F.cluster", "F.boot")
    # calculate ratio
    if (p_iv == 1) {
      ratio <- rep(NA, 1)
      ratio[1] <- IV.Coef/OLS.Coef
    } else {
      ratio <- rep(NA, p_iv + 1)
      names(ratio) <- c("together", Z)
      ratio[1] <- IV.Coef/OLS.Coef
      # run individual IV regressions
      for (i in 1:p_iv) {
        iv.coef.tmp <- IV(data = d0,  D=D, Y=Y, Z=Z[i], X=controls, FE=FE, cl=cl, weights=weights)$coef
        ratio[(i+1)]  <- iv.coef.tmp/OLS.Coef
      }
    }
    # save results
    if (public == TRUE) { # do not publish rho and ratio
      output <- list(
        # OLS and IV results
        est_ols =  round(est_ols, prec),
        est_2sls = round(est_2sls, prec),
        # bootstrap F stat
        F_stat = round(F_stat, prec),
        # number of instruments
        p_iv = p_iv,
        # number of observations
        N = n,
        # number of clusters
        N_cl = ncl
        )
    } else {
      output <- list(
        # OLS and IV results
        est_ols =  round(est_ols, prec),
        est_2sls = round(est_2sls, prec),
        # bootstrap F stat
        F_stat = round(F_stat, prec),
        # number of instruments
        p_iv = p_iv,
        # number of observations
        N = n,
        # number of clusters
        N_cl = ncl,
        # first stage correlation coefficients
        rho_DZ = round(rho, prec),
        # ratio
        ratio = round(ratio, prec),
        # sensitivity analysis
        Bias_threshold_pt_est = OLS.Coef  * (sig_x/sig_epsi) * rho,
        Bias_thresh_ci_lb     = ols_ci[1] * (sig_x/sig_epsi) * rho,
        ## constituent parts of bias computation
        sig_x = sig_x,
        sig_epsi = sig_epsi
        )
    }
    return(output)
}


# %%
