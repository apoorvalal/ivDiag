#' Ensemble of weak instrument tests
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
#' @importFrom glue glue
#' @importFrom fixest feols
#' @importFrom ggplot2 ggplot
#' @export
first_stage_tests <- function(data, Y, D, Z, controls = NULL, FE = NULL, cl = NULL,
    weights = NULL, boot = TRUE, nboots = 1000, parallel = TRUE, seed = 94305, cores = NULL,
    prec = 4) {
    ############################################################
    # data prep
    ############################################################
    # drop missingness
    data <- data[, c(Y, D, Z, controls, FE, cl, weights)]
    d0 <- data[complete.cases(data), ]
    n <- nrow(d0); p_iv <- length(Z)
    ############################################################
    # scatterplot
    ############################################################
    # predict Dhat using first stage
    if(!is.null(weights)){
      dmod = fixest::feols(formula_fixest(D, X = controls, W = Z, D = FE),
        data = df, weights = df[, weights])
    } else{
      dmod = fixest::feols(formula_fixest(D, X = controls, W = Z, D = FE),
        data = df)
    }
    dhat = predict(dmod)
    dtilde = partialer(D, c(Z, controls), FE = FE, weights = weights, data = df)
    scatter = ggplot(data.frame(dtilde, dhat), aes(dtilde, dhat)) +
            geom_point() + geom_smooth() +
            labs(x = "Treatment (residualised)", y = "Predicted Treatment")
    ############################################################
    # fit first stage and store
    out0 <- first_stage_coefs(data = d0, D = D, Z = Z, X = controls, FE = FE, weights = weights)
    rho  <- first_stage_rho(  data = d0, D = D, Z = Z, X = controls, FE = FE, weights = weights)
    AR_res = AR_test(data = d0, Y = Y,   D = D, Z = Z, X = controls, FE = FE)
    fs_coefs0 <- matrix(out0, p_iv, 1)
    if (is.null(cl)==FALSE) { # find clusters
        d0 <- d0[order(d0[,cl]),]
        clusters <- unique(d0[,cl])
        id.list <- split(1:n,d0[,cl])
        ncl <- length(clusters)
    } else {
        ncl <- NULL
    }
    # bootstrap call for first stage coefs
    if (boot){
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
          res <- rep(NA,  p_iv) # last p_iv are first stage IV coefs
          # first stage coefficients
          res <- first_stage_coefs(data=s, D=D, Z=Z, X=controls, FE = FE, weights = weights)
          return(res)
      }
      # in case there's an error
      one.boot <- function() {
        est <- try(boot.core(), silent = TRUE)
        if ('try-error' %in% class(est)) {
          est0 <- rep(NA, p_iv)
          return(est0)
        } else {
          return(est)
        }
      }
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
          expfun <- c("first_stage_coefs", "formula_lfe", "robustify")
          coefs.boot <- foreach (i=1:nboots,.combine=rbind, .inorder=FALSE,
              .export = expfun,
              .packages = c("lfe")) %dopar% {
              return(one.boot())
          }
          doParallel::stopImplicitCluster()
      }
      # post-bootstrap processing
      # drop NAs
      coefs.boot <- coefs.boot[complete.cases(coefs.boot),]
      # OLS and IV CIs
      # Calculate F'
      FStat_bootout <- coefs.boot
      F.boot = c((t(fs_coefs0) %*% solve(var(FStat_bootout)) %*% fs_coefs0)/p_iv)
      names(F.boot) <- "F.boot"
    } else{
      F.boot = NA
    }
    ##############################
    # prep output
    ##############################
    # F.stats
    vcov <- get.vcov(data, D, Y, Z, X=controls, FE=FE, cl=cl, weights=weights)
    F.standard <- c((t(fs_coefs0) %*% solve(vcov[[1]]) %*% fs_coefs0)/p_iv)
    F.robust   <- c((t(fs_coefs0) %*% solve(vcov[[2]]) %*% fs_coefs0)/p_iv)
    if (is.null(cl)==FALSE) {
      F.cluster <- c((t(fs_coefs0) %*% solve(vcov[[3]]) %*% fs_coefs0)/p_iv)
    } else {
      F.cluster <- NA
    }
    F_stat <- c(F.standard, F.robust, F.cluster, F.boot, AR_res$Fstat)
    names(F_stat) <- c("F.standard", "F.robust", "F.cluster", "F.boot", "F.AR")
    output <- list(
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
        # AR Results
        AR.ci_inf = AR_res$ci.info,
        AR.ci = AR_res$ci,
        # scatterplot
        scatterplot = scatter
      )
    return(output)
}
