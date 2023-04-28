#' @param out an ivDiag object
#' @param ols.methods inferential methods for OLS
#' @param iv.methods inferential methods for 2SLS
#' @param main title of the plot
#' @param ylab ylab of the plot
#' @param grid whether to put on grids
#' @param stats whether to show Statistics
#' @param ylim y-range
#' @export
#' 
plot_coef <- function(out, 
  ols.methods = c("analy","bootc","boott"),
  iv.methods = c("analy","bootc","boott","ar","tf"),
  main = NULL,
  ylab = "Coefficient",
  grid = TRUE,
  stats = TRUE,
  ylim = NULL
) {

  # check
  if (is.null(out)==FALSE) {
    if (inherits(out, "ivDiag") == FALSE) {stop("\"out\" needs to be a \"ltz\" object.")}
  }

  # title
  if (is.null(main)==TRUE) {
    main <- "OLS and 2SLS Estimates with 95% CIs"
  }

  # stats
  for (s in ols.methods) {
    if (!s %in% ols.methods) {stop("\"ols.methods\" should only include \"analy\", \"bootc\", \"boott\".")}
  }
  for (s in iv.methods) {
    if (!s %in% iv.methods) {stop("\"iv.methods\" should only include \"analy\", \"bootc\", \"boott\", \"ar\", \"tf\".")}
  }
  
  # collect coefficients and CIs
  coef.ols <- out$est_ols[1,1]
  coef.iv <- out$est_2sls[1,1]
  coef <- c(rep(coef.ols, 3), rep(coef.iv, 4))
  lower_ci <- c(out$est_ols[,4], out$est_2sls[, 4], out$AR$ci[1])
  upper_ci <- c(out$est_ols[,5], out$est_2sls[, 5], out$AR$ci[2])
  method <- c("Analytic", "Boot-c", "Boot-t", "Analytic", "Boot-c", "Boot-t", "AR")
  row.names <- c("ols-analy", "ols-bootc", "ols-boott", "iv-analy", "iv-bootc", "iv-boott", "iv-ar")

  # is AR bounded?
  if (out$AR$bounded == FALSE) {
    lower_ci[7] <- upper_ci[7] <- NA
  }

  # is tF available?
  if (is.null(out$tF) == FALSE) {
    coef <- c(coef, coef.iv)
    lower_ci <- c(lower_ci, out$tF[6])
    upper_ci <- c(upper_ci, out$tF[7])
    method <- c(method, "tF")
    row.names <- c(row.names, "iv-tf")
  } else {
    if ("tf" %in% iv.methods) {
      iv.methods <- setdiff(iv.methods, "tf")
    }
  }

  # update stats names
  ols.methods2 <- paste0("ols-", ols.methods)
  iv.methods2 <- paste0("iv-", iv.methods)

  # Data for the plot
  data <- data.frame(method = method, coef = coef, lower_ci = lower_ci, upper_ci = upper_ci)
  rownames(data) <- row.names
  data <- data[c(ols.methods2, iv.methods2), ]
 
  rg <- range(data[,c(3,4)], na.rm = TRUE)
  adj <- rg[2] - rg[1]
  if (is.null(ylim) == TRUE) {
    ylim  <- c(min(0, rg[1] - 0.3*adj), max(0, rg[2] + 0.35*adj))
  }
  adj2 <- ylim[2] - ylim[1] 
  
  # Set up the plot
  ncoefs.ols <- length(ols.methods2)
  ncoefs.iv <- length(iv.methods2)
  ncoefs <- ncoefs.ols + ncoefs.iv
  plot(1: ncoefs, data$coef, xlim = c(0.5, ncoefs + 0.5), ylim = ylim,
      ylab = "", xlab = "", main = main, 
      axes = FALSE, xaxt = "n", yaxt = "n", type = "n")
  axis(1, at = 1: ncoefs, labels = data$method, las = 2, cex.axis = 0.8)
  axis(2)
  mtext(ylab, 2, line = 2.5)
  if (grid == TRUE) {
    abline(h = axTicks(2), lty = "dotted", col = "gray50")
    abline(v = c(0.5, c(1: ncoefs) + 0.5), lty = "dotted", col = "gray50") # horizontal grid
  }
  if (ncoefs.ols > 0) {
    rect(0.5, ylim[1]-0.2*adj2, ncoefs.ols + 0.5, ylim[1]+0.05*adj2, col = "#AAAAAA80", border = "white") # label at bottom
    text((ncoefs.ols + 1)/2, ylim[1], "OLS", cex = 1)
  }
  if (ncoefs.iv > 0) {  
    rect(ncoefs.ols + 0.5, ylim[1]-0.2*adj2, ncoefs + 0.5, ylim[1]+0.05*adj2, col = "#77777780", border = "white") # label at bottom
    text((ncoefs.iv + 1)/2 + ncoefs.ols, ylim[1], "2SLS", cex = 1)
  } 
  abline(h = 0, col = "red", lwd = 2, lty = "solid")
  segments(y0 = data$lower_ci, x0 = c(1: ncoefs), y1 = data$upper_ci, x1 = c(1: ncoefs), lwd = 2) #CI
  points(1: ncoefs, data$coef, pch = c(rep(16, ncoefs.ols), rep(17, ncoefs.iv)), col = "blue") #point coefs
  if ("ar" %in% iv.methods && out$AR$bounded == FALSE) {
    if (length(out$AR$ci) == 1) {
      ar.lab <- "Empty CI"
    } else {
      ar.lab <- "Unbounded CI"
    }
    text(which(data$method == "AR"),  -0.05 * adj2 * sign(coef.iv), ar.lab, col = "red", cex = 0.8)
  }
  if (stats == TRUE) {
    if (is.null(out$N_cl) == TRUE) {
      legend("topleft", c(paste0("Effective F = ", sprintf("%.1f", out$F_stat[5])),
        paste0("N = ", out$N)), bty="n")
    } else {
      legend("topleft", c(paste0("Effective F = ", sprintf("%.1f", out$F_stat[5])),
        paste0("N = ", out$N, ", Ncl = ", out$N_cl)), bty="n")
    }
    
  }    
  box()
}



