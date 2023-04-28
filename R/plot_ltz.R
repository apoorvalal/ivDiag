#' Visualise scalar IV coefficient sampling distribution
#' @param out output from ltz
#' @param iv_est a two-element vector of IV estimates (coef & se)
#' @param ltz_est a two-element of local-to-zero estimates (coef & se)
#' @param prior prior mean and se
#' @param xlim
#' @importFrom ggfortify ggdistribution
#' @export
plot_ltz = function(out = NULL, iv_est = NULL, ltz_est = NULL, prior = NULL, xlim = NULL) {

  # check validity
  if (is.null(out)==FALSE) {
    if (inherits(out, "ltz") == FALSE) {stop("\"out\" needs to be a \"ltz\" object.")}
  } else {
    if (is.null(iv_est)==TRUE) {stop("\"iv_est\" is missing")}
    if (is.null(ltz_est)==TRUE) {stop("\"ltz_est\" is missing")}
    if (is.null(prior)==TRUE) {stop("\"mu\" is missing")}
    if (length(iv_est) != 2) {stop("\"iv_est\" should be of length 2")}
    if (length(ltz_est) != 2) {stop("\"iv_est\" should be of length 2")}
    if (length(prior) != 2) {stop("\"prior\" should be of length 2")}
  }

  # extract elements
  if (is.null(out)==FALSE) { # from "out"
    prior <- out$prior
    iv_est <- out$iv[1:2] # coef, se
    ltz_est <- out$ltz[1:2] # coef, se
  } 
  iv_ci <- qnorm(c(0.025, 0.975), iv_est[1], iv_est[2])
  iv_est <- c(iv_est, iv_ci)
  ltz_ci <- qnorm(c(0.025, 0.975), ltz_est[1], ltz_est[2])
  ltz_est <- c(ltz_est, ltz_ci)    
  
  # define common bounds for figure
  if (is.null(xlim) == TRUE) {
    adj <- iv_est[2] * 0.5
    ub = iv_est[1] + 4 * iv_est[2]
    lb = iv_est[1] - 4 * iv_est[2]
    if (ltz_est[3] < lb) lb = ltz_est[3] - adj
    if (ltz_est[4] > ub) ub = ltz_est[4] + adj
    if (lb > 0) lb = - adj # always include 0 in figure
    if (ub < 0) ub = adj # always include 0 in figure
  } else {
    lb <- xlim[1]
    ub <- xlim[2]
  }
  
  ##############################################
  # figures
  ##############################################
  # first: prior on mu
  p0 = ggdistribution(dnorm, seq(lb, ub, 0.01),
    mean = prior[1], sd = prior[2], colour = 'red'
    ) +
    ggplot2::xlim(c(lb, ub)) +
    geom_vline(
      xintercept = prior[1],
      linetype = 'solid', size = 1.5, alpha = 0.5, colour = "red"
    ) +
    geom_vline(
      xintercept = c(prior[1] - 1.96 * prior[2], prior[1] + 1.96 * prior[2]),
      linetype = 'dotted', size = 1.5, colour = "red"
    ) +
    geom_vline(xintercept = 0, alpha = 0.5) +
    labs(title = "Prior on the direct effect", x = "", y = "")  
  # second: original IV
  p1 = ggdistribution(dnorm, seq(lb, ub, 0.01),
    mean = iv_est[1], sd = iv_est[2], colour = 'blue'
    ) +
    ggplot2::xlim(c(lb, ub)) +
    geom_vline(
      xintercept = iv_est[1],
      linetype = 'solid', size = 1.5, alpha = 0.5, colour = "blue"
    ) +
    geom_vline(
      xintercept = c(iv_est[3], iv_est[4]),
      linetype = 'dotted', size = 1.5, colour = "blue"
    ) +
    geom_vline(xintercept = 0, alpha = 0.5) +
    labs(title = "2SLS coefficient", x = "", y = "")  
  # third: local-to-zero
  p2 = ggdistribution(dnorm, seq(lb, ub, 0.01),
    mean = ltz_est[1], sd = ltz_est[2], colour = 'olivedrab4'
    ) +
    ggplot2::xlim(c(lb, ub)) +
    geom_vline(
      xintercept = ltz_est[1],
      linetype = 'solid', size = 1.5, colour = "olivedrab4"
    ) +
    geom_vline(
      xintercept = c(ltz_est[3], ltz_est[4]),
      linetype = 'dotted', size = 1.5, colour = "olivedrab4"
    ) +
    geom_vline(xintercept = 0, alpha = 0.8) +
    labs(title = "2SLS coefficient w/ LTZ adjustment", x = "", y = "")
  # wrapped in patchwork
  p <- patchwork::wrap_plots(list(p0, p1, p2), nrow = 3)
  return(p)
}
