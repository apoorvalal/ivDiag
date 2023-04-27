#' Visualise scalar IV coefficient sampling distribution
#' @param x output from ivDiag
#' @param ltz_est 4-vector of local-to-zero estimates produced by `ivDiag::ltz`
#' @importFrom ggfortify ggdistribution
#' @export
viz_iv_dists = function(out, ltz_est, xlim = NULL) {
  # extract elements from x
  analytic_est <- out$est_2sls[1, 1:2]
  boot_est <- out$est_2sls[2, c(1, 2, 4, 5)] # coef, se, ci_l, ci_u

  # define common bounds for figure
  if (is.null(xlim) == TRUE) {
    adj <- analytic_est[2] * 0.5
    ub = analytic_est[1] + 4 * analytic_est[2]
    lb = analytic_est[1] - 4 * analytic_est[2]
    if (ltz_est[1, 3] < lb) lb = ltz_est[1, 3] - adj
    if (ltz_est[1, 4] > ub) ub = ltz_est[1, 4] + adj
    if (lb > 0) lb = - adj # always include 0 in figure
    if (ub < 0) ub = adj # always include 0 in figure
  } else {
    lb <- xlim[1]
    ub <- xlim[2]
  }
  
  ##############################################
  # figures
  ##############################################
  # first: analytic
  p0 = ggdistribution(dnorm, seq(lb, ub, 0.01),
    mean = analytic_est[1], sd = analytic_est[2], colour = 'red'
    ) +
    ggplot2::xlim(c(lb, ub)) +
    geom_vline(
      xintercept = analytic_est[1],
      linetype = 'solid', size = 1.5, alpha = 0.5, colour = "red"
    ) +
    geom_vline(
      xintercept = c(
        analytic_est[1] - 1.96 * analytic_est[2],
        analytic_est[1] + 1.96 * analytic_est[2]
      ),
      linetype = 'dotted', size = 1.5, colour = "red"
    ) +
    geom_vline(xintercept = 0, alpha = 0.5) +
    labs(title = "2SLS Coef. w/ Analytic CI", x = "", y = "")
  # second: bootstrap-c
  p1 = ggdistribution(dnorm, seq(lb, ub, 0.01),
    mean = boot_est[1], sd = boot_est[2], colour = 'blue'
    ) +
    ggplot2::xlim(c(lb, ub)) +
    geom_vline(
      xintercept = boot_est[1],
      linetype = 'solid', size = 1.5, alpha = 0.5, colour = "blue"
    ) +
    geom_vline(
      xintercept = c(boot_est[3], boot_est[4]),
      linetype = 'dotted', size = 1.5, colour = "blue"
    ) +
    geom_vline(xintercept = 0, alpha = 0.5) +
    labs(title = "2SLS Coef. w/ Bootstrapped CI", x = "", y = "")
  # third: local-to-zero
  p2 = ggdistribution(dnorm, seq(lb, ub, 0.01),
    mean = ltz_est[1, 1], sd = ltz_est[1, 2], colour = 'olivedrab4'
    ) +
    ggplot2::xlim(c(lb, ub)) +
    geom_vline(
      xintercept = ltz_est[1, 1],
      linetype = 'solid', size = 1.5, colour = "olivedrab4"
    ) +
    geom_vline(
      xintercept = c(ltz_est[1, 3], ltz_est[1, 4]),
      linetype = 'dotted', size = 1.5, colour = "olivedrab4"
    ) +
    geom_vline(xintercept = 0, alpha = 0.8) +
    labs(title = "Local-to-Zero Adjustment", x = "", y = "")
  # returns list of ggplot figures - can be plotted separately or wrapped in patchwork
  p <- patchwork::wrap_plots(list(p0, p1, p2), nrow = 3)
  return(p)
}
