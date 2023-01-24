#' Visualise scalar IV coefficient sampling distribution
#' @param analytic_est 2-vector of IV coefficient and SE
#' @param boot_est k-vector of bootstrap estimates of the IV coefficient
#' @param ltz_est 4-vector of local-to-zero estimates produced by `ivDiag::ltz`
#' @importFrom ggfortify ggdistribution
#' @export
viz_iv_dists = function(analytic_est, boot_est, ltz_est) {
  # define common bounds for figure
  ub = analytic_est[1] + 3 * analytic_est[2]
  lb = analytic_est[1] - 3 * analytic_est[2]
  if (min(boot_est) < lb) lb = min(boot_est) - 1
  if (max(boot_est) > ub) ub = max(boot_est) + 1
  if (ltz_est$ci[1] < lb) lb = ltz_est$ci[2] - 1
  if (ltz_est$ci[2] > ub) ub = ltz_est$ci[1] + 1
  if (lb > 0) lb = -1 # always include 0 in figure
  lb = floor(lb); ub = ceiling(ub) # round x axis values
  ##############################################
  # figures
  ##############################################
  p0 = ggdistribution(dnorm, seq(lb, ub, 0.01),
    mean = analytic_est[1], sd = analytic_est[2], colour = 'red'
  ) +
    xlim(c(lb, ub)) +
    geom_vline(
      xintercept = analytic_est[1],
      linetype = 'solid', size = 1.5, alpha = 0.5, colour = "red"
    ) +
    geom_vline(
      xintercept = c(
        analytic_est[1] - 1.96 * analytic_est[2],
        analytic_est[1] + 1.96 * analytic_est[2]
      ),
      linetype = 'dashed', size = 1.5, colour = "red"
    ) +
    geom_vline(xintercept = 0, alpha = 0.5) +
    labs(title = "Conventional 2SLS", x = "", y = "")
  # second - bootstrap
  p1 = ggplot(data.frame(ivb = boot_est), aes(x = ivb)) +
    geom_density(colour = 'blue') +
    xlim(c(lb, ub)) +
    geom_vline(
      xintercept = quantile(boot_est, c(0.025, 0.975)),
      size = 1.5, linetype = 'dotted', colour = 'blue'
    ) +
    geom_vline(xintercept = quantile(boot_est, c(0.5)), size = 1.5, colour = 'blue') +
    geom_vline(xintercept = 0, alpha = 0.5) +
    labs(title = "Bootstrap", x = "", y = "")
  # third - local-to-zero
  p2 = ggdistribution(dnorm, seq(lb, ub, 0.01),
    mean = ltz_est$est, sd = ltz_est$se, colour = 'olivedrab4'
  ) +
    xlim(c(lb, ub)) +
    geom_vline(
      xintercept = ltz_est$est[1, 1],
      linetype = 'dashed', size = 1.5, colour = "olivedrab4"
    ) +
    geom_vline(
      xintercept = c(ltz_est$ci[1], ltz_est$ci[2]),
      linetype = 'dotted', size = 1.5, colour = "olivedrab4"
    ) +
    geom_vline(xintercept = 0, alpha = 0.8) +
    labs(title = "Local-to-zero", x = "", y = "")
  # returns list of ggplot figures - can be plotted separately or wrapped in patchwork
  return(list(p0, p1, p2))
}
