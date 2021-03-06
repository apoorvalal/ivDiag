% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/viz_iv_dists.R
\name{viz_iv_dists}
\alias{viz_iv_dists}
\title{Visualise scalar IV coefficient sampling distribution}
\usage{
viz_iv_dists(analytic_est, boot_est, ltz_est)
}
\arguments{
\item{analytic_est}{2-vector of IV coefficient and SE}

\item{boot_est}{k-vector of bootstrap estimates of the IV coefficient}

\item{ltz_est}{4-vector of local-to-zero estimates produced by \code{ivDiag::ltz}}
}
\description{
Visualise scalar IV coefficient sampling distribution
}
