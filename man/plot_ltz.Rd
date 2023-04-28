\name{plot_ltz}
\alias{plot_ltz}
\title{Visualizing Local-to-Zero Adjustment}
\description{Visualise approximate sampling distributions for scalar IV coefficient with local-to-zero adjustment.}
\usage{
plot_ltz(out = NULL, iv_est = NULL, ltz_est = NULL, prior = NULL, xlim = NULL)
}
\arguments{
  \item{out}{output from \code{ivDiag::ltz}.}
  \item{iv_est}{a two-element vector of IV estimate and standard error.}
  \item{ltz_est}{a two-element vector of local-to-zero estimate and standard error.}
  \item{prior}{a two-element vector of prior mean and standard deviation.}
  \item{xlim}{a two-element vector specifying the range of the x-axis.}
}
\value{A ggplot2 object.}
\seealso{
  \code{\link{ltz}}
}
\references{
  Conley, Timothy G, Christian B Hansen, and Peter E Rossi. 2012. "Plausibly Exogenous." Review of Economics and Statistics 94 (1): 260–72.
}
