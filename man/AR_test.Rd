\name{AR_test}
\alias{AR_test}
\title{Anderson Rubin Test}
\description{Performs the Anderson Rubin test, which is robust to weak instruments.}
\usage{
AR_test(data, Y, D, Z, controls, FE = NULL, cl = NULL, weights = NULL, 
  prec = 4, CI = TRUE, alpha = 0.05, parallel = NULL, cores = NULL)
}
\arguments{
  \item{data}{name of a dataframe.}
  \item{Y}{a string indicating the outcome variable.}
  \item{D}{a string indicating the treatment variable.}
  \item{Z}{a vector of strings indicating the instrumental variables. A factor or character instrument enters as dummy variables for all its levels but the first (base) level, named \code{<instrument>_<level>}, so a factor with K levels counts as K - 1 instruments.}
  \item{controls}{a vector of strings indicating the control variables.}
  \item{FE}{a vector of strings indicating the fixed effects variables.}
  \item{cl}{a string indicating the clustering variable.}
  \item{weights}{a string indicating the variable that stores weights.}
  \item{CI}{a logical flag controlling whether to calcualte the confidence interval using the inversion method.}
  \item{prec}{precision of results (4 by default).}
  \item{alpha}{level of statitical significance; the default is 0.05.}
  \item{parallel}{no longer used (the confidence set no longer needs a grid of regressions); kept for backward compatibility.}
  \item{cores}{no longer used; kept for backward compatibility.}
}
\details{
  The confidence set collects every value \eqn{\beta_0} that the AR test does not reject, on the whole real line. The reduced-form coefficients on the instruments are linear in \eqn{\beta_0}, and their robust or cluster-robust variance is quadratic in \eqn{\beta_0}; both are recovered exactly from three regressions and checked at a fourth. With one instrument the set solves a quadratic inequality: a bounded interval, two unbounded rays, the whole real line, or (rarely) empty. With several instruments the AR statistic is evaluated over the whole real line (writing \eqn{\beta_0 = \hat\beta_{2SLS} + SE \cdot \tan\theta}) and every boundary is refined by root finding; the set can then also be a union of several intervals. If the variance is not quadratic in \eqn{\beta_0} (for example, a multiway cluster-robust variance repaired to be positive definite), the regression is run at each point instead.

  Versions up to 1.0.6 inverted the test on a grid of \eqn{\hat\beta_{2SLS} \pm 10} SEs and could miss pieces of the set outside the grid.
}
\value{
  \item{Fstat}{F statistic, degrees of freedoms, and p-value.}
  \item{ci.print}{Confidence interval via intervsion (printed version).}
  \item{ci}{Confidence set via inversion (numeric version): the end points in increasing order, in pairs, with \code{-Inf}/\code{Inf} for unbounded pieces; e.g. \code{c(-Inf, a, b, Inf)} is \eqn{(-\infty, a] \cup [b, \infty)}. \code{NA} if empty.}
  \item{bounded}{If the confidence set is bounded and not empty.}
}
\seealso{
  \code{\link{ivDiag}}
}
\examples{
data(ivDiag)
AR.out <- AR_test(data = mcclendon, Y = "intended", D = "openedesteem", 
    Z = "esteem", controls = NULL)
AR.out$ci.print
}
\references{
  Chernozhukov, Victor, and Christian Hansen. 2008. "The Reduced Form: A Simple Approach to Inference with Weak Instruments." Economics Letters 100 (1): 68–71.
}