\name{tF}
\alias{tF}
\title{Valid t-Ratio Procedure}
\description{Performs the valid t-ratio procedure.}
\usage{
tF(coef, se, Fstat, prec = 4)
}
\arguments{
\item{coef}{a 2SLS coefficient.}
\item{se}{a standard error estimate for the estimated 2SLS coefficient.}
\item{Fstat}{a first-stage partial F statistic.}
\item{prec}{precision of results (4 by default).}
}
\details{
  The critical value \code{cF} is interpolated from the table in Lee et al. (2022) by the
  square root of \code{Fstat}. It is 1.96 for \code{Fstat} of 106.09 (10.3 squared) and above.
  Below \code{Fstat} = 4, where the table starts, \code{cF} is \code{Inf}: the valid critical
  value is infinite for \code{Fstat} up to 3.84 (1.96 squared) and above 18.66 between 3.84
  and 4. The confidence interval is then the whole real line and the p-value is 1.
}
\value{
  A named vector: the F statistic (\code{F}), the critical value (\code{cF}), the coefficient,
  its SE, the t statistic, the confidence interval (\code{Coef} plus or minus \code{cF} times
  \code{SE}) and the adjusted p-value.
}
\examples{
tF(coef = 0.3149, se = 0.0890, Fstat = 207.18)
}
\references{
  Lee, David S, Justin McCrary, Marcelo J Moreira, and Jack Porter. 2022. "Valid t-Ratio Inference for IV." American Economic Review 112 (10): 3260–90.
}