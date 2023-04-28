\name{tF}
\alias{tF}
\title{Valid t-Ratio Procedure}
\description{Perform the valid t-ratio procedure.}
\usage{
tF(coef, se, Fstat)
}
\arguments{
\item{coef}{a 2SLS coefficient.}
\item{se}{a standard error estimate for the estimated 2SLS coefficient.}
\item{Fstat}{a first-stage partial F statistic.}
}
\value{
  Results from a valid t-ratio test given the first-stage F statistic.
}
\references{
  Lee, David S, Justin McCrary, Marcelo J Moreira, and Jack Porter. 2022. "Valid t-Ratio Inference for IV." American Economic Review 112 (10): 3260–90.
}