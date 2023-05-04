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
\value{
  Results from a valid t-ratio test given the first-stage F statistic.
}
\examples{
tf.out <- tF(coef = -0.9835, se = 0.1540, Fstat = 8598)
library(testthat)
test_that("Check tF cF", {
  expect_equal(as.numeric(tf.out[2]), 1.96)
})
}
\references{
  Lee, David S, Justin McCrary, Marcelo J Moreira, and Jack Porter. 2022. "Valid t-Ratio Inference for IV." American Economic Review 112 (10): 3260–90.
}