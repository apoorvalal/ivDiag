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
  \item{Z}{a vector of strings indicating the instrumental variables.}
  \item{controls}{a vector of strings indicating the control variables.}
  \item{FE}{a vector of strings indicating the fixed effects variables.}
  \item{cl}{a string indicating the clustering variable.}
  \item{weights}{a string indicating the variable that stores weights.}
  \item{CI}{a logical flag controlling whether to calcualte the confidence interval using the inversion method.}
  \item{prec}{precision of CI in string; the default is 4.}
  \item{alpha}{level of statitical significance; the default is 0.05.}
  \item{parallel}{a logical flag controlling parallel computing.}
  \item{cores}{setting the number of cores.}
}
\value{
  \item{Fstat}{F statistic, degrees of freedoms, and p-value.}
  \item{ci.print}{Confidence interval via intervsion (printed version).}
  \item{ci}{Confidence interval via intervsion (numeric version).}
  \item{bounded}{If the confidence interval is bounded.}
}
\seealso{
  \code{\link{ivDiag}}
}
\examples{
data(ivDiag)
AR.out <- AR_test(data = rueda, Y = "e_vote_buying", D = "lm_pob_mesa", 
    Z = "lz_pob_mesa_f", controls = c("lpopulation", "lpotencial"),
    cl = "muni_code", CI = FALSE)
library(testthat)    
test_that("Check AR results", {
  expect_equal(as.numeric(AR.out$Fstat[1]), 50.5097)
})
}
\references{
  Chernozhukov, Victor, and Christian Hansen. 2008. "The Reduced Form: A Simple Approach to Inference with Weak Instruments." Economics Letters 100 (1): 68–71.
}