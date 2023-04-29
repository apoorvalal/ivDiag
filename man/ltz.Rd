\name{ltz}
\alias{ltz}
\title{Local-to-Zero Test}
\description{Estimates Local-to-Zero IV coefficients and SEs for a single instrument.}
\usage{
ltz(data, Y, D, Z, controls, FE = NULL, cl = NULL, weights = NULL, prior, prec = 4)
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
  \item{prior}{prior mean and standard deviation of the direct effect of instrument on outcome.}
  \item{prec}{precision of CI in string; the default is 4.}
}
\value{
  \item{iv}{results from a 2SLS regression.}
  \item{ltz}{results after local-to-zerio adjustment.}
  \item{prior}{prior mean and standard deviation}  
}
\examples{
data(ivDiag)
controls <- c('altitudine', 'escursione', 'costal', 'nearsea', 'population', 
    'pop2', 'gini_land', 'gini_income')
ltz_out <- ltz(data = gsz, Y = "totassoc_p", D = "libero_comune_allnord", 
    Z = "bishopcity", controls = controls, weights = "population", 
    prior = c(0.178, 0.137))
library(testthat)    
test_that("Check local-to-zero adjustment", {
  expect_equal(as.numeric(ltz_out$ltz[1]), 3.6088)
})
}
\references{
  Conley, Timothy G, Christian B Hansen, and Peter E Rossi. 2012. "Plausibly Exogenous." Review of Economics and Statistics 94 (1): 260–72.
}