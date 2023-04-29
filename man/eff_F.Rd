\name{eff_F}
\alias{eff_F}
\title{Effective F}
\description{Computes the effective F statistic.}
\usage{
eff_F(data, Y, D, Z, controls = NULL, FE = NULL, cl = NULL, weights = NULL)
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
}
\value{
  the effective F statistic. 
}
\seealso{
  \code{\link{ivDiag}}
}
\examples{
effF <- eff_F(data = rueda, Y = "e_vote_buying", D = "lm_pob_mesa", 
    Z = "lz_pob_mesa_f", controls = c("lpopulation", "lpotencial"), 
    cl = "muni_code")
library(testthat)
test_that("Check effective F", {
  expect_equal(floor(as.numeric(effF)), 8598)
})
}
\references{
  Olea, José Luis Montiel, and Carolin Pflueger. 2013. "A Robust Test for Weak Instruments."" Journal of Business & Economic Statistics 31 (3): 358–69.
}
