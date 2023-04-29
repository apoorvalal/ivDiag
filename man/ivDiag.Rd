\docType{package}
\name{ivDiag}
\alias{ivDiag}
\title{Omnibus Function for IV Estimation and Diagnostics}
\description{Conducts various estimation and diagnostic procedure for instrumental variable designs in one shot.}
\usage{
ivDiag(data, Y, D, Z, controls = NULL, FE = NULL, cl = NULL, weights = NULL,
  nboots = 1000, parallel = TRUE, cores = NULL, 
  seed = 94305, prec = 4, debug = FALSE)
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
  \item{nboots}{a numeric value indicating the number of bootstrap runs.}
  \item{parallel}{a logical flag controlling parallel computing.}
  \item{cores}{setting the number of cores.}
  \item{prec}{precision of CI in string; the default is 4.}
  \item{seed}{setting seed.}
  \item{debug}{for debugging purposes.}
}
\value{
  \item{est_ols}{results from an OLS regression.}
  \item{est_2sls}{results from a 2SLS regression.}
  \item{AR}{results from an Anderson-Rubin test}
  \item{F_stat}{various F statistics.}
  \item{rho}{Pearson correlation coefficient between the treatment and predicted treatment from the first stage regression (all covariates are partialled out).}
  \item{tF}{results from the tF procedure based on Lee et al. (2022)}
  \item{est_rf}{results from the first stage regression.}
  \item{est_fs}{results from the reduced form regression.}
  \item{p_iv}{the number of instruments.}
  \item{N}{the number of observations.}
  \item{N_cl}{ the number of clusters.}
  \item{df}{the degree of freedom left from the 2SLS regression}
  \item{nvalues}{the unique values the outcome Y, the treatment D, and each instrument in Z in the 2SLS regression.}
}
\author{
  Apoorva Lal; Yiqing Xu
}
\references{
Lal, Apoorva, Mackenzie William Lockhart, Yiqing Xu, and Ziwen Zu. 2023. "How Much Should We Trust Instrumental Variable Estimates in Political Science? Practical Advice Based on 67 Replicated Studies." Available at: \url{https://yiqingxu.org/papers/english/2021_iv/LLXZ.pdf}

Lee, David S, Justin McCrary, Marcelo J Moreira, and Jack Porter. 2022. "Valid t-Ratio Inference for IV." American Economic Review 112 (10): 3260–90.
}
\seealso{
  \code{\link{plot_coef}}
  \code{\link{eff_F}}
  \code{\link{AR_test}}
  \code{\link{tF}}
}
\examples{
data(ivDiag)
g <- ivDiag(data = rueda, Y = "e_vote_buying", D = "lm_pob_mesa", 
    Z = "lz_pob_mesa_f", controls = c("lpopulation", "lpotencial"), 
    cl = "muni_code", parallel = FALSE, nboots = 10)
plot_coef(g)    
    
library(testthat)
test_that("Check omnibus function", {
  expect_equal(as.numeric(g$est_2sls[1, 1]), -0.9835)
  expect_equal(as.numeric(g$AR$Fstat[1]), 50.5097)
})
}


