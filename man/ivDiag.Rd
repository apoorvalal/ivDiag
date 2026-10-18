\docType{package}
\name{ivDiag}
\alias{ivDiag}
\title{Omnibus Function for IV Estimation and Diagnostics}
\description{Conducts various estimation and diagnostic procedure for instrumental variable designs in one shot.}
\usage{
ivDiag(data, Y, D, Z = NULL, controls = NULL, FE = NULL, cl = NULL, weights = NULL,
  bootstrap = TRUE, run.AR = TRUE,
  nboots = 1000, parallel = TRUE, cores = NULL, 
  seed = 94305, prec = 4, debug = FALSE,
  running = NULL, cutoff = 0, rd.poly = 1, bandwidth = NULL,
  kernel = "triangular")
}
\arguments{
  \item{data}{name of a dataframe.}
  \item{Y}{a string indicating the outcome variable.}
  \item{D}{a string indicating the treatment variable.}
  \item{Z}{a vector of strings indicating the instrumental variables. May be \code{NULL} when \code{running} is given; the instrument is then the indicator for crossing the cutoff.}
  \item{controls}{a vector of strings indicating the control variables.}
  \item{FE}{a vector of strings indicating the fixed effects variables.}
  \item{cl}{a string indicating the clustering variable.}
  \item{weights}{a string indicating the variable that stores weights.}
  \item{bootstrap}{whether to turn on bootstrap (TRUE by default).}
  \item{run.AR}{whether to run AR test (TRUE by default).}  
  \item{nboots}{a numeric value indicating the number of bootstrap runs.}
  \item{parallel}{a logical flag controlling parallel computing.}
  \item{cores}{setting the number of cores.}
  \item{prec}{precision of CI in string (4 by default).}
  \item{seed}{setting seed.}
  \item{debug}{for debugging purposes.}
  \item{running}{a string indicating the running variable of a fuzzy regression discontinuity (RD) design. When it is given, the instrument is the indicator for crossing the cutoff (built when \code{Z = NULL}), and a polynomial in the distance to the cutoff, separate on each side, is added to \code{controls}, so that only the jump at the cutoff identifies the effect. A treatment predicted from the running variable (for example, from an assignment rule) is not a valid instrument on its own: it also varies smoothly with the running variable, and that variation overstates the first stage. If such a \code{Z} is supplied, it is used with the RD controls and a message is printed.}
  \item{cutoff}{the cutoff, or several cutoffs (each unit is then compared with its nearest cutoff, and cutoff fixed effects are added). 0 by default.}
  \item{rd.poly}{order of the polynomial in the running variable, 0 to 4 (1 by default).}
  \item{bandwidth}{keep units within this distance of the cutoff; \code{NULL} (the whole sample) by default. See \code{rdrobust::rdbwselect} for a data-driven choice.}
  \item{kernel}{\code{"triangular"} (default) or \code{"uniform"}; used only with a bandwidth. Kernel weights multiply \code{weights}.}
}
\value{
  \item{est_ols}{results from an OLS regression.}
  \item{est_2sls}{results from a 2SLS regression.}
  \item{AR}{results from an Anderson-Rubin test}
  \item{F_stat}{various F statistics.}
  \item{rho}{Pearson correlation coefficient between the treatment and predicted treatment from the first stage regression (all covariates are partialled out).}
  \item{tF}{results from the tF procedure based on Lee et al. (2022)}
  \item{est_rf}{results from the reduced form regression.}
  \item{est_fs}{results from the first stage regression.}
  \item{p_iv}{the number of instruments.}
  \item{N}{the number of observations.}
  \item{N_cl}{ the number of clusters.}
  \item{df}{the degree of freedom left from the 2SLS regression}
  \item{nvalues}{the unique values the outcome Y, the treatment D, and each instrument in Z in the 2SLS regression.}
  \item{rd}{with \code{running} only: the RD setup (running variable, cutoffs, polynomial order, bandwidth, kernel) and the number of units on each side of the cutoff.}
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
g <- ivDiag(data = mcclendon, Y = "intended", D = "openedesteem", 
    Z = "esteem", bootstrap = FALSE)
plot_coef(g)
library(testthat)    
test_that("Check ivDiag output", {
  expect_equal(as.numeric(g$est_2sls[1,1]), 0.3149)
})
}


