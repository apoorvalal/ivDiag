\name{ivDiag}
\alias{ivDiag}
\title{Omnibus function to conduct IV estimation and diagnostics}
\usage{
ivDiag(
  data,
  Y,
  D,
  Z,
  controls = NULL,
  FE = NULL,
  cl = NULL,
  weights = NULL,
  nboots = 1000,
  parallel = TRUE,
  cores = NULL,
  seed = 94305,
  prec = 4,
  debug = FALSE
)
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
\item{seed}{set seed.}
\item{debug}{for debugging purposes.}
}
\value{
list of results
}
\description{
Conduct various estimation and diagnostic procedure for IV designs in one shot.
}
