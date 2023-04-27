\name{AR_test}
\alias{AR_test}
\title{Weak-IV Robust Anderson Rubin Test}
\usage{
AR_test(
  data,
  Y,
  D,
  Z,
  controls,
  FE = NULL,
  cl = NULL,
  weights = NULL,
  prec = 4,
  alpha = 0.05,
  parallel = NULL,
  cores = NULL
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
\item{prec}{precision of CI in string; the default is 4.}
\item{alpha}{level of statitical significance; the default is 0.05.}
\item{parallel}{a logical flag controlling parallel computing.}
\item{cores}{setting the number of cores.}
}
\description{
Performs the Anderson Rubin test, which is robust to weak IVs. 
}
