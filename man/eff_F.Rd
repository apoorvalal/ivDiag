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
\references{
  Olea, José Luis Montiel, and Carolin Pflueger. 2013. "A Robust Test for Weak Instruments."" Journal of Business & Economic Statistics 31 (3): 358–69.
}
