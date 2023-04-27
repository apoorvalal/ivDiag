\name{ltz}
\alias{ltz}
\title{Estimates Local-to-Zero IV coefficients and SEs for a single instrument}
\usage{
ltz(data, Y, D, Z, controls, FE = NULL, cl = NULL, weights = NULL, mu, Sig)
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
\item{mu}{prior point estimate of direct effect of instrument on outcome.}
\item{Sig}{variance of above prior.}
}
\description{
Estimates Local-to-Zero IV coefficients and SEs for a single instrument.
}
