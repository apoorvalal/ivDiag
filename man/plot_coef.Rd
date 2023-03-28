\name{plot_coef}
\alias{plot_coef}
\title{Visualise point estimates and CIs of OLS and IV estimates}
\usage{
plot_coef(out)
}
\arguments{
\item{out}{output from \code{ivDiag}}
\item{ols.methods}{a vector specifying the inferential methods for OLS to be shown. The default is \code{c("analy","bootc","boott")}}
\item{iv.methods}{a vector specifying inferential methods for 2SLS to be shown. The default is \code{c("analy","bootc","boott","ar","tf")}}
\item{main}{a string specifying the title of the plot}
\item{ylab}{a string specifying ylab of the plot}
\item{grid}{a logical flag indicating whether to show the grids}
\item{stats}{a logical flag indicating  whether to show the statistics, including the effective F, the number of observations, and the number of clusters (if applicable)}
\item{ylim}{a two-element vector specifying the range of the y-axis}
}
\description{
Visualise point estimates and confidence intervals of OLS and IV estimates
}
