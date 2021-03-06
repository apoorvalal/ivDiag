% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/first_stage_tests.R
\name{first_stage_tests}
\alias{first_stage_tests}
\title{Ensemble of weak instrument tests}
\usage{
first_stage_tests(
  data,
  Y,
  D,
  Z,
  controls = NULL,
  FE = NULL,
  cl = NULL,
  weights = NULL,
  boot = TRUE,
  nboots = 1000,
  parallel = TRUE,
  seed = 94305,
  cores = NULL,
  prec = 4,
  scatter = F
)
}
\arguments{
\item{data}{dataframe}

\item{Y}{outcome (string)}

\item{D}{treatment (string)}

\item{Z}{instrument (string)}

\item{controls}{control variables (character vector)}

\item{FE}{fixed effects (character vector)}

\item{cl}{clustering column for SE (character vector)}

\item{weights}{weights name (string)}

\item{nboots}{number of bootstrap reps (1000 by default)}

\item{parallel}{boolean for parallel bootstrap (on by default)}

\item{seed}{seed}

\item{cores}{number of cores to parallelise across (defaults to all - 1)}

\item{prec}{precision of summary (4 by default)}

\item{scatter}{scatterplot of treatment against predicted treatment (to diagnose functional form)}
}
\value{
list of results
}
\description{
Ensemble of weak instrument tests
}
