% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_and_pred.R
\name{ht.par.boot}
\alias{ht.par.boot}
\title{Bootstrapped HT parameter estimation}
\usage{
ht.par.boot(data, q, n.boot)
}
\arguments{
\item{data}{2 column matrix, data to fit initial model to}

\item{q}{threshold exceedance quantile to fit intial model above}

\item{n.boot}{number of parametric bootstraps to perform}
}
\value{
list of parametric bootstrapped theta and residuals
}
\description{
performs parametric bootstrap for ht (2004) model with keef (2013) constraints
}
