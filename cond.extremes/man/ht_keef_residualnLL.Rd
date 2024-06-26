% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_and_pred.R
\name{ht_keef_residualnLL}
\alias{ht_keef_residualnLL}
\title{computes the negative residual log-likelihood for H&T model under G
normality assumption with keef (2013) constraints}
\usage{
ht_keef_residualnLL(Ymnsi, Yi, theta)
}
\arguments{
\item{Ymnsi}{vector for dependent variable}

\item{Yi}{vector of conditioned variable}

\item{theta}{vector of H&T parameters}
}
\value{
float, negative log-likelihood
}
\description{
computes the negative residual log-likelihood for H&T model under G
normality assumption with keef (2013) constraints
}
\keyword{internal}
