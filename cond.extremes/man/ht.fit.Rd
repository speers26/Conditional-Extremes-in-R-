% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_and_pred.R
\name{ht.fit}
\alias{ht.fit}
\title{Fitting for conditional extremes}
\usage{
ht.fit(Yi, Ymnsi, qu, Y = T, keef = F, theta0, plot = T)
}
\arguments{
\item{Yi}{vector, data for conditioning variable}

\item{Ymnsi}{vector, data for dependent variable}

\item{Y}{bool, true if we are modeling Y|X, false if X|Y}

\item{keef}{bool, true if we are using Keef (2013) constraints}

\item{theta0}{vector of initial H&T parameter values to be used in optimisation}

\item{plot}{bool, set to true to see visualisation}

\item{q}{float, threshold quantile to fit H&T model above}
}
\value{
list, contains fitted parameters and residuals
}
\description{
fits the H&T (2004) model using Keef's (2013) constraints by optimising
for alpha first, then beta
}
