% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/est.zidm.EM.R
\name{est.zidm.EM}
\alias{est.zidm.EM}
\title{Parameter estimation for Zero-inflated Dirichlet-multinomial (ZIDM) distribution.}
\usage{
est.zidm.EM(X, init.a = NULL, init.pi = NULL, iter = 100, conv = 1e-06)
}
\arguments{
\item{X}{a data frame or matrix containing the count data. Rows of the matrix represent observations and columns are the components.}

\item{init.a}{the initial value of DM part for the interative algorithm; if not given, we would like to use moment estimation.}

\item{init.pi}{he initial value of zero-inflated part for the interative algorithm;}

\item{iter}{the number of the interative steps.}

\item{conv}{the precison for stop the estimation.}
}
\value{
a list with the estimated parameters and log-likelihood value.
}
\description{
Using expectation-maximization to estimate the parameter of ZIDM.
}
