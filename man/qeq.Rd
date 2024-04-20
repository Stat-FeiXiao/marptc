\name{geq}
\alias{geq}
\title{Estimation the covariate coefficients using the QIF approach}
\usage{
geq(Status, baseF, X, beta, id, eps, corstr)
}
\arguments{
\item{Status}{right censored data which is the follow up time.}
\item{baseF}{estimation of baseline cumulative distribution function obtained by Bernstein polynomial.}
\item{X}{a matrix of covariates that that may have effect on the failure times.}
\item{beta}{the initial value of the covariate coefficients.}
\item{id}{a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.}
\item{eps}{tolerance for convergence. The default is \code{eps = 1e-6}. Iteration stops once the relative change in deviance is less than \code{eps}.}
\item{corstr}{a character string specifying the correlation structure. The following are permitted: \code{independence}, \code{exchangeable} and \code{AR1}.}
}
\description{
Estimation the covariate coefficients using the QIF approach.
}
