\name{vargee}
\alias{vargee}
\title{Variance estimate with sandwich formula based on the GEE approach}
\usage{
vargee(Time, Status, X, id, beta, baseF, gam, N, corstr)
}
\arguments{
\item{Time}{right censored data which is the follow up time.}
\item{Status}{the censoring indicator, normally 0 = event of interest happens, and 0 = censoring.}
\item{X}{a matrix of covariates that may have effect on the failure times.}
\item{id}{a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.}
\item{beta}{estimation of covariate coefficients based on the GEE approach.}
\item{baseF}{estimation of baseline cumulative distribution function based on the Bernstein polynomial.}
\item{gam}{estimation of the  Bernstein polynomials' coefficients.}
\item{N}{the degree of Bernstein polynomial used to approximate the baseline distribution function in model.}
\item{corstr}{a character string specifying the correlation structure. The following are permitted: \code{independence}, \code{exchangeable} and \code{AR1}.}
}
\description{
Calculate the variance estimates using the sandwich formula based on the GEE approach.
}
