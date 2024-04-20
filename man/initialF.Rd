\name{initialF}
\alias{initialF}
\title{Initial value of the baseline cumulative distribution function in model}
\usage{
initialF(Time, Status, X, id, betaest)
}
\arguments{
\item{Time}{right censored data which is the follow up time.}
\item{Status}{the censoring indicator, normally 0 = event of interest happens, and 0 = censoring}
\item{X}{a matrix of covariates that may have effect on the failure times.}
\item{id}{a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.}
}
\description{
Obtain the initial value of the baseline cumulative distribution function in model.
}
