\name{gbeta}
\alias{gbeta}
\title{Estimate the covariate coefficients using the GEE approach and baseline cumulative distribution function in model using Bernstein polynomial}
\usage{
gbeta(Time, Status, X, stad, Ibeta, id, itermax, eps, N, tau, corstr)
}
\arguments{
\item{Time}{right censored data which is the follow up time.}
\item{Status}{the censoring indicator, normally 0 = event of interest happens, and 0 = censoring.}
\item{X}{a matrix of covariates that may have effect on the failure times.}
\item{stad}{if it is TRUE, all the covariates in the \code{formula} are standardized. By default, \code{stdz = FALSE}.}
\item{Ibeta}{the initial value of the covariate coefficients.}
\item{id}{a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.}
\item{itermax}{specifies the maximum iteration number. If the convergence criterion is not met, the iteration will be stopped after \code{itermax} iterations and the estimates will be based on the last iteration. The default \code{itermax = 100}.}
\item{eps}{tolerance for convergence. The default is \code{eps = 1e-6}. Iteration stops once the relative change in deviance is less than \code{eps}.}
\item{N}{the degree of Bernstein polynomial used to approximate the baseline distribution function in model.}
\item{tau}{the cure threshold. Individuals with survival time greater than \code{tau} are considered cured.}
\item{corstr}{a character string specifying the correlation structure. The following are permitted: \code{independence}, \code{exchangeable} and \code{AR1}.}
}
\description{
Iteration between the estimation of covariate coefficients and baseline cumulative distribution funtion. The iteration will be stoped once the relative change in deviance is less than \code{eps} or the number of iterations reach \code{itermax}.
}
