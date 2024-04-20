\name{geeptc}
\alias{geeptc}
\title{Fit the marginal semiparametric promotion time cure model using GEE method}
\usage{
geeptc(
  formula,
  data,
  id,
  Var = TRUE,
  Ibeta = NULL,
  boots = FALSE,
  stad = FALSE,
  nboot = 100,
  N,
  tau = NULL,
  corstr = "independence",
  itermax = 100,
  eps = 1e-06
)
}
\arguments{
\item{formula}{a formula expression, of the form \code{response ~ predictors}. The \code{response} is a \code{Surv} object with right censoring. The expression to the right of the "~" specifies the effect of covariates on the failure time.}
\item{data}{a data frame in which to interpret the variables named in the \code{formula}.}
\item{id}{a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.}
\item{Var}{if it is TRUE, the program returns Std.Error. By default, \code{Var = TRUE}.}
\item{Ibeta}{the initial value of the covariate coefficients.}
\item{boots}{if it is TRUE, the program returns Std.Error by the bootstrap method. By default, \code{boots = FALSE}}
\item{stad}{if it is TRUE, all the covariates in the \code{formula} are standardized. By default, \code{stdz = FALSE}.}
\item{nboot}{the number of bootstrap samples. The default is \code{nboot = 100}.}
\item{N}{the degree of Bernstein polynomial used to approximate the baseline distribution function in model.}
\item{tau}{the cure threshold. Individuals with survival time greater than \code{tau} are considered cured.}
\item{corstr}{a character string specifying the correlation structure. The following are permitted: \code{independence}, \code{exchangeable} and \code{AR1}.}
\item{itermax}{specifies the maximum iteration number. If the convergence criterion is not met, the iteration will be stopped after \code{esmax} iterations and the estimates will be based on the last iteration. The default \code{itermax = 100}.}
\item{eps}{tolerance for convergence. The default is \code{eps = 1e-6}. Iteration stops once the relative change in deviance is less than \code{eps}.}
}
\description{
Fit marginal promotion time cure (PTC) model based on the general estimation equations (GEE) method. We consider three common correlation structures in this funciton. The baseline cumulative distribution function in PTC model is appromximated by the Bernstein polynomial.
}
