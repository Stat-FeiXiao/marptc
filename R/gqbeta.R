#' Title Estimate the covariate coefficients using the GEE or QIF approach
#'
#' @description Estimate the covariate coefficients using the GEE or QIF approach. For the GEE approach, an estimation of correlation coefficient
#' in the working correlation structure is given.
#'
#' @param Time right censored data which is the follow up time.
#' @param Status the censoring indicator, normally 0 = event of interest happens, and 0 = censoring
#' @param X a matrix of covariates that may have effect on the failure times.
#' @param stad if it is TRUE, all the covariates in the \code{formula} are standardized.
#' @param Ibeta the initial value of the covariate coefficients.
#' @param id a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.
#' @param IC if \code{TRUE}, the function returns the QIC for the GEE method or the BIQIF for the QIF method. Default is \code{FALSE}.
#' @param itermax specifies the maximum iteration number. If the convergence criterion is not met, the iteration will be stopped after \code{itermax} iterations and
#' the estimates will be based on the last iteration.
#' @param eps tolerance for convergence. Iteration stops once the relative change in deviance is less than \code{eps}.
#' @param method using GEE or QIF method to estimate the covariate coefficients.
#' @param corstr a character string specifying the correlation structure. The following are permitted: \code{independence}, \code{exchangeable} and \code{AR1}.
#'
#' @export gqbeta

gqbeta <- function(Time, Status, X, stad, Ibeta, id, IC=FALSE, itermax, eps, method, corstr) {
  
  stadX <- X
  
  if (stad && ncol(X) > 1) {
    cov_X <- X[, -1, drop = FALSE]
    scaled_X <- scale(cov_X)
    mu_X <- attr(scaled_X, "scaled:center")
    sd_X <- attr(scaled_X, "scaled:scale")
    
    stadX[, -1] <- scaled_X
    
    Ibeta_new <- Ibeta
    Ibeta_new[1] <- Ibeta[1] + sum(mu_X * Ibeta[-1])
    Ibeta_new[-1] <- Ibeta[-1] * sd_X
    Ibeta <- Ibeta_new
  }
  
  if (method == "GEE") {
    befit <- geq(Time, Status, X = stadX, beta = Ibeta, id = id, IC=IC,
                 eps = eps, corstr = corstr, itermax = itermax)
  } else if (method == "QIF") {
    befit <- qeq(Time, Status, X = stadX, beta = Ibeta, id = id, IC=IC,
                 eps = eps, corstr = corstr, itermax = itermax)
  } else {
    stop("Invalid method. Please choose 'GEE' or 'QIF'.")
  }
  
  beta1 <- befit$beta
  
  if (stad && ncol(X) > 1) {
    beta2 <- beta1
    beta2[-1] <- beta1[-1] / sd_X
    beta2[1] <- beta1[1] - sum(beta1[-1] * mu_X / sd_X)
  } else {
    beta2 <- beta1
  }
  
  baseF <- befit$baseF
  convergence <- befit$convergence

  if(method=="GEE") {
    rho <- befit$rho
    pphi <- befit$pphi
    V_Ind <- befit$V_Ind
    Q_beta <- befit$Q_beta
    Quad <- NULL
  }else if(method=="QIF"){
    Quad <- befit$Quad
    rho <- pphi <- V_Ind <- Q_beta <- NULL
  }

  return(  list(beta1 = beta1, beta2=beta2, rho = rho, pphi = pphi, Quad = Quad, 
                V_Ind = V_Ind, Q_beta = Q_beta, baseF = baseF, convergence = convergence))
}

#' Title An auxiliary function for solving lambda given beta
#'
#' @description An auxiliary function for solving lambda given beta
#'
#' @param lam unknown parameter corresponding to the largrane multiplier.
#' @param dRbet a matrix needed in this function.
#' @param delta the censoring indicator, normally 0 = event of interest happens, and 0 = censoring.
#'
#' @export getlambda
getlambda <- function(lam, dRbet, delta = NULL) {
  if (!is.null(delta)) {
    dRbet <- dRbet[delta == 1]
  }
  return(sum(1 / (dRbet - lam)) - 1)
}