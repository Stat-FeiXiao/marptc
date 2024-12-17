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
#' @param itermax specifies the maximum iteration number. If the convergence criterion is not met, the iteration will be stopped after \code{itermax} iterations and
#' the estimates will be based on the last iteration.
#' @param eps tolerance for convergence. Iteration stops once the relative change in deviance is less than \code{eps}.
#' @param method using GEE or QIF method to estimate the covariate coefficients.
#' @param corstr a character string specifying the correlation structure. The following are permitted: \code{independence}, \code{exchangeable} and \code{AR1}.
#'
#' @export gqbeta

gqbeta <- function(Time, Status, X, stad, Ibeta, id, itermax, eps, method, corstr) {

  Kn <- length(id)
  K <- length(unique(id))
  n <- as.vector(table(id))
  stadX <- X

  if(stad){
    Ibeta1 <- rep(0,length(Ibeta))
    cumIbe <- 0
    for(lln in 2:length(Ibeta)){
      cumIbe <- mean(X[,lln]) * Ibeta[lln]
      Ibeta1[lln] =  sd(X[,lln ]) * Ibeta[lln]
    }
    Ibeta1[1] <- Ibeta[1] + cumIbe
    Ibeta <- Ibeta1

   for (i in 2:ncol(X)) {
                stadX[, i] <- (X[, i] - mean(X[, i]))/sd(X[, i])
            }
  }

  beta2 <- Ibeta


  iter <- 1

     repeat{
       if(method=="GEE"){
         befit <- geq(Time,Status, X=stadX, beta = beta2,  id, eps, corstr, itermax)
       }else if(method=="QIF")  befit <- qeq(Time,Status, X=stadX, beta = beta2,  id, eps, corstr, itermax)
      beta1 <- befit$beta
      if ( (any(abs(beta1 - beta2) > eps))   & (iter<= itermax)) {
        beta2 <- beta1
        iter <- iter + 1
      } else  break
     }

  if(stad){
    beta2<- c(beta1[1] - sum((beta1[-1] * apply(X[, -1, drop = FALSE], 2, mean)/apply(X[, -1, drop = FALSE], 2, sd))),
              beta1[-1]/apply(X[,  -1, drop = FALSE], 2, sd))
  }else{ beta2 <- beta1 }

  if(method=="GEE") {
    rho <- befit$rho
  }else rho <- NULL

   convergence <- befit$convergence & (iter<= itermax)


  list(beta1 = beta1, beta2=beta2,rho = rho, convergence = convergence)
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

getlambda <- function(lam,dRbet,delta){

  val <- sum(1/(dRbet[delta==1]-lam)) - 1
  return(val)

}
