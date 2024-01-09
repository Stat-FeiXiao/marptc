#' Title Iteration Process of GEE method improved by quadratic inference functions
#'
#' @param Time  right censored data which is the follow up time.
#' @param Status the censoring indicator, 0 = event of interest happens, and 0 = censoring.
#' @param X a matrix of covariates that that may have effect on the failure times.
#' @param id a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.
#' @param esmax specifies the maximum iteration number. If the convergence criterion is not met, the iteration will be stopped after \code{esmax} iterations and the estimates will be based on the last iteration. The default \code{esmax = 100}.
#' @param eps tolerance for convergence. The default is \code{eps = 1e-6}. Iteration stops once the relative change in deviance is less than \code{eps}.
#'
#' @return estimation of the cumulative baseline hazard function, covariate coefficients and parameters in the working matrix.
#' @export

esqif <- function(Time, Status, X,stad,  id, esmax, eps,formula, data) {
  Kn <- length(id)
  K <- length(unique(id))
  n <- as.vector(table(id))
  t11 <- sort(Time)
  c11 <- Status[order(Time)]
  kk <- length(table(t11[c11 == 1]))
  SK1<-1
  X1 <- X
  if(stad){
    stadX <-  std(X) #std(X)
    beta2 <- as.numeric(coxph(formula=formula, data=data)$coef) 
    Lambda <- initialF(Time, Status, stadX, id,beta2)$Lambda
  }else{
    stadX <- X
  beta2 <- as.numeric(coxph(formula=formula, data=data)$coef) 
  Lambda <- initialF(Time, Status, stadX, id,beta2)$Lambda
  } 
  beta=rep(10000,dim(X)[2])
 # beta2 <- tbeta(Time, Status, X,  id)$tbeta
 # beta2 <- coxph(formula, data,data)
 

  gSSS1 <- rep(0,kk)
  repeat{
    SK2<-0
    repeat {

      beta1 <- qbeta(Status, Lambda, X=stadX, beta = beta2,  id,eps)$beta
 
      if ((any(abs(beta1 - beta2) > eps)) & (SK2 <= esmax) ) {
       beta2 <- beta1
       SK2<-SK2+1
      } else  break
    }

    gSS3 <- baseF(Time, Status, X = stadX, beta = beta1)$gSS3
    gSS1 <-  baseF(Time, Status, X = stadX, beta = beta1)$gSS1
    if (((any(abs(beta1 - beta) > eps))) & (SK1 <= esmax)) {    #|| any(abs(gSS1-gSSS1)>eps) 
      beta <- beta1
      gSSS1<-gSS1
      Lambda<-gSS3
      SK1<-SK1+1
    } else  break

  }
  if(stad){
    beta2 <- beta1/attr(stadX,"scale")
  }else{ beta2 <- beta1 }


  gS <- baseF(Time, Status, X=stadX , beta = beta1)$gS
  convergence <- sum((beta1 - beta)^2)
  list(beta = beta1,beta2 = beta2, gSS3 = gSS3, gS = gS, tau = convergence)
}
