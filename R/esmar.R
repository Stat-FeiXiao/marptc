#
#' Title Iteration Process of GEE method
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

es <- function(Time, Status, X,stad,  id, esmax, eps,struc,formula, data) {
  Kn <- length(id)
  K <- length(unique(id))
  n <- as.vector(table(id))
  struc<-struc
  X1 <- X
  if(stad){
    stadX <-  scale(X) #std(X)
  }else{stadX <- X}
  beta=10000
  t11 <- sort(Time)
  c11 <- Status[order(Time)]
  kk <- length(table(t11[c11 == 1]))
  beta2 <-     as.numeric(coxph(formula=formula, data=data)$coef)
  gSSS1<- rep(0,kk)
  Lambda <- initialF(Time, Status,stadX, id,beta2)$Lambda
  repeat{
    SK2<-0
    repeat {

      beta1 <- mbeta(Status, Lambda, X=stadX, beta = beta2,  id,struc)$beta
      if ((any(abs(beta1 - beta2) > eps)) & (SK2 <= esmax) ) {
        beta2 <- beta1
        SK2<-SK2+1
      } else  break
    }

    gSS3 <- baseF(Time, Status, X = stadX, beta = beta1)$gSS3
    gSS1 <-  baseF(Time, Status, X = stadX, beta = beta1)$gSS1
    if ((any(abs(beta1 - beta) > eps))|| any(abs(gSS1-gSSS1)>eps)) {
      beta <- beta1
      gSSS1<-gSS1
      Lambda<-gSS3
    } else  break

  }
  if(stad){
    beta2 <- beta1/attr(stadX,"scaled:scale")
  }else{beta2 <- beta1 }
  gSS3 <- baseF(Time, Status, X = stadX, beta = beta1)$gSS3
  gS <- baseF(Time, Status, stadX , beta = beta1)$gS
 # gSS31 <-     exp(-sum(beta2*attr(stadX,"scaled:center")))*gSS3
  gSS31 <- baseF(Time, Status, X , beta = beta2)$gSS3
  gS1 <- baseF(Time, Status, X , beta = beta2)$gS
  rho <- mbeta(Status, Lambda=gSS31, X, beta = beta2, id,struc)$rho 
  pphi <- mbeta(Status, Lambda=gSS31, X, beta = beta2, id,struc)$pphi

 
  convergence <- sum((beta1 - beta)^2)
  es <- list(beta = beta1, beta2=beta2,gSS3 = gSS3,PgSS3 = gSS31, rho = rho,
             pphi = pphi, gS = gS, PgS = gS1,tau = convergence)
}
