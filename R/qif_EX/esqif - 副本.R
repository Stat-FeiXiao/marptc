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
 stadX <- X
  if(stad){
   for (i in 2:ncol(X1)) {
                stadX[, i] <- (X[, i] - mean(X[, i]))/sd(X[, i])
            }
  #  stadX <-  std(X) #std(X)
    beta2 <-  rep(0,dim(X)[2])#rep(0,dim(X)[2]) #as.numeric(coxph(formula=formula, data=data)$coef) 
    Lambda <- initialF(Time, Status, stadX, id,beta2)$Lambda
  }else{
    stadX <- X      
    beta2 <- rep(0,dim(X)[2])#rep(0,dim(X)[2]) # as.numeric(coxph(formula=formula, data=data)$coef) 
    Lambda <- initialF(Time, Status, stadX, id,beta2)$Lambda
  } 
Lambda1 <- Lambda
  beta=10000
  t11 <- sort(Time)
  c11 <- Status[order(Time)]
  kk <- length(table(t11[c11 == 1]))
  gSSS1<- rep(0,kk)
  alpha_esti <- rep( 1/Kn,kk-1)
POOL1 <- POOL2<- list()
Qic <- rep(1,6)
 for(ii in 2:7){
N <- ii
Lambda <- Lambda1
 beta2 <- rep(0,dim(X)[2])
  repeat{
    SK2<-0
    repeat {

 QQQ <- qbeta(Status, Lambda, X=stadX, beta = beta2,  id,eps)
      beta1 <- QQQ$beta
    
 
      if ((any(abs(beta1 - beta2) > eps)) & (SK2 <= esmax) ) {
       beta2 <- beta1
       SK2<-SK2+1
      } else  break
    }
BBB<- baseF(Time, Status, X = stadX, beta = beta1,N)
  gSS3 <- BBB$gSS3
    if (((any(abs(beta1 - beta) > eps))) & (SK1 <= esmax)) {    #|| any(abs(gSS1-gSSS1)>eps) 
      beta <- beta1
      Lambda<-gSS3
    } else  break
  }
Qic[ii-1] <-BBB$qic
POOL1[[ii-1]] <- BBB
POOL2[[ii-1]] <- QQQ
}
BBB <- POOL1[[which.min(Qic)]]
QQQ<- POOL2[[which.min(Qic)]]
 gSS3 <- BBB$gSS3
  gS <- BBB$gS
  gSS31 <- NULL #    exp(-sum(beta2*attr(stadX,"scaled:center")))*gSS3
  gSS31 <- NULL #BBB$gSS3
  gS1 <- NULL #  BBB$gS
beta1 <- QQQ$beta


   if(stad){
 beta2<- c(beta1[1] - sum((beta1[-1] * apply(X1[, -1, drop = FALSE], 2, mean)/apply(X1[, -1, drop = FALSE], 2, sd))), beta1[-1]/apply(X1[,  -1, drop = FALSE], 2, sd)) 
  }else{ beta2 <- beta1 }


  convergence <- sum((beta1 - beta)^2)
  list(beta = beta1,beta2 = beta2, gSS3 = gSS3, gS = gS, tau = convergence)
}
