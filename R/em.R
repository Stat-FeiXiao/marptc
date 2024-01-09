#' Title Iteration Process of Frailty method
#'
#' @param Time  right censored data which is the follow up time.
#' @param Status the censoring indicator, 0 = event of interest happens, and 0 = censoring.
#' @param X a matrix of covariates that that may have effect on the failure times.
#' @param id a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.
#' @param emmax specifies the maximum iteration number. If the convergence criterion is not met, the iteration will be stopped after \code{esmax} iterations and the estimates will be based on the last iteration. The default \code{esmax = 100}.
#' @param eps tolerance for convergence. The default is \code{eps = 1e-6}. Iteration stops once the relative change in deviance is less than \code{eps}.
#'
#' @return estimation of the cumulative baseline hazard function, covariate coefficients and parameters in the working matrix.
#' @export

em <- function(Time, Status,X,id, emmax, eps) {
  Kn <- length(id)
  K <- length(unique(id))
  n <- as.vector(table(id))
  SK2<-1
  X1 <- X
  id<-id
  c1<-Status
  xxxx<-X
  varrest<-0.05
  beta2 <- tbeta(Time, Status, X,  id)$tbeta
  Lambda <- initialF(Time, Status, X, id)$gs
  gSS3 <- initialF(Time, Status, X, id)$gSS3
  dy<-stats::deriv(expression(digamma(x)),"x",func=T)

  de<-function(x) as.numeric(attributes(dy(x))$gradient)
  repeat {
      beta1 <- fbeta(Status, Lambda, X=X1, beta = beta2,varrest=varrest,id=id)$beta
      varrest <- fbeta(Status, Lambda, X=X1, beta = beta2,varrest=varrest,id=id)$varrest
      if ((any(abs(beta1 - beta2) > eps)) & (SK2 <= emmax) ) {
        beta2 <- beta1
        SK2<-SK2+1
      } else  break
  }
  W<-fbeta(Status, Lambda, X=X1, beta = beta1,varrest,id)$W
  stat<- (var(W)<0.1)
  if(stat){
    betae <- tbeta(Time, Status, X,  id)$tbeta
    repeat{
      repeat{

        Pi=1-exp(-exp(betae[s]*xxxx))
        y1_hat=c1+(1-c1)*( exp(-exp(betae[s]*xxxx)*gSS3) +  Pi-1   )*(exp(-exp(betae[s]*xxxx)*gSS3)) ^(-1)
        repeat{
        fun<-rep(0,dim(xxxx)[2])
        fun1<-rep(0,dim(xxxx)[2])
        for(s in 1:dim(xxxx)[2]){
          fun1=exp(xxxx%*%betae)*xxxx[,s]*(y1_hat-1-c1*gSS3)
          fun2=c1*xxxx[,s]
          fun3=(1-c1)*(exp(-exp(xxxx%*%betae)*gSS3)*(-exp(xxxx%*%betae)*xxxx[,s]*gSS3)+exp(-exp(xxxx%*%betae))*exp(xxxx%*%betae)*xxxx[,s] )/(exp(-exp(xxxx%*%betae)*gSS3))
          fun[s]=sum(fun1)+sum(fun2)+sum(fun3)
        }
        for(s in 1:dim(xxxx)[2]){
          fu1=exp(xxxx%*%betae)*xxxx[,s]*xxxx[,s]*(y1_hat-1-c1*gSS3)
          fu2=(c1-1)*(exp(xxxx%*%betae)*gSS3*xxxx[,s]*xxxx[,s])
          fu3=(1-c1)*(exp(-exp(xxxx%*%betae)*(1-gSS3)))*exp(xxxx%*%betae)*(exp(xxxx%*%betae)*xxxx[,s]*xxxx[,s]*(gSS3-1)+xxxx[,s]*xxxx[,s])
          fu[s]=sum(fu1)+sum(fu2)+sum(fu3)
        }

        oldbeta=betae
        betae=betae-fun/fu
        if( all(abs( betae- oldbeta)<1e-6)) break
      }
      beta1 <- betae
      if ((any(abs(beta1 - betae) > eps)) & (SK2 <= emmax) ) {
        betae <- beta1
        SK2<-SK2+1
      } else  break
      }
      gSS3<- baseF(Time, Status, X = X1, beta = beta1)$gSS3
      gSS1 <-  baseF(Time, Status, X = X1, beta = beta1)$gSS1
      if ((any(abs(beta1 - beta) > eps))|| any(abs(gSS1-gSSS1)>eps)) {
        beta2 <- beta1
        gSSS1<-gSS1
      } else  break
    }
  }

  gS <- baseF(Time, Status, X=X1 , beta = beta1)$gS
  convergence <- sum((beta1 - beta2)^2)
  list(beta = beta1, varrest=varrest,gSS3 = gSS3,Lambda=Lambda, gS = gS,stat=stat,tau = convergence)
}
