#' Title Estimation of the covariate coefficients using Frailty method
#'
#' @param Status the censoring indicator, 1 = event of interest happens, and 0 = censoring.
#' @param Lambda initial cumulative baseline hazard function.
#' @param X a matrix of covariates that that may have effect on the failure times.
#' @param beta initial values of covariates coefficients.
#' @param varrest the estimation of variance of frailty factors
#' @param id a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.
#'
#' @return estimation of the covariate coefficients and variance of frailty factors.
#' @export


fbeta <- function(Status, Lambda, X, beta,varrest, id) {

  K <- length(unique(id))
  n <- as.vector(table(id))
  c1<-Status
  newY1 <- c1/Lambda
  W1 <- diag(Lambda)
  xxxx<-as.matrix(X)
  betae<-beta
  mu <- exp(xxxx %*% betae )
  dy<-stats::deriv(expression(digamma(x)),"x",func=T)

  de<-function(x) as.numeric(attributes(dy(x))$gradient)

  Pi=1-exp(- mu)
  y1_hat=c1+(1-c1)*( exp(- mu*Lambda) +  Pi-1   )*(exp(- mu*Lambda)) ^(-1)

  B<-rep(0,K)
  D<-rep(0,K)
  for (i in 1:K) {
    B[i]=1/varrest+sum(c1[id==i])

  }

  for (i in 1:K) {
    D[i]=1/varrest+sum(((Lambda*c1+1-y1_hat)* mu)[id==i]  )
  }

  W=B/D

repeat{
  fun<-rep(0,dim(xxxx)[2])
  fu<-rep(0,dim(xxxx)[2])
  for(s in 1:dim(xxxx)[2]){
    fun1=W[id]* mu*xxxx[,s]*(y1_hat-1-c1*Lambda)
    fun2=c1*xxxx[,s]
    fun3=(1-c1)*(exp(-W[id]* mu*Lambda)*(-W[id]* mu*xxxx[,s]*Lambda)+exp(-W[id]* mu)*W[id]* mu*xxxx[,s] )/(exp(-W[id]* mu*Lambda))
    fun[s]=sum(fun1)+sum(fun2)+sum(fun3)
}
  for(s in 1:dim(xxxx)[2]){
    fu1=W[id]* mu*xxxx[,s]*xxxx[,s]*(y1_hat-1-c1*Lambda)
    fu2=(c1-1)*(W[id]* mu*Lambda*xxxx[,s]*xxxx[,s])
    fu3=(1-c1)*(exp(-W[id]* mu*(1-Lambda)))*W[id]* mu*(W[id]* mu*xxxx[,s]*xxxx[,s]*(Lambda-1)+xxxx[,s]*xxxx[,s])
    fu[s]=sum(fu1)+sum(fu2)+sum(fu3)
  }


    oldbeta=betae
    betae=betae-fun/fu
    if( all(abs( betae- oldbeta)<1e-6)) break
    mu <- exp(xxxx %*% betae )

}

  fun2=function(V){
    V^(-2)*(K*digamma(V^(-1))+K*log(V)-K+sum(-(digamma(B)-log(D))+W))
  }
  fun22=function(V){
    V^(-3)*(-2*(K*digamma(1/V)+K*log(V)-K+sum(-(digamma(B)-log(D))+W))+K*de(1/V)*(-1/V)+K)
  }
  repeat{
    oldvarrestt=varrest
    varrest=varrest-fun2(varrest)/fun22(varrest)
    if( abs( varrest- oldvarrestt)<1e-6) break
  }



  list(beta = betae, varrest=varrest,W=W)
}
