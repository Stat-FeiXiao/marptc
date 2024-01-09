#' Title Estimation of the variance of the covariate coefficients and baseline cumulative distribution function of QIF method
#'
#' @param Time right censored data which is the follow up time.
#' @param Status the censoring indicator, 1 = event of interest happens, and 0 = censoring.
#' @param X a matrix of covariates that that may have effect on the failure times.
#' @param id a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.
#' @param beta initial value of covariates coefficients.
#' @param gS the parameters in the baseline cumulative distribution function.
#' @param Lambda initial value of baseline cumulative distribution function.
#' @param gSS3 estimation of baseline cumulative distribution function.
#' @param stat  the indepence indicator,1 =indepence, and 0 = non-zero variance of the fragile term
#' @param varrest the estimation of the variance of frailty factors.
#' @return Estimation of the variance of the covariate coefficients and baseline cumulative distribution function of QIF method.
#' @export

varrfra <- function(Time, Status, X, id, beta, gS, Lambda, gSS3,stat,varrest) {
  Kn<- length(id)
  K <- length(unique(id))
  n <- as.vector(table(id))
  xxxx<-X
  t2 <- Time
  c1 <- Status
  t11 <- sort(Time)
  c11 <- Status[order(Time)]
  tt1 <- unique(t11[c11 == 1])
  kk <- length(table(t11[c11 == 1]))
  betaest<-beta
  mu=exp(xxxx%*%betaest)
  gg1 <- rep(1,Kn)
  dy<-stats::deriv(expression(digamma(x)),"x",func=T)

  de<-function(x) as.numeric(attributes(dy(x))$gradient)

  if(stat){
    Pi=1-exp(-exp(xxxx%*%betaest))
    y1_hat=c1+(1-c1)*( exp(-exp(xxxx%*%betaest)*gSS3) +  Pi-1   )*(exp(-exp(xxxx%*%betaest)*gSS3)) ^(-1)


    BC1=matrix(0,dim(xxxx)[2],kk)
    for(v in 1:dim(xxxx)[2]){
    for(s in 1:(kk))
    {
      elem=0
      for(i in 1:K)
      {
        xxx1=xxxx[id==i,v]
        mu1=mu[id==i]
        t21=t2[id==i]
        c21=c1[id==i]
        gSS31=gSS3[id==i]
        for(j in 1:n)
        {
          elem=elem-c21[j]*sum(as.numeric(tt1<=t21[j])*gS)/(( sum(gS))^(2))*(mu1[j]*xxx1[j])-(1-c21[j])*sum(as.numeric(tt1<=t21[j])*gS)/(( sum(gS))^(2))*exp(-mu1[j]*(1-gSS31[j]))*(mu1[j])^2*xxx1[j]
          if(t21[j]>=tt1[s])
            elem=elem+c21[j]*1/(( sum(gS)))*(mu1[j]*xxx1[j])+(1-c21[j])*1/(( sum(gS)))*exp(-mu1[j]*(1-gSS31[j]))*(mu1[j])^2*xxx1[j]
        }
      }
      BC1[v,s]=elem
    }
    }

    BBC <- matrix(0, kk, dim( xxxx)[2])
    for (j in 1:dim( xxxx)[2]) {
      for (s in 1:(kk)) {
        BCm <- gS[s] * exp(xxxx[(c1 == 1) & (t2 == tt1[s]), ,drop = FALSE] %*% betaest)
        BBC[s, j] <- sum(exp(xxxx[(c1 == 1) & (t2 == tt1[s]), ,drop = FALSE] %*% betaest) * (exp(-BCm) + BCm * exp(-BCm) - 1)/((1 - exp(-BCm))^2) * xxxx[(c1 == 1) & (t2 == tt1[s]), j]) + sum(gg1[t2 >=
                                                                                                                                                                                                tt1[s]] * exp(xxxx[t2 >= tt1[s], ,drop = FALSE] %*% betaest) * xxxx[t2 >= tt1[s], j])
      }
    }
    CCC <- rep(0, (kk))
    for (s in 1:(kk)) {
      CCm <- gS[s] * exp(xxxx[(c1 == 1) & (t2 == tt1[s]), ,drop = FALSE] %*% betaest)
      CCC[s] <- sum(exp(2 * (xxxx[(c1 == 1) & (t2 == tt1[s]), ,drop = FALSE] %*% betaest) - CCm)/(1 - exp(-CCm))^2)
    }

    Ubb<-matrix(0,dim(xxxx)[2],dim(xxxx)[2])
    for(v in 1:dim(xxxx)[2]){
      for(s in 1:dim(xxxx)[2]){

        ela=exp(betaest*xxxx[,s])*xxxx[,s]*xxxx[,s]*(y1_hat-1-c1*gSS3)
        elb=(c1-1)*(exp(betaest*xxxx[,s])*gSS3*xxxx[,s]*xxxx[,s])
        elc=(1-c1)*(exp(-exp(betaest*xxxx[,s])*(1-gSS3)))*exp(betaest*xxxx[,s])*(exp(betaest*xxxx[,s])*xxxx[,s]*xxxx[,s]*(gSS3-1)+xxxx[,s]*xxxx[,s])

        Ubb[v,s]<- sum(ela)+sum(ela)+sum(ela)
      }
    }

    Q= rbind(   cbind(Ubb,t(BC1)),
                cbind(BBC,diag(CCC)))

    fdm=matrix(0,dim(xxxx)[2]+kk+1,dim(xxxx)[2]+kk+1)
    funest1=rep(0,dim(xxxx)[2])
    for(i in 1:K)
    {

      xxxi=xxxx[id==i,,drop = FALSE]
      c1i=c1[id==i]
      mu22 <- mu[id == i]
      y1_hati=y1_hat[id==i]
      gSS3i=gSS3[id==i]
      t22 <- t2[id == i]
      funest1=rep(0,dim(xxxx)[2])
      for(s in 1:dim(xxxx)[2]){
        fun1=exp(xxxi%*%betaest)*xxxi[,s]*(y1_hati-1-c1i*gSS3i)
        fun2=c1i*xxxi[,s]
        fun3=(1-c1i)*(exp(-exp(xxxi%*%betaest)*gSS3i)*(-exp(xxxi%*%betaest)*xxxi[,s]*gSS3i)+exp(-exp(xxxi%*%betaest))*exp(xxxi%*%betaest)*xxxi[,s] )/(exp(-exp(xxxi%*%betaest)*gSS3i))
        funest1[s]=sum(fun1)+sum(fun2)+sum(fun3)
      }
      eqalpha <- rep(0, kk)
      for (j in 1:kk) {
        Aalpha <- (1:n[i])[t22 == tt1[j]]
        Balpha <- (1:n[i])[t22 >= tt1[j]]
        if (length(Balpha) == 0) {
          eqalpha[j] <- 0
        }
        if (length(Balpha) != 0 & length(Aalpha) == 0) {
          eqalpha[j] <- -sum(mu22[Balpha])
        } else eqalpha[j] <- sum(mu22[Aalpha]/(1 - exp(-gS[j] * mu22[Aalpha]))) - sum(mu22[Balpha])
      }



      fdm=fdm+t(t(c(funest1,eqalpha)))%*%t(c(funest1,eqalpha))

    }

    vcm<-  ( MASS::ginv(Q)%*%fdm%*%t(MASS::ginv(Q))  )
  }
  else{
    Pi=1-exp(-exp(xxxx%*%betaest))
    y1_hat=c1+(1-c1)*( exp(-exp(xxxx%*%betaest)*Lambda) +  Pi-1   )*(exp(-exp(xxxx%*%betaest)*Lambda)) ^(-1)
    B<-rep(0,K)
    D<-rep(0,K)
    for (i in 1:K) {
      B[i]=1/varrest+sum(c1[id==i])

    }


    for (i in 1:K) {
      D[i]=1/varrest+sum(((Lambda*c1+1-y1_hat)*exp(xxxx%*%betaest))[id==i]  )
    }
    W<-B/D


    Ubb<-matrix(0,dim(xxxx)[2],dim(xxxx)[2])
    for(v in 1:dim(xxxx)[2]){
    for(s in 1:dim(xxxx)[2]){
      fu1=W[id]*exp(xxxx%*%betaest)*xxxx[,v]*xxxx[,s]*(y1_hat-1-c1*Lambda)
      fu2=(c1-1)*(W[id]*exp(xxxx%*%betaest)*Lambda*xxxx[,v]*xxxx[,s])
      fu3=(1-c1)*(exp(-W[id]*exp(xxxx%*%betaest)*(1-Lambda)))*W[id]*exp(xxxx%*%betaest)*(W[id]*exp(xxxx%*%betaest)*xxxx[,v]*xxxx[,s]*(Lambda-1)+xxxx[,v]*xxxx[,s])
      Ubb[v,s]=sum(fu1)+sum(fu2)+sum(fu3)
    }
    }

    Uyy= varrest^(-3)*(-2*(K*digamma(1/varrest)+K*log(varrest)-K+sum(-(digamma(B)-log(D))+W))+K*de(1/varrest)*(-1/varrest)+K)

    Q= rbind(   cbind(Ubb,as.matrix(rep(0,dim(xxxx)[2]))) ,
                cbind(t(as.matrix(rep(0,dim(xxxx)[2]))),Uyy) )
    fdm=matrix(0,dim(xxxx)[2]+1,dim(xxxx)[2]+1)

    for(i in 1:K)
    {

      xxxi=xxxx[id==i,,drop = FALSE]
      c1i=c1[id==i]
      Wi=W[i]
      Bi=B[i]
      Di=D[i]
      y1_hati=y1_hat[id==i]
      gSS3i=Lambda[id==i]
      t2i=t2[id==i]
      funest1=rep(0,dim(xxxx)[2])
      for(s in 1:dim(xxxx)[2]){
        fun1=Wi*exp(xxxi%*%betaest)*xxxi[,s]*(y1_hati-1-c1i*gSS3i)
        fun2=c1i*xxxi[,s]
        fun3=(1-c1i)*(exp(-Wi*exp(xxxi%*%betaest)*gSS3i)*(-Wi*exp(xxxi%*%betaest)*xxxi[,s]*gSS3i)+exp(-Wi*exp(xxxi%*%betaest))*Wi*exp(xxxi%*%betaest)*xxxi[,s] )/(exp(-Wi*exp(xxxi%*%betaest)*gSS3i))
        funest1[s]=sum(fun1)+sum(fun2)+sum(fun3)
      }
      Uy= varrest^(-2)*(K*digamma(varrest^(-1))+K*log(varrest)-K+sum(-(digamma(Bi)-log(Di))+Wi))
      fdm=fdm+t(t(c(funest1,Uy)))%*%t(c(funest1,Uy))
    }

    vcm<-  ( MASS::ginv(Q)%*%fdm%*%t(MASS::ginv(Q))  )
  }

  var_beta <- diag(vcm)[1: dim(xxxx)[2]]
  sd_beta <- sqrt(var_beta)

  list(var_beta = var_beta, sd_beta = sd_beta)
}


