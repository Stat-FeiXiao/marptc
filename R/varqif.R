#' Title  Estimate Variance of coefficients estimation from the QIF approach with a sandwich formula
#'
#' @description Calculate the variance estimates using the sandwich formula.
#'
#' @param Time right censored data which is the follow up time.
#' @param Status the censoring indicator, normally 0 = event of interest happens, and 0 = censoring
#' @param X a matrix of covariates that may have effect on the failure times.
#' @param id a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.
#' @param beta covariate coefficients estimation from QIF method.
#' @param corstr a character string specifying the correlation structure. The following are permitted: \code{independence}, \code{exchangeable} and \code{AR1}.
#'
#' @export varqif

varqif <- function(Time, Status, X, id, beta,corstr) {
  Kn<- length(id)
  K <- length(unique(id))
  n <- as.vector(table(id))
  t2 <- Time
  tt1 <- t2[Status==1]
  t11 <- sort(tt1)
  c1 <- Status
  kk <- length(table(t2[c1 == 1]))
  mu <- exp(X %*% beta )

  dRbet <- sapply(1:kk,function(i){
    sum( mu*(Time>=t11[i]) )
  })
  Rbet.min <- min(dRbet)
  interval <- c(Rbet.min-sum(Status),Rbet.min-1)
  lambet <- uniroot(getlambda,interval,tol = .Machine$double.eps^0.75,
                    dRbet=dRbet)$root


  baseF <- sapply(Time,function(tmj){
    sum( (t11<=tmj)/(dRbet-lambet) )
  })
  baseF1 <- baseF
  baseF1[baseF1==0] <- 1e-8
  S1 <- (Status - mu*baseF)

  be <- beta
  mu=exp(X%*%beta)
  g1 <- rep(1,Kn)
  xxxx <- X
  gg1 <- g1


  if(corstr=="independence"){
    VA1 <- matrix(0, dim(X)[2], dim(X)[2])
    G=matrix(0,dim(X)[2],1)
    C=matrix(0,dim(X)[2],dim(X)[2])
  }else{
    VA1 <- matrix(0, 2*dim(X)[2], dim(X)[2])
    G=matrix(0,2*dim(X)[2],1)
    C=matrix(0,2*dim(X)[2],2*dim(X)[2])
  }


  for (i in 1:K) {
    M1 <-  diag(n[i])
    if(corstr=="exchangeable"){
      M2 <-	matrix(1,n[i],n[i])
      diag(M2) <- 0
    }else if(corstr=="AR1"){
      M2 <-	matrix(0,n[i],n[i])
      M2[1,2]<-1
      for(o in 2:n[i]-1){
        M2[o,o+1] <-1
        M2[o,o-1] <-1
      }
      M2[n[i],n[i]-1]<-1
    }

    W1=diag(baseF[id==i])
    D1=diag(mu[id==i])%*%X[id==i,]
    V1=MASS::ginv( sqrt(diag(mu[id==i])) )%*%M1%*%MASS::ginv( sqrt(diag(mu[id==i])) )
    G_1=t(D1)%*%V1%*%S1[id==i]
    G11=-t(D1)%*%V1%*%W1%*%D1

    if(corstr!="independence"){
      V2=MASS::ginv( sqrt(diag(mu[id==i])) )%*%M2%*%MASS::ginv( sqrt(diag(mu[id==i])) )
      G_1=rbind(G_1,t(D1)%*%V2%*%S1[id==i])
      G11=rbind(G11,-t(D1)%*%V2%*%W1%*%D1)
    }

    G=G+G_1
    C=C+G_1%*%t(G_1)
    VA1<-VA1+G11
  }

  M= -t(VA1)%*% ginv(C)%*%VA1

  fdv <- rep(0, dim(X)[2])
  fdm <- matrix(0,dim(X)[2], dim(X)[2] )
  for (i in 1:K) {
    z22 <- X[id == i, ]
    mu22 <- mu[id == i]
    mu22m <- diag(mu22, n[i], n[i])
    baseF01 <- baseF[id == i]
    G2 <- diag( baseF01)
    c22 <- c1[id == i]
    g11 <- g1[id == i]
    C2 <- (c22) - mu22 *baseF01
    M1 <- diag(n[i])
    if(corstr=="exchangeable"){
      M2 <- matrix(1,n[i],n[i])
      diag(M2) <- 0
    }else if(corstr=="AR1"){
      M2 <- matrix(0,n[i],n[i])
      M2[1,2]<-1
      for(o in 2:n[i]-1){
        M2[o,o+1] <-1
        M2[o,o-1] <-1
      }
      M2[n[i],n[i]-1]<-1
    }
    if(corstr=="independence"){
      fdv0= t(mu22m%*%z22)%*%MASS::ginv(sqrt(mu22m))%*%M1%*%MASS::ginv(sqrt(mu22m))%*%C2
    }else{
      fdv0= rbind(t(mu22m%*%z22)%*%MASS::ginv(sqrt(mu22m))%*%M1%*%MASS::ginv(sqrt(mu22m))%*%C2 ,
                  t(mu22m%*%z22)%*%MASS::ginv(sqrt(mu22m))%*%M2%*%MASS::ginv(sqrt(mu22m))%*%C2 )
    }

    fdv <- t(VA1)%*%ginv(C)%*%fdv0
    fdm1 <- fdv %*% t(fdv)
    fdm <- fdm + fdm1
  }
  vcm <- ginv(M) %*% fdm %*% t(ginv(M))
  var_beta <- diag( as.matrix(vcm)  )
  sd_beta <- sqrt(var_beta)
  list(var_beta = var_beta, sd_beta = sd_beta)
}
