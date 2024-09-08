#' Title Estimation the covariate coefficients using the QIF approach
#'
#' @description Estimation the covariate coefficients using the QIF approach.
#'
#' @param Time right censored data which is the follow up time.
#' @param Status right censored data which is the follow up time.
#' @param X a matrix of covariates that that may have effect on the failure times.
#' @param beta the initial value of the covariate coefficients.
#' @param id a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.
#' @param eps tolerance for convergence. The default is \code{eps = 1e-6}. Iteration stops once the relative change in deviance is less than \code{eps}.
#' @param corstr a character string specifying the correlation structure. The following are permitted: \code{independence}, \code{exchangeable} and \code{AR1}.
#' @param itermax specifies the maximum iteration number. If the convergence criterion is not met, the iteration will be stopped after \code{itermax} iterations and
#' the estimates will be based on the last iteration.
#'
#' @export qeq

qeq <- function(Time, Status, X, beta, id, eps, corstr, itermax) {


  t11 <- sort(Time[Status==1])
  kk <- sum(Status)
  K <- length(unique(id))
  n <- as.vector(table(id))
  qbeta <- beta1 <-beta

  mu <- exp(X %*% qbeta )

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
  SK1 <- 1
  repeat{
    mu <- exp(X %*% qbeta )
    S1 <- Status - mu * baseF
    if(corstr=="independence"){
      G=matrix(0,dim(X)[2],1)
      C=matrix(0,dim(X)[2],dim(X)[2])
      G1=matrix(0,dim(X)[2],dim(X)[2])
    }else{
      G=matrix(0,2*dim(X)[2],1)
      C=matrix(0,2*dim(X)[2],2*dim(X)[2])
      G1=matrix(0,2*dim(X)[2],dim(X)[2])
    }
    for (i in 1:K) {
      M1 <-  diag(n[i])
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

      W1=diag(baseF[id==i])
      D1=diag(mu[id==i])%*%(as.matrix(X[id==i,]))

      V1=MASS::ginv( sqrt(diag(mu[id==i])) )%*%M1%*%MASS::ginv( sqrt(diag(mu[id==i])) )
      G_1=t(D1)%*%V1%*%S1[id==i]
      VA1 <- -t(D1)%*%V1%*%W1%*%D1
      if( corstr!="independence"){
        V2=MASS::ginv( sqrt(diag(mu[id==i])) )%*%M2%*%MASS::ginv( sqrt(diag(mu[id==i])) )
        G_1=rbind(G_1,t(D1)%*%V2%*%S1[id==i])
        VA1=rbind( VA1,-t(D1)%*%V2%*%W1%*%D1)
      }
      G=G+G_1
      C=C+G_1%*%t(G_1)
      G1<-G1 +  VA1
    }


    Q1=t(G1)%*%ginv(C)%*%G
    Q2=t(G1)%*%ginv(C)%*%G1


    qbeta=qbeta-ginv(Q2)%*%Q1

    mu <- exp(X %*% qbeta )

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

    S1 <- (Status - mu*baseF)

    if( any(abs( qbeta-beta1)>eps) & SK1 <= itermax) {
      beta1 <- qbeta
    }else break
  }

  convergence <- (SK1 <= itermax)
  list(beta = qbeta, convergence=convergence)
}

