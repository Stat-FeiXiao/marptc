#' Title Estimate the covariate coefficients using the QIF approach and baseline cumulative distribution function in model using Bernstein polynomial
#'
#' @description Iteration between the estimation of covariate coefficients and baseline cumulative distribution funtion. The iteration will be stoped once the relative change
#' in deviance is less than \code{eps} or the number of iterations reach \code{itermax}.
#'
#' @param Time right censored data which is the follow up time.
#' @param Status the censoring indicator, normally 0 = event of interest happens, and 0 = censoring
#' @param X a matrix of covariates that may have effect on the failure times.
#' @param stad if it is TRUE, all the covariates in the \code{formula} are standardized. By default, \code{stdz = FALSE}.
#' @param Ibeta the initial value of the covariate coefficients.
#' @param id a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.
#' @param itermax specifies the maximum iteration number. If the convergence criterion is not met, the iteration will be stopped after \code{esmax} iterations and
#' the estimates will be based on the last iteration. The default \code{itermax = 100}.
#' @param eps tolerance for convergence. The default is \code{eps = 1e-6}. Iteration stops once the relative change in deviance is less than \code{eps}.
#' @param N the degree of Bernstein polynomial used to approximate the baseline distribution function in model.
#' @param tau the cure threshold. Individuals with survival time greater than \code{tau} are considered cured.
#' @param corstr a character string specifying the correlation structure. The following are permitted: \code{independence}, \code{exchangeable} and \code{AR1}.
#'
#' @export qbeta

qbeta <- function(Time, Status, X,stad,Ibeta,  id, itermax, eps,N,tau,corstr) {
  Kn <- length(id)
  K <- length(unique(id))
  n <- as.vector(table(id))
 stadX <- X
  if(stad){
   for (i in 2:ncol(X1)) {
                stadX[, i] <- (X[, i] - mean(X[, i]))/sd(X[, i])
            }
  }else{
    stadX <- X
  }
 beta2 <- beta <- Ibeta
  Fhat <- initialF(Time, Status, stadX, id)$Fhat

  iter <- 0
  repeat{
    iter1<-0
    repeat {
      qbeta <- qeq(Status, Fhat, X=stadX, beta = beta2,  id,eps,corstr)
      beta1 <- qbeta$beta
      if ((any(abs(beta1 - beta2) > eps)) & (iter1 <= itermax) ) {
       beta2 <- beta1
       iter1<-iter1+1
      } else  break
    }
    basecdf<- baseF(Time, Status, X = stadX, beta = beta1,N,tau)
    Fhat <- basecdf$Fhat
    if ((any(abs(beta1 - beta) > eps)) & (iter<= itermax)) {
      beta <- beta1
      iter<-iter+1
    } else  break
  }

gam <- basecdf$gam

   if(stad){
 beta2<- c(beta1[1] - sum((beta1[-1] * apply(X1[, -1, drop = FALSE], 2, mean)/apply(X1[, -1, drop = FALSE], 2, sd))), beta1[-1]/apply(X1[,  -1, drop = FALSE], 2, sd))
  }else{ beta2 <- beta1 }

  convergence <- sqrt(sum(((beta1 - beta)^2)))
  list(beta1 = beta1,beta2 = beta2, Fhat = Fhat, gam=gam, convergence = convergence)
}

#' Title Estimation the covariate coefficients using the QIF approach
#'
#' @description Estimation the covariate coefficients using the QIF approach.
#'
#' @param Status right censored data which is the follow up time.
#' @param baseF estimation of baseline cumulative distribution function obtained by Bernstein polynomial.
#' @param X a matrix of covariates that that may have effect on the failure times.
#' @param beta the initial value of the covariate coefficients
#' @param id a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.
#' @param eps tolerance for convergence. The default is \code{eps = 1e-6}. Iteration stops once the relative change in deviance is less than \code{eps}.
#' @param corstr a character string specifying the correlation structure. The following are permitted: \code{independence}, \code{exchangeable} and \code{AR1}.
#'
#' @export qeq
qeq <- function(Status, baseF, X, beta, id,eps,corstr) {
  K <- length(unique(id))
  n <- as.vector(table(id))
  newY1 <- Status/baseF
  W1 <- diag(baseF)
  qbeta <- beta1 <- beta

  repeat{
    mu <- exp(X %*% qbeta )
    S1 <- newY1 - mu
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
      }else{
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
      G_1=t(D1)%*%V1%*%W1%*%S1[id==i]
      VA1 <- -t(D1)%*%V1%*%W1%*%D1
      if( corstr!="independence"){
        V2=MASS::ginv( sqrt(diag(mu[id==i])) )%*%M2%*%MASS::ginv( sqrt(diag(mu[id==i])) )
        G_1=rbind(G_1,t(D1)%*%V2%*%W1%*%S1[id==i])
        VA1=rbind( VA1,-t(D1)%*%V2%*%W1%*%D1)
      }
      G=G+G_1
      C=C+G_1%*%t(G_1)
      G1<-G1 +  VA1
    }

    a <- C
    storage.mode(a)<-"double"
    lda<-as.integer(nrow(a))
    npp<-as.integer(ncol(a))
    z <- .Fortran("finv", a=a, lda, npp)
    Cinv <- z$a

    Q1=t(G1)%*%Cinv%*%G
    Q2=t(G1)%*%Cinv%*%G1

    a <- Q2
    storage.mode(a)<-"double"
    lda<-as.integer(nrow(a))
    npp<-as.integer(ncol(a))
    z <- .Fortran("finv", a=a, lda, npp)
    Q2inv <- z$a

    qbeta=qbeta-Q2inv%*%Q1

    if( any(abs( qbeta-beta1)>eps)) {
      beta1 <- qbeta
    }else break
  }

  list(beta = qbeta)
}

