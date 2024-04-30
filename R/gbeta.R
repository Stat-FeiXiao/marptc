#' Title Estimate the covariate coefficients using the GEE approach and baseline cumulative distribution function in model using Bernstein polynomial
#'
#' @description Iteration between the estimation of covariate coefficients and baseline cumulative distribution funtion. The iteration will be stoped once the relative change
#' in deviance is less than \code{eps} or the number of iterations reach \code{itermax}.
#'
#' @param Time right censored data which is the follow up time.
#' @param Status the censoring indicator, normally 0 = event of interest happens, and 0 = censoring.
#' @param X a matrix of covariates that may have effect on the failure times.
#' @param stad if it is TRUE, all the covariates in the \code{formula} are standardized. By default, \code{stdz = FALSE}.
#' @param Ibeta the initial value of the covariate coefficients.
#' @param id a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.
#' @param itermax specifies the maximum iteration number. If the convergence criterion is not met, the iteration will be stopped after \code{itermax} iterations and
#' the estimates will be based on the last iteration. The default \code{itermax = 100}.
#' @param eps tolerance for convergence. The default is \code{eps = 1e-6}. Iteration stops once the relative change in deviance is less than \code{eps}.
#' @param N the degree of Bernstein polynomial used to approximate the baseline distribution function in model.
#' @param tau the cure threshold. Individuals with survival time greater than \code{tau} are considered cured.
#' @param corstr a character string specifying the correlation structure. The following are permitted: \code{independence}, \code{exchangeable} and \code{AR1}.
#'
#' @export gbeta
gbeta <- function(Time, Status, X,stad, Ibeta, id, itermax, eps,N,tau,corstr) {
  Kn <- length(id)
  K <- length(unique(id))
  n <- as.vector(table(id))
  stadX <- X
  if(stad){
   for (i in 2:ncol(X1)) {
                stadX[, i] <- (X[, i] - mean(X[, i]))/sd(X[, i])
            }
  }
  
  beta2 <- beta <- Ibeta
  Fhat <- initialF(Time, Status, stadX, id)$Fhat

  iter <- 0
  repeat{
    iter1<-0
    repeat {
      gbeta <- geq(Status, Fhat, X=stadX, beta = beta2,  id,eps,corstr)
      beta1 <- gbeta$beta
      if ((any(abs(beta1 - beta2) > eps)) & (iter1<= itermax) ) {
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
 rho <- gbeta$rho
 pphi <- gbeta$pphi

   if(stad){
 beta2<- c(beta1[1] - sum((beta1[-1] * apply(X[, -1, drop = FALSE], 2, mean)/apply(X[, -1, drop = FALSE], 2, sd))), beta1[-1]/apply(X[,  -1, drop = FALSE], 2, sd))
  }else{ beta2 <- beta1 }

  convergence <- sqrt(sum(((beta1 - beta)^2)))
  list(beta1 = beta1, beta2=beta2,Fhat = Fhat, rho = rho,
             pphi = pphi,  gam=gam,convergence = convergence)
}


#' Title General estimation equations for the covariate coefficients
#'
#' @description Estimation the covariate coefficients using the GEE approach.
#'
#' @param Status right censored data which is the follow up time.
#' @param baseF estimation of baseline cumulative distribution function obtained by Bernstein polynomial.
#' @param X a matrix of covariates that that may have effect on the failure times.
#' @param beta the initial value of the covariate coefficients.
#' @param id a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.
#' @param eps tolerance for convergence. The default is \code{eps = 1e-6}. Iteration stops once the relative change in deviance is less than \code{eps}.
#' @param corstr a character string specifying the correlation structure. The following are permitted: \code{independence}, \code{exchangeable} and \code{AR1}.
#'
#' @export geq
geq <- function(Status, baseF, X, beta, id,eps,corstr) {
  K <- length(unique(id))
  n <- as.vector(table(id))
  newY1 <- Status/baseF
  W1 <- diag( baseF )

  gbeta <- beta1 <-beta

  mu <- exp(X %*% gbeta )

  res <- as.vector((newY1 - mu)/sqrt(mu))
  rres <- 0
  pphi <- sum(res^2)/(sum(n) - dim(X)[2] )
  resm <- matrix(0, ncol = K, nrow = max(n))
  for (i in 1:K) {
    resm[1:n[i], i] <- res[id == i]
  }
  res <- resm
  res <- t(res)
  if(corstr=="independence"){
    rho <- 0
  }else if(corstr=="exchangeable"){
    for (i in 1:K) {
      if (n[i] == 1){
        rres <- rres + res[i, 1] }else {
          for (j in 1:(n[i] - 1)) rres <- rres + res[i, j] * sum(res[i, (j + 1):n[i]])
        }
    }
    rho <- (pphi^(-1)) * rres/(sum(n * (n - 1))/2 -dim(X)[2]  )

  }else if(corstr=="AR1"){
    for (i in 1:K) {
      if (n[i] == 1)
        rres <- rres + res[i, 1] else {
          for (j in 1:(n[i] - 1)) rres <- rres + res[i, j] * res[i, j+1]
        }
    }
    rho <- (pphi^(-1)) * rres/(sum((n - 1)) - dim(X)[2]  )

  }

  repeat{
    mu <- exp(X %*% gbeta )
    S1 <- newY1 - mu
    if(corstr=="independence"){
      QC <- diag(n[1])
    }else if(corstr=="exchangeable"){
      QC=matrix(rho,n[1],n[1])
      diag(QC)=1
    }else if(corstr=="AR1"){
      exponent <- abs(matrix(1:n[1] - 1, nrow = n[1], ncol = n[1], byrow = TRUE) - (1:n[1] - 1))
      QC=rho^exponent
    }
    W1=diag(baseF[id==1])
    D1=diag(mu[id==1])%*%X[id==1,]
    V1=ginv( sqrt(diag(mu[id==1])) )%*%(pphi^(-1)*ginv(QC))%*%ginv( sqrt(diag(mu[id==1])) )
    G=t(D1)%*%V1%*%W1%*%S1[id==1]
    G1=-t(D1)%*%V1%*%W1%*%D1


    for(i in 2:K)
    {
      if(corstr=="independence"){
        QC <- diag(n[i])
      }else if(corstr=="exchangeable"){
        QC=matrix(rho,n[i],n[i])
        diag(QC)=1
      }else if(corstr=="AR1"){
        exponent <- abs(matrix(1:n[i] - 1, nrow = n[i], ncol = n[i], byrow = TRUE) - (1:n[i] - 1))
        QC=rho^exponent
      }

      W1=diag(baseF[id==i])
      D1=diag(mu[id==i])%*%diag(rep(1,n[i]))%*%(X[id==i,])
      V1=ginv( sqrt(diag(mu[id==i])) )%*%(pphi^(-1)*ginv(QC))%*%ginv( sqrt(diag(mu[id==i])) )

      G_1=t(D1)%*%V1%*%W1%*%S1[id==i]
      G=G+G_1


      G11=-t(D1)%*%V1%*%W1%*%D1
      G1=G1+G11

    }


    a <- G1
    storage.mode(a)<-"double"
    lda<-as.integer(nrow(a))
    npp<-as.integer(ncol(a))
    z <- .Fortran("finv", a=a, lda, npp)
    G1inv <- z$a

    gbeta=gbeta-G1inv%*%G

    if( any(abs( gbeta-beta1)>eps) ) {
      beta1 <- gbeta
    }else break
  }

  list(beta = gbeta, rho = rho, pphi = pphi)
}
