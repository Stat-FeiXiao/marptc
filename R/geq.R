#' Title General estimation equations for the covariate coefficients estimation
#'
#' @description Estimation the covariate coefficients using the GEE approach.
#'
#' @param Time right censored data which is the follow up time.
#' @param Status right censored data which is the follow up time.
#' @param X a matrix of covariates that that may have effect on the failure times.
#' @param beta the initial value of the covariate coefficients.
#' @param id a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.
#' @param eps tolerance for convergence. Iteration stops once the relative change in deviance is less than \code{eps}.
#' @param corstr a character string specifying the correlation structure. The following are permitted: \code{independence}, \code{exchangeable} and \code{AR1}.
#' @param itermax specifies the maximum iteration number. If the convergence criterion is not met, the iteration will be stopped after \code{itermax} iterations and
#' the estimates will be based on the last iteration.
#'
#' @export geq

geq <- function(Time, Status, X, beta, id, eps, corstr, itermax) {


  t11 <- sort(Time[Status==1])
  kk <- sum(Status)
  K <- length(unique(id))
  n <- as.vector(table(id))
  gbeta <- beta1 <-beta

  mu <- exp(X %*% gbeta )

  dRbet <- sapply(1:sum(n),function(i){
    sum( mu*(Time>=Time[i]) )
  })
  Rbet.min <- min(dRbet[Status==1])
  interval <- c(Rbet.min-sum(Status),Rbet.min-1)
  lambet <- uniroot(getlambda,interval,tol = .Machine$double.eps^0.75,
                    dRbet=dRbet,delta=Status)$root

  baseF <- sapply(Time,function(tmj){
    sum( (Time<=tmj)[Status==1]/(dRbet[Status==1]-lambet) )
  })
  baseF1 <- baseF
  baseF1[baseF1==0] <- 1e-08

  newY1 <- Status/baseF1
  W1 <- diag( baseF )



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

  S1 <- (Status - mu*baseF)
  SK1 <- 1
  repeat{

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
    G=t(D1)%*%V1%*%S1[id==1]
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
      D1=diag(mu[id==i])%*%(X[id==i,])
      V1=ginv( sqrt(diag(mu[id==i])) )%*%(pphi^(-1)*ginv(QC))%*%ginv( sqrt(diag(mu[id==i])) )

      G_1=t(D1)%*%V1%*%S1[id==i]
      G=G+G_1


      G11=-t(D1)%*%V1%*%W1%*%D1
      G1=G1+G11

    }



    gbeta=gbeta-ginv(G1)%*%G

    mu <- exp(X %*% gbeta )

    dRbet <- sapply(1:sum(n),function(i){
      sum( mu*(Time>=Time[i]) )
    })
    Rbet.min <- min(dRbet[Status==1])
    interval <- c(Rbet.min-sum(Status),Rbet.min-1)
    lambet <- uniroot(getlambda,interval,tol = .Machine$double.eps^0.75,
                      dRbet=dRbet,delta=Status)$root

    baseF <- sapply(Time,function(tmj){
      sum( (Time<=tmj)[Status==1]/(dRbet[Status==1]-lambet) )
    })

    S1 <- (Status - mu*baseF)

    if( any(abs( gbeta-beta1)>eps) & SK1 <= itermax ) {
      beta1 <- gbeta
      SK1 <- SK1 + 1
    }else break
  }
  convergence <- (SK1 <= itermax)
  list(beta = gbeta, rho = rho, convergence=convergence)
}
