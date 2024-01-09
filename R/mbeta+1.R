#' Title Estimation of the covariate coefficients using GEE method
#'
#' @param Status the censoring indicator, 1 = event of interest happens, and 0 = censoring.
#' @param Lambda initial cumulative baseline hazard function.
#' @param X a matrix of covariates that that may have effect on the failure times.
#' @param beta initial values of covariates coefficients.
#' @param id a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.
#'
#' @return estimation of the covariate coefficients and parameters in the working matrix.
#' @export

mbeta <- function(Status, Lambda, X, beta, id,struc) {
  K <- length(unique(id))
  n <- as.vector(table(id))
  newY1 <- Status/Lambda
  W1 <- diag( Lambda)
  xxxx<-X
  mbeta<-beta
  mu <- exp(xxxx %*% mbeta )
  SK1 <- 1
  res <- as.vector((newY1 - mu)/sqrt(mu))
  rres <- 0
  pphi <- sum(res^2)/(sum(n) - dim(X)[2])
  resm <- matrix(0, ncol = K, nrow = max(n))
  for (i in 1:K) {
    resm[1:n[i], i] <- res[id == i]
  }
  res <- resm
  res <- t(res)
  for (i in 1:K) {
    if (n[i] == 1)
      rres <- rres + res[i, 1] else {
        for (j in 1:(n[i] - 1)) rres <- rres + res[i, j] * sum(res[i, (j + 1):n[i]])
      }
  }
  rho <- (pphi^(-1)) * rres/(sum(n * (n - 1))/2 - (dim(X)[2]+1))


  repeat{
    mu <- exp(xxxx %*% mbeta )
    S1 <- newY1 - mu
    if(struc=="exchangeable"){
      QC=matrix(rho,n[1],n[1])
      diag(QC)=1
    }else if(struc=="AR1"){
      exponent <- abs(matrix(1:n[1] - 1, nrow = n[1], ncol = n[1], byrow = TRUE) - (1:n[1] - 1))
      QC=rho^exponent
    }
  W1=diag(Lambda[id==1])
  D1=diag(mu[id==1])%*%diag(rep(1, n[1])) %*%(xxxx[id==1,])
  V1=ginv( sqrt(diag(mu[id==1])) )%*%(pphi^(-1)*ginv(QC))%*%ginv( sqrt(diag(mu[id==1])) )
  G=t(D1)%*%V1%*%W1%*%S1[id==1]


 
  G1=-t(D1)%*%V1%*%W1%*%D1


  for(i in 2:K)
  {
    if(struc=="exchangeable"){
      QC=matrix(rho,n[i],n[i])
      diag(QC)=1
    }else if(struc=="AR1"){
      exponent <- abs(matrix(1:n[i] - 1, nrow = n[i], ncol = n[i], byrow = TRUE) - (1:n[i] - 1))
      QC=rho^exponent
    }

    W1=diag(Lambda[id==i])
    D1=diag(mu[id==i])%*%diag(rep(1,n[i]))%*%(xxxx[id==i,])
    V1=solve( sqrt(diag(mu[id==i])) )%*%(pphi^(-1)*ginv(QC))%*%solve( sqrt(diag(mu[id==i])) )

    G_1=t(D1)%*%V1%*%W1%*%S1[id==i]
    G=G+G_1

   
    G11=-t(D1)%*%V1%*%W1%*%D1
    G1=G1+G11

  }


  VA1 <- matrix(0, dim(xxxx)[2], dim(xxxx)[2])
  newG11<-matrix(0, dim(xxxx)[2], 1)
  for (v in 1:dim(xxxx)[2]) {

    for (i in 1:K) {
      if(struc=="exchangeable"){
        QC=matrix(rho,n[i],n[i])
        diag(QC)=1
      }else if(struc=="AR1"){
        exponent <- abs(matrix(1:n[i] - 1, nrow = n[i], ncol = n[i], byrow = TRUE) - (1:n[i] - 1))
        QC=rho^exponent
      }

      W1=diag(Lambda[id==i])
      D1=diag(mu[id==i])%*%diag(rep(1,n[i]))%*%(xxxx[id==i,])
      V1=MASS::ginv( sqrt(diag(mu[id==i])) )%*%(pphi^(-1)*MASS::ginv(QC))%*%MASS::ginv( sqrt(diag(mu[id==i])) )


      
      newD1=diag(mu[id==i])%*%diag(rep(1,n[i]))%*%(xxxx[id==i,v])
      G11=-t(D1)%*%V1%*%W1%*%newD1

      newG11<-newG11+G11
    }
    VA1[, v] <- newG11
    newG11<-matrix(0, dim(xxxx)[2], 1)

  }
  Q1=G
  Q2=VA1

  oldbeta=mbeta
  mbeta=mbeta-ginv(Q2)%*%Q1

  if( all(abs( mbeta-oldbeta)<1e-6)) break
  }


  list(beta = mbeta, rho = rho, pphi = pphi)
}
