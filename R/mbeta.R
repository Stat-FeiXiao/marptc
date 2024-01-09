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
oldbeta1 <- oldbeta  <- rep(1000,length(as.double(mbeta)))

  mu <- exp(xxxx %*% mbeta )
  SK1 <- 1
  res <- as.vector((newY1 - mu)/sqrt(mu))
  rres <- 0
  pphi <- sum(res^2)/(sum(n) - dim(X)[2] ) 
  resm <- matrix(0, ncol = K, nrow = max(n))
  for (i in 1:K) {
    resm[1:n[i], i] <- res[id == i]
  }
  res <- resm
  res <- t(res)

  if(struc=="exchangeable"){
  for (i in 1:K) {
    if (n[i] == 1){
      rres <- rres + res[i, 1] }else {
        for (j in 1:(n[i] - 1)) rres <- rres + res[i, j] * sum(res[i, (j + 1):n[i]])
      }
  }
  rho <- (pphi^(-1)) * rres/(sum(n * (n - 1))/2 -dim(X)[2] )

  }else if(struc=="AR1"){
    for (i in 1:K) {
      if (n[i] == 1)
        rres <- rres + res[i, 1] else {
          for (j in 1:(n[i] - 1)) rres <- rres + res[i, j] * res[i, j+1]
        }
    }
    rho <- (pphi^(-1)) * rres/(sum((n - 1)) - dim(X)[2]  )
    
  }

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
  V1=ginv( sqrt(diag(mu[id==1])) )%*%(pphi^(-1)*ginv(QC))%*%ginv( sqrt(diag(mu[id==1])) )    #
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
    V1=ginv( sqrt(diag(mu[id==i])) )%*%(pphi^(-1)*ginv(QC))%*%ginv( sqrt(diag(mu[id==i])) )    #

    G_1=t(D1)%*%V1%*%W1%*%S1[id==i]
    G=G+G_1

   
    G11=-t(D1)%*%V1%*%W1%*%D1
    G1=G1+G11

  }


  ABC1 <- rep(0, K)
  VA1 <- matrix(0, dim(X)[2], dim(X)[2])
        for (v in 1:dim(X)[2]) {
            for (w1 in 1:dim(X)[2]) {
                for (i in 1:K) {
                
  if(struc=="exchangeable"){
      Q1=matrix(rho,n[i],n[i])
      diag(Q1)=1
    }else if(struc=="AR1"){
      exponent <- abs(matrix(1:n[i] - 1, nrow = n[i], ncol = n[i], byrow = TRUE) - (1:n[i] - 1))
      Q1=rho^exponent
    }
                  IQ1 <- ginv(Q1)
                  B2 <- matrix(0, n[i], n[i])
                  z22 <- matrix(X[id == i, ], nrow = n[i], )
                  A2 <- t(z22[, v])
                  c22 <- Status[id == i]
                  mu22 <- mu[id == i]
                  Lambda22 <- Lambda[id == i]
                  BB1 <- (mu22^(1/2)) %*% ((t(mu22))^(-1/2)) * IQ1
                  for (s in 1:n[i]) {
                    for (l in 1:n[i]) {
                      B2[s, l] <- (1/2) * (z22[s, w1] - z22[l, w1]) * BB1[s, l]
                    }
                  }
                  C2 <- c22/Lambda22 - mu22
                  D2 <- BB1
                  E2 <- z22[, w1] * mu22
                  G2 <- diag(Lambda22)
                  ABC1[i] <- A2 %*% (B2 %*% G2 %*% C2 - D2 %*% G2 %*% E2)
                }
                VA1[v, w1] <- sum(ABC1) * (pphi^(-1))
                ABC1 <- rep(0, K)
            }
        }

  Q1=G
  Q2=G1#  VA1#

  mbeta=mbeta-ginv(Q2)%*%Q1

  if( all(abs( mbeta-oldbeta)<1e-6)) {
oldbeta <- mbeta
}else break
  }

  list(beta = mbeta, rho = rho, pphi = pphi)
}
