#' Title Estimation of the covariate coefficients using GEE method improved by quadratic inference functions
#'
#' @param Status the censoring indicator, 1 = event of interest happens, and 0 = censoring.
#' @param Lambda initial cumulative baseline hazard function.
#' @param X a matrix of covariates that that may have effect on the failure times.
#' @param beta initial values of covariates coefficients.
#' @param id a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.
#'
#' @return estimation of the covariate coefficients and parameters in the working matrix.
#' @export

qbeta <- function(Status, Lambda, X, beta, id,eps) {
  K <- length(unique(id))
  n <- as.vector(table(id))
  newY1 <- Status/Lambda
  W1 <- diag(Lambda)
  xxxx<-X
  qbeta<-beta
 oldbeta1 <- oldbeta <- rep(1000,length(as.double(beta)))
  mu <- exp(xxxx %*% qbeta )
  S1 <- newY1 - mu
  SK1 <- 1
repeat{
 mu <- exp(xxxx %*% qbeta )
 res <- as.vector((newY1 - mu)/sqrt(mu))


  repeat{
    mu <- exp(xxxx %*% qbeta )
    S1 <- newY1 - mu

  G=matrix(0,2*dim(xxxx)[2],1)
  C=matrix(0,2*dim(xxxx)[2],2*dim(xxxx)[2])
  G1=matrix(0,2*dim(xxxx)[2],dim(xxxx)[2])
  for (i in 1:K) {
    M1 <- diag(n[i])
    M2 <-	matrix(0,n[i],n[i])
    M2[1,2]<-1
    for(o in 2:n[i]-1){
      M2[o,o+1] <-1
      M2[o,o-1] <-1
    }
    M2[n[i],n[i]-1]<-1
   
    W1=diag(Lambda[id==i])
    D1=diag(mu[id==i])%*%(as.matrix(xxxx[id==i,]))
  #  V1=MASS::ginv( sqrt(diag(mu[id==i])) )%*%M1%*%MASS::ginv( sqrt(diag(mu[id==i])) )
  #  V2=MASS::ginv( sqrt(diag(mu[id==i])) )%*%M2%*%MASS::ginv( sqrt(diag(mu[id==i])) )
    V1=( sqrt(diag(1/mu[id==i])) )%*%M1%*%( sqrt(diag(1/mu[id==i])) )
    V2=( sqrt(diag(1/mu[id==i])) )%*%M2%*%( sqrt(diag(1/mu[id==i])) )
    
    G_1=t(D1)%*%V1%*%W1%*%S1[id==i]
    G_1=rbind(G_1,t(D1)%*%V2%*%W1%*%S1[id==i])
    G=G+G_1
    C=C+G_1%*%t(G_1)

    VA1 <- -t(D1)%*%V1%*%W1%*%D1
    VA1=rbind( VA1,-t(D1)%*%V2%*%W1%*%D1)
    G1<-G1 +  VA1
  }

 # a <- C
 # storage.mode(a)<-"double"
 # lda<-as.integer(nrow(a))
 # npp<-as.integer(ncol(a))
  #print("a:")
  #print(a)
 # z <- .Fortran("finv", a=a, lda, npp)
 # Cinv <- z$a
  #print("a inv:")
  #print(arcinv)

  Q1=t(G1)%*%MASS::ginv(C)%*%G
  Q2=t(G1)%*%MASS::ginv(C)%*%G1
#  Q1=t(G1)%*%Cinv%*%G
# Q2=t(G1)%*%Cinv%*%G1

#  a <- Q2
 # storage.mode(a)<-"double"
 # lda<-as.integer(nrow(a))
 # npp<-as.integer(ncol(a))
  #print("a:")
  #print(a)
 # z <- .Fortran("finv", a=a, lda, npp)
 # Q2inv <- z$a
  #print("a inv:")
  #print(arcinv)
  
   qbeta=qbeta-MASS::ginv(Q2)%*%Q1
 # qbeta=qbeta-Q2inv%*%Q1

  if( all(abs( qbeta-oldbeta)<1e-6)) {
oldbeta <- qbeta
}else break
  }

  if( all(abs( qbeta-oldbeta1)<eps) ){
  oldbeta1=qbeta
}else  break
}
  list(beta = qbeta)
}
