#' Title  Variance estimate with sandwich formula based on the QIF approach
#'
#' @description Calculate the variance estimates using the sandwich formula based on the QIF approach.
#'
#' @param Time right censored data which is the follow up time.
#' @param Status the censoring indicator, normally 0 = event of interest happens, and 0 = censoring
#' @param X a matrix of covariates that may have effect on the failure times.
#' @param id a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.
#' @param beta estimation of covariate coefficients based on the QIF approach.
#' @param baseF estimation of baseline cumulative distribution function based on the Bernstein polynomial.
#' @param gam estimation of the  Bernstein polynomials' coefficients
#' @param N the degree of Bernstein polynomial used to approximate the baseline distribution function in model.
#' @param corstr a character string specifying the correlation structure. The following are permitted: \code{independence}, \code{exchangeable} and \code{AR1}.
#'
#' @export varqif
varqif <- function(Time, Status, X, id, beta, baseF,tau,gam,N,corstr) {
  Kn<- length(id)
  K <- length(unique(id))
  n <- as.vector(table(id))
  t2 <- Time
  c1 <- Status
  Y1=Status/baseF
  mu=exp(X%*%beta)
  S1=Y1-mu

  b.kN <- function(t,k,N){
    y <- choose(N,k) * (t/tau)^k * (1 - t/tau)^(N-k)
    return(y)
  }

  D_b.kN = function(t,k,N){
    result = (t/tau)^(k - 1) * (k * (1/tau)) * (1 - t/tau)^(N - k) - (t/tau)^k * ((1 - t/tau)^((N - k) - 1) * ((N - k) * (1/tau)))
    y = choose(N,k) * result
    return(y)
  }

  A <- matrix(0, nrow = length(Time), ncol = N+1)
  for(k in 0:N){
    A[,k+1] <- b.kN(Time,k,N)
  }

  A.derive <- matrix(0, nrow = length(Time), ncol = N+1)
  for(k in 0:N){
    A.derive[,k+1] <- D_b.kN(Time,k,N)
  }

  gamma0 <- gam
  exp_ga <- exp(gamma0)
  PPHI <- cumsum(exp_ga)/sum(exp_ga)
  DPPHI <- U_gam_gam <- matrix(0,N+1,N+1)
  DDPPHI <- array(0,c(N+1,N+1,N+1) )

  for(s in 1:(N+1)){
    for(j in 1:(N+1)){
      if(j >= s) {
        DPPHI[j,s] <- exp_ga[s]* (sum(exp_ga)- sum(exp_ga[1:j]))/((sum(exp_ga))^2)}else{
          DPPHI[j,s] <-  exp_ga[s]*  (- sum(exp_ga[1:j]) )/((sum(exp_ga))^2)  }
    }
  }

  for(l in 1:(N+1)){
    for(s in 1:(N+1)){
      for(k in 1:(N+1)){

        if(l >= s) {
          if(k!=s&l>=k){
            DDPPHI[l,s,k] <-  (  -2*sum(exp_ga)  *exp_ga[s]*exp_ga[k]*( sum(exp_ga) -  sum(exp_ga[1:l])  )      )       /((sum(exp_ga) )^4)
          }else if(k!=s&l<k){
            DDPPHI[l,s,k] <-  exp_ga[s] * ( exp_ga[k]* (sum(exp_ga))^2  -2*sum(exp_ga) * exp_ga[k]* ( sum(exp_ga) -  sum(exp_ga[1:l])  )    )             /((sum(exp_ga) )^4)
          }else{
            DDPPHI[l,s,k] <- ( ( sum(exp_ga) ^2)*( exp_ga[s] * ( sum(exp_ga) -  sum(exp_ga[1:l])  )   )     -2*(sum(exp_ga)) *exp_ga[s] * exp_ga[s] * ( sum(exp_ga) -  sum(exp_ga[1:l])  )      )   /((sum(exp_ga) )^4)
          }
        }else{

          if(k!=s&l>=k){
            DDPPHI[l,s,k] <-  (   exp_ga[s]  * (- exp_ga[k] )  *( sum(exp_ga)  ) ^2 - 2* (sum(exp_ga) )* exp_ga[s] * exp_ga[k] *( sum(exp_ga) -  sum(exp_ga[1:l])  )                )       /((sum(exp_ga) )^4)
          }else if(k!=s&l<k){
            DDPPHI[l,s,k] <-  (   - exp_ga[s] *(  -  sum(exp_ga[1:l])    ) * 2*(sum(exp_ga) )   * exp_ga[k]          )       /((sum(exp_ga) )^4)
          }else{
            DDPPHI[l,s,k] <-  (   exp_ga[s] *(- sum(exp_ga[1:l]) )  * (sum(exp_ga) ) ^2-2*sum(exp_ga) *exp_ga[s] *exp_ga[s] *(-sum(exp_ga[1:l]) )      )       /((sum(exp_ga) )^4)
          }
        }
      }
    }
  }

  for(s in 1:(N+1)){
    for(j in 1:(N+1)){

      U_gam_gam[s,j] <-  sum( Status * ( (  ( A.derive %*% as.double( DDPPHI[,s,j]) )   * ( A.derive %*% PPHI )      - (A.derive %*% as.double(DPPHI[,s]) ) *( A.derive %*% as.double(DPPHI[,j]) )     )  / (   A.derive %*% PPHI   )^2  ) ) - sum(  exp(X %*% beta) * (  A %*% as.double( DDPPHI[,s,j]  )   )  )

    }
  }

  U_gam_beta  <- matrix(0, (N+1), dim( X)[2])
  for(s in 1:(N+1)){
    for(j in 1:dim(X)[2] ){
      U_gam_beta[s,j] <-  - sum ( exp(X %*% beta) *  ( A %*% as.double( DPPHI[,s] )  )  *   as.double(X[,j]) )

    }
  }



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
    }else{
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
    G_1=t(D1)%*%V1%*%W1%*%S1[id==i]
    G11=-t(D1)%*%V1%*%W1%*%D1

    if(corstr!="independence"){
      V2=MASS::ginv( sqrt(diag(mu[id==i])) )%*%M2%*%MASS::ginv( sqrt(diag(mu[id==i])) )
      G_1=rbind(G_1,t(D1)%*%V2%*%W1%*%S1[id==i])
      G11=rbind(G11,-t(D1)%*%V2%*%W1%*%D1)
    }

    G=G+G_1
    C=C+G_1%*%t(G_1)
    VA1<-VA1+G11
  }


  BC1 <- matrix(0, dim(X)[2], N+1)
  BC2 <- matrix(0, dim(X)[2], N+1)

  for (s in 1:(N+1)) {

    Lamb <- A %*% as.double( DPPHI[,s] )
    elem <-  rep(0,dim(X)[2])
    for(i in 1:K){
      M1 <- diag(n[i])
      W1=diag(Lamb[id==i])
      D1=diag(mu[id==i])%*%diag(rep(1, n[i])) %*%(X[id==i,])
      V1=ginv( sqrt(diag(mu[id==i])) )%*%(M1)%*%ginv( sqrt(diag(mu[id==i])) )
      BC1[, s] <- elem + t(D1)%*%V1%*%W1%*%S1[id==i]
    }
  }

  for (s in 1:(N+1)) {

    Lamb <- A %*% as.double( DPPHI[,s] )
    elem <-  rep(0,dim(X)[2])
    for(i in 1:K){
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
      W1=diag(Lamb[id==i])
      D1=diag(mu[id==i])%*%X[id==i,]
      V2=ginv( sqrt(diag(mu[id==i])) )%*%(M2)%*%ginv( sqrt(diag(mu[id==i])) )
      BC2[, s] <- elem + t(D1)%*%V2%*%W1%*%S1[id==i]
    }
  }


  a <- C
  storage.mode(a)<-"double"
  lda<-as.integer(nrow(a))
  npp<-as.integer(ncol(a))
  z <- .Fortran("finv", a=a, lda, npp)
  Cinv <- z$a

  if(corstr=="independence"){
    BC=t( VA1)%*% Cinv%*%BC1
  }else{
    BC=t( VA1)%*% Cinv%*%rbind(BC1,BC2)
  }

  M22= -t(VA1)%*% Cinv%*%VA1
  M23=-BC
  M32=-U_gam_beta
  M33=-U_gam_gam

  M=rbind(cbind(M22,M23),cbind(M32,M33))
  fdv <- rep(0, dim(X)[2] + N+1)
  fdm <- matrix(0,dim(X)[2] + N+1, dim(X)[2] + N+1)
  for (i in 1:K) {
    z22 <- X[id == i, ]
    mu22 <- mu[id == i]
    mu22m <- diag(mu22, n[i], n[i])
    baseF01 <- baseF[id == i]
    G2 <- diag( baseF01)
    c22 <- c1[id == i]
    t22 <- t2[id == i]
    C2 <- (c22/baseF01) - mu22
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
      fdv0= t(mu22m%*%z22)%*%MASS::ginv(sqrt(mu22m))%*%M1%*%MASS::ginv(sqrt(mu22m))%*%G2%*%C2
    }else{
      fdv0= rbind(t(mu22m%*%z22)%*%MASS::ginv(sqrt(mu22m))%*%M1%*%MASS::ginv(sqrt(mu22m))%*%G2%*%C2 ,
                  t(mu22m%*%z22)%*%MASS::ginv(sqrt(mu22m))%*%M2%*%MASS::ginv(sqrt(mu22m))%*%G2%*%C2 )
    }


    fdv[1:dim(X)[2]] <-t(VA1)%*%Cinv%*%fdv0

    t22 <- t2[id == i]

    eqalpha <- rep(0, N+1)
    for (j in 1:(N+1)) {
      eqalpha[j] <-  sum(  c22 * ( A.derive[id==i,] %*%  DPPHI[,j] ) / ( A.derive[id==i,] %*%  PPHI ) ) - sum (  mu22 * A[id==i,] %*% DPPHI[,j] )
    }
    fdv[ (dim(X)[2] + 1):( dim(X)[2] + N+1)] <- eqalpha
    fdm1 <- fdv %*% t(fdv)
    fdm <- fdm + fdm1
  }
  vcm <- ginv(M) %*% fdm %*% t(ginv(M))
  var_beta <- as.matrix(vcm[1: dim(X)[2],1: dim(X)[2]]  )
  sd_beta <- sqrt(diag(vcm)[1: dim(X)[2]])
  var_F<-(vcm)[ (dim(X)[2] + 1):( dim(X)[2] + N+1),(dim(X)[2] + 1):( dim(X)[2] + N+1)]
  sd_F<- sqrt(diag(var_F))

  list(var_beta = var_beta, sd_beta = sd_beta,var_F = var_F, sd_F = sd_F)
}
