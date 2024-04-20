#' Title  Variance estimate with sandwich formula based on the GEE approach
#'
#' @description Calculate the variance estimates using the sandwich formula based on the GEE approach.
#'
#' @param Time right censored data which is the follow up time.
#' @param Status the censoring indicator, normally 0 = event of interest happens, and 0 = censoring
#' @param X a matrix of covariates that may have effect on the failure times.
#' @param id a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.
#' @param beta estimation of covariate coefficients based on the GEE approach.
#' @param baseF estimation of baseline cumulative distribution function based on the Bernstein polynomial.
#' @param gam estimation of the  Bernstein polynomials' coefficients.
#' @param N the degree of Bernstein polynomial used to approximate the baseline distribution function in model.
#' @param corstr a character string specifying the correlation structure. The following are permitted: \code{independence}, \code{exchangeable} and \code{AR1}.
#'
#' @export vargee
vargee <- function(Time, Status, X, id, beta, baseF,tau,gam,N,corstr) {
  Kn<- length(id)
  K <- length(unique(id))
  n <- as.vector(table(id))
  t2 <- Time
  c1 <- Status
  Y1 <- Status/baseF
  mu=exp(X%*%beta)
  S1=Y1-mu
  res <- as.vector((Y1 - mu)/sqrt(mu))
  rres <- 0
  pphi <- sum(res^2)/(sum(n) - dim(X)[2])

  resm <- matrix(0, ncol = K, nrow = max(n))
  for (i in 1:K) {
    resm[1:n[i], i] <- res[id == i]
  }
  res <- resm
  res <- t(res)
  if(corstr=="exchangeable"){
    for (i in 1:K) {
      if (n[i] == 1)
        rres <- rres + res[i, 1] else {
          for (j in 1:(n[i] - 1)) rres <- rres + res[i, j] * sum(res[i, (j + 1):n[i]])
        }
    }
    rho <- (pphi^(-1)) * rres/(sum(n * (n - 1))/2 -dim(X)[2] )
  }else if(corstr=="AR1"){
    for (i in 1:K) {
      if (n[i] == 1)
        rres <- rres + res[i, 1] else {
          for (j in 1:(n[i] - 1)) rres <- rres + res[i, j] * res[i, j+1]
        }
    }
    rho <- (pphi^(-1)) * rres/(sum((n - 1)) - dim(X)[2]  )

  }

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


  U_gam_beta  <- matrix(0, (N+1), dim(X)[2])
  for(s in 1:(N+1)){
    for(j in 1:dim(X)[2] ){
      U_gam_beta[s,j] <-  - sum ( exp(X %*% beta) *  ( A %*% as.double( DPPHI[,s] )  )  *   as.double(X[,j]) )

    }
  }




  BC <- matrix(0, dim(X)[2], N+1)

  for (s in 1:(N+1)) {

    Lamb <- A %*% as.double( DPPHI[,s] )
    elem <-  rep(0,dim(X)[2])
    for(i in 1:K){
      if(corstr=="independence"){
        QC <- diag(n[i])
      }else if(corstr=="exchangeable"){
        QC=matrix(rho,n[i],n[i])
        diag(QC)=1
      }else if(corstr=="AR1"){
        exponent <- abs(matrix(1:n[i] - 1, nrow = n[i], ncol = n[i], byrow = TRUE) - (1:n[i] - 1))
        QC=rho^exponent
      }
      W1=diag(Lamb[id==i])
      D1=diag(mu[id==i])%*%diag(rep(1, n[i])) %*%(X[id==i,])
      V1=ginv( sqrt(diag(mu[id==i])) )%*%(pphi^(-1)*ginv(QC))%*%ginv( sqrt(diag(mu[id==i])) )    #
      BC[, s] <- elem + t(D1)%*%V1%*%W1%*%S1[id==i]
    }
  }


  VA1 <- matrix(0, dim(X)[2], dim(X)[2])


  for (i in 1:K) {
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
    V1=MASS::ginv( sqrt(diag(mu[id==i])) )%*%(pphi^(-1)*MASS::ginv(QC))%*%MASS::ginv( sqrt(diag(mu[id==i])) )


    G11=-t(D1)%*%V1%*%W1%*%D1

    VA1<-VA1+G11
  }



  M22=-VA1
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
    C2 <- (c22/baseF01) - mu22
    if(corstr=="independence"){
      Q1 <- diag(n[i])
    }else  if(corstr=="exchangeable"){
      Q1=matrix(rho,n[i],n[i])
      diag(Q1)=1
    }else if(corstr=="AR1"){
      exponent <- abs(matrix(1:n[i] - 1, nrow = n[i], ncol = n[i], byrow = TRUE) - (1:n[i] - 1))
      Q1=rho^exponent
    }
    fdv[1:dim(X)[2]] <- t(mu22m %*% z22) %*% ginv(sqrt(mu22m) %*% Q1 %*% sqrt(mu22m) * pphi) %*% G2 %*% C2

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
