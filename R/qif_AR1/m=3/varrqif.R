#' Title Estimation of the variance of the covariate coefficients and baseline cumulative distribution function of QIF method
#'
#' @param Time right censored data which is the follow up time.
#' @param Status the censoring indicator, 1 = event of interest happens, and 0 = censoring.
#' @param X a matrix of covariates that that may have effect on the failure times.
#' @param id a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.
#' @param beta initial value of covariates coefficients.
#' @param gS the parameters in the baseline cumulative distribution function.
#' @param Lambda initial value of baseline cumulative distribution function.
#'
#' @return Estimation of the variance of the covariate coefficients and baseline cumulative distribution function of QIF method.
#' @export
varrqif <- function(Time, Status, X, id, beta, gS, Lambda) {
  Kn<- length(id)
  K <- length(unique(id))
  n <- as.vector(table(id))
  xxxx<-X
  t2 <- Time
  c1 <- Status
  t11 <- sort(Time)
  c11 <- Status[order(Time)]
  tt1 <- unique(t11[c11 == 1])
  kk <- length(table(t11[c11 == 1]))
  Lambda<-Lambda
  Y1=Status/Lambda
  be=as.matrix(beta)
  mu=exp(xxxx%*%be)
  S1=Y1-mu
  gg1 <- rep(1,Kn)


  BBC <- matrix(0, kk, dim( xxxx)[2])
  for (j in 1:dim( xxxx)[2]) {
    for (s in 1:(kk)) {
      BCm <- gS[s] * exp(xxxx[(c1 == 1) & (t2 == tt1[s]), ,drop = FALSE] %*% be)
      BBC[s, j] <- sum(exp(xxxx[(c1 == 1) & (t2 == tt1[s]), ,drop = FALSE] %*% be) * (exp(-BCm) + BCm * exp(-BCm) - 1)/((1 - exp(-BCm))^2) * xxxx[(c1 == 1) & (t2 == tt1[s]), j]) + sum(gg1[t2 >=
                                                                                                                                                                                              tt1[s]] * exp(xxxx[t2 >= tt1[s], ,drop = FALSE] %*% be) * xxxx[t2 >= tt1[s], j])
    }
  }
  CCC <- rep(0, (kk))
  for (s in 1:(kk)) {
    CCm <- gS[s] * exp(xxxx[(c1 == 1) & (t2 == tt1[s]), ,drop = FALSE] %*% be)
    CCC[s] <- sum(exp(2 * (xxxx[(c1 == 1) & (t2 == tt1[s]), ,drop = FALSE] %*% be) - CCm)/(1 - exp(-CCm))^2)
  }


  VA1 <- matrix(0, 3*dim(xxxx)[2], dim(xxxx)[2])
  newG11<-matrix(0, 3*dim(xxxx)[2], 1)
  for (v in 1:dim(xxxx)[2]) {

      for (i in 1:K) {
        M1 <- diag(n[i])
        M2 <-	matrix(0,n[i],n[i])
        M2[1,2]<-1
        for(o in 2:n[i]-1){
          M2[o,o+1] <-1
          M2[o,o-1] <-1
        }
        M2[n[i],n[i]-1]<-1
        M3<-	matrix(0,n[i],n[i])
        M3[1,1]<-1;M3[n[i],n[i]]<-1

        W1=diag(Lambda[id==i])
        D1=diag(mu[id==i])%*%diag(rep(1,n[i]))%*%(xxxx[id==i,])
        V1=MASS::ginv( sqrt(diag(mu[id==i])) )%*%M1%*%MASS::ginv( sqrt(diag(mu[id==i])) )
        V2=MASS::ginv( sqrt(diag(mu[id==i])) )%*%M2%*%MASS::ginv( sqrt(diag(mu[id==i])) )
        V3=MASS::ginv( sqrt(diag(mu[id==i])) )%*%M3%*%MASS::ginv( sqrt(diag(mu[id==i])) )

        G_1=t(D1)%*%V1%*%W1%*%S1[id==i]
        G_1=rbind(G_1,t(D1)%*%V2%*%W1%*%S1[id==i])
        G_1=rbind(G_1,t(D1)%*%V3%*%W1%*%S1[id==i])

        D11=t(D1)%*%diag(rep(1,n[i]))%*%diag(xxxx[id==i,v])
        newD1=diag(mu[id==i])%*%diag(rep(1,n[i]))%*%(xxxx[id==i,v])
        G11=D11%*%V1%*%W1%*%S1[id==i]-t(D1)%*%V1%*%W1%*%newD1-t(D1)%*%diag(xxxx[id==i,v])%*%V1%*%W1%*%S1[id==i]
        G11=rbind(G11,D11%*%V2%*%W1%*%S1[id==i]-t(D1)%*%V2%*%W1%*%newD1-t(D1)%*%diag(xxxx[id==i,v])%*%V2%*%W1%*%S1[id==i])
        G11=rbind(G11,D11%*%V3%*%W1%*%S1[id==i]-t(D1)%*%V3%*%W1%*%newD1-t(D1)%*%diag(xxxx[id==i,v])%*%V3%*%W1%*%S1[id==i])
        newG11<-newG11+G11
      }
      VA1[, v] <- newG11
      newG11<-matrix(0, 3*dim(xxxx)[2], 1)

  }


  G=matrix(0,3*dim(xxxx)[2],1)
  C=matrix(0,3*dim(xxxx)[2],3*dim(xxxx)[2])
  for (i in 1:K) {
    M1 <- diag(n[i])
    M2 <-	matrix(0,n[i],n[i])
    M2[1,2]<-1
    for(o in 2:n[i]-1){
      M2[o,o+1] <-1
      M2[o,o-1] <-1
    }
    M2[n[i],n[i]-1]<-1
    M3<-	matrix(0,n[i],n[i])
    M3[1,1]<-1;M3[n[i],n[i]]<-1

    W1=diag(Lambda[id==i])
    D1=diag(mu[id==i])%*%diag(rep(1,n[i]))%*%(xxxx[id==i,])
    V1=MASS::ginv( sqrt(diag(mu[id==i])) )%*%M1%*%MASS::ginv( sqrt(diag(mu[id==i])) )
    V2=MASS::ginv( sqrt(diag(mu[id==i])) )%*%M2%*%MASS::ginv( sqrt(diag(mu[id==i])) )
    V3=MASS::ginv( sqrt(diag(mu[id==i])) )%*%M3%*%MASS::ginv( sqrt(diag(mu[id==i])) )


    G_1=t(D1)%*%V1%*%W1%*%S1[id==i]
    G_1=rbind(G_1,t(D1)%*%V2%*%W1%*%S1[id==i])
    G_1=rbind(G_1,t(D1)%*%V3%*%W1%*%S1[id==i])
    G=G+G_1
    C=C+G_1%*%t(G_1)


  }


  BC1 <- matrix(0, dim(xxxx)[2], kk)
  BC2 <- matrix(0, dim(xxxx)[2], kk)
  BC3 <- matrix(0, dim(xxxx)[2], kk)
  for (r in 1:dim(xxxx)[2]) {
    for (s in 1:(kk)) {
      elem <- 0
      for (i in 1:K) {
        mu22 <- mu[id == i]
        xxx1 <- xxxx[id == i, r]
        t21 <- t2[id == i]

        M1 <- diag(n[i])

        for(j in 1:n[i])
        {
          if(t21[j]>=tt1[s])
            elem=elem+sum(xxx1*((mu22)^(1/2))*((mu22[j])^(-1/2))*M1[,j])*mu22[j]
        }
      }
      BC1[r, s] <- -elem
    }
  }

  for (r in 1:dim(xxxx)[2]) {
    for (s in 1:(kk)) {
      elem <- 0
      for (i in 1:K) {
        mu22 <- mu[id == i]
        xxx1 <- xxxx[id == i, r]
        t21 <- t2[id == i]

        M2 <-	matrix(0,n[i],n[i])
        M2[1,2]<-1
        for(o in 2:n[i]-1){
          M2[o,o+1] <-1
          M2[o,o-1] <-1
        }
        M2[n[i],n[i]-1]<-1

        for(j in 1:n[i])
        {
          if(t21[j]>=tt1[s])
            elem=elem+sum(xxx1*((mu22)^(1/2))*((mu22[j])^(-1/2))*M2[,j])*mu22[j]
        }
      }
      BC2[r, s] <- -elem
    }
  }

  for (r in 1:dim(xxxx)[2]) {
    for (s in 1:(kk)) {
      elem <- 0
      for (i in 1:K) {
        mu22 <- mu[id == i]
        xxx1 <- xxxx[id == i, r]
        t21 <- t2[id == i]

        M3<-	matrix(0,n[i],n[i])
        M3[1,1]<-1;M3[n[i],n[i]]<-1

        for(j in 1:n[i])
        {
          if(t21[j]>=tt1[s])
            elem=elem+sum(xxx1*((mu22)^(1/2))*((mu22[j])^(-1/2))*M3[,j])*mu22[j]
        }
      }
      BC3[r, s] <- -elem
    }
  }

  BC=-t( VA1)%*%MASS::ginv(C)%*%rbind(BC1,BC2,BC3)

  M22= -t(VA1)%*%MASS::ginv(C)%*%VA1
  M23=BC
  M32=t(t(BBC))
  M33=diag(CCC)

  M=rbind(cbind(M22,M23),cbind(M32,M33))

  fdv <- rep(0, dim(xxxx)[2] + kk)
  fdm <- matrix(0,dim(xxxx)[2] + kk, dim(xxxx)[2] + kk)
  for (i in 1:K) {
    z22 <- xxxx[id == i, ]
    mu22 <- mu[id == i]
    mu22m <- diag(mu22, n[i], n[i])
    Lambda01 <- Lambda[id == i]
    G2 <- diag( Lambda01)
    c22 <- c1[id == i]
    t22 <- t2[id == i]
    C2 <- (c22/Lambda01) - mu22
    M1 <- diag(n[i])
    M2 <-	matrix(0,n[i],n[i])
    M2[1,2]<-1
    for(o in 2:n[i]-1){
      M2[o,o+1] <-1
      M2[o,o-1] <-1
    }
    M2[n[i],n[i]-1]<-1
    M3<-	matrix(0,n[i],n[i])
    M3[1,1]<-1;M3[n[i],n[i]]<-1
    fdv0= rbind(t(mu22m%*%z22)%*%MASS::ginv(sqrt(mu22m))%*%M1%*%MASS::ginv(sqrt(mu22m))%*%G2%*%C2 ,
                t(mu22m%*%z22)%*%MASS::ginv(sqrt(mu22m))%*%M2%*%MASS::ginv(sqrt(mu22m))%*%G2%*%C2 ,
                t(mu22m%*%z22)%*%MASS::ginv(sqrt(mu22m))%*%M3%*%MASS::ginv(sqrt(mu22m))%*%G2%*%C2 )

    fdv[1:dim(xxxx)[2]] <-t(VA1)%*%MASS::ginv(C)%*%fdv0

    eqalpha <- rep(0, kk)
    for (j in 1:kk) {
      Aalpha <- (1:n[i])[t22 == tt1[j]]
      Balpha <- (1:n[i])[t22 >= tt1[j]]
      if (length(Balpha) == 0) {
        eqalpha[j] <- 0
      }
      if (length(Balpha) != 0 & length(Aalpha) == 0) {
        eqalpha[j] <- -sum(mu22[Balpha])
      } else eqalpha[j] <- sum(mu22[Aalpha]/(1 - exp(-gS[j] * mu22[Aalpha]))) - sum(mu22[Balpha])
    }
    fdv[ (dim(xxxx)[2] + 1):( dim(xxxx)[2] + kk)] <- eqalpha
    fdm1 <- fdv %*% t(fdv)
    fdm <- fdm + fdm1
  }
  vcm <- MASS::ginv(M) %*% fdm %*% t(MASS::ginv(M))
  var_beta <- diag(vcm)[1: dim(xxxx)[2]]
  sd_beta <- sqrt(var_beta)
  var_F<-(vcm)[ (dim(xxxx)[2] + 1):( dim(xxxx)[2] + kk),(dim(xxxx)[2] + 1):( dim(xxxx)[2] + kk)]
  sd_F<- sqrt(diag(var_F))

  list(var_beta = var_beta, sd_beta = sd_beta,var_F = var_F, sd_F = sd_F)
}
