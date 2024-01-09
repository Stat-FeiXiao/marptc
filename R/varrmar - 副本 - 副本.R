#' Title Estimation of the variance of the covariate coefficients and baseline cumulative distribution function of GEE method
#'
#' @param Time right censored data which is the follow up time.
#' @param Status the censoring indicator, 1 = event of interest happens, and 0 = censoring.
#' @param X a matrix of covariates that that may have effect on the failure times.
#' @param id a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.
#' @param beta initial value of covariates coefficients.
#' @param gS the parameters in the baseline cumulative distribution function.
#' @param Lambda initial value of baseline cumulative distribution function.
#'
#' @return Estimation of the variance of the covariate coefficients and baseline cumulative distribution function of GEE method.
#' @export

varrmar <- function(Time, Status, X, id, beta, gS, Lambda,struc) {
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
  be=beta
  mu=exp(xxxx%*%be)
  S1=Y1-mu
  gg1 <- rep(1,Kn)
  res <- as.vector((Y1 - mu)/sqrt(mu))
  rres <- 0
  pphi <- sum(res^2)/(sum(n) - dim(X)[2])
  resm <- matrix(0, ncol = K, nrow = max(n))
  for (i in 1:K) {
    resm[1:n[i], i] <- res[id == i]
  }
  res <- resm
  res <- t(res)
  if(struc=="exchangeable"){
    for (i in 1:K) {
      if (n[i] == 1)
        rres <- rres + res[i, 1] else {
          for (j in 1:(n[i] - 1)) rres <- rres + res[i, j] * sum(res[i, (j + 1):n[i]])
        }
    }
    rho <- (pphi^(-1)) * rres/(sum(n * (n - 1))/2 - dim(X)[2] )
  }else{
    for (i in 1:K) {
      if (n[i] == 1)
        rres <- rres + res[i, 1] else {
          for (j in 1:(n[i] - 1)) rres <- rres + res[i, j] * res[i, j+1]
        }
    }
    rho <- (pphi^(-1)) * rres/(sum((n - 1))/2 - dim(X)[2] )
    
  }


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

  BC <- matrix(0, dim(xxxx)[2], kk)
  for (r in 1:dim(xxxx)[2]) {
    for (s in 1:(kk)) {
      elem <- 0
      for (i in 1:K) {
        mu22 <- mu[id == i]
        xxx1 <- xxxx[id == i, r]
        t21 <- t2[id == i]

        if(struc=="exchangeable"){
          QC=matrix(rho,n[i],n[i])
          diag(QC)=1
        }else if(struc=="AR1"){
          exponent <- abs(matrix(1:n[i] - 1, nrow = n[i], ncol = n[i], byrow = TRUE) - (1:n[i] - 1))
          QC=rho^exponent
        }

        M1=pphi^(-1)*ginv(QC)
        for(j in 1:n[i])
        {
        
          elem=elem-sum(xxx1*((mu22)^(1/2))*((mu22[j])^(-1/2))*M1[,j])*sum(as.numeric(tt1<=t21[j])*gS)/(( sum(gS))^(2))*mu22[j]
          if(t21[j]>=tt1[s])
            elem=elem+sum(xxx1*((mu22)^(1/2))*((mu22[j])^(-1/2))*M1[,j])*sum(as.numeric(tt1>t21[j])*gS)/(( sum(gS))^(2))*mu22[j]
          
        }
      }
      BC[r, s] <- elem
    }
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

  M22=-VA1
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
    C2 <- (c22/Lambda01) - mu22
    if(struc=="exchangeable"){
      Q1=matrix(rho,n[i],n[i])
      diag(Q1)=1
    }else if(struc=="AR1"){
      exponent <- abs(matrix(1:n[i] - 1, nrow = n[i], ncol = n[i], byrow = TRUE) - (1:n[i] - 1))
      Q1=rho^exponent
    }
    fdv[1:dim(xxxx)[2]] <- t(mu22m %*% z22) %*% ginv(sqrt(mu22m) %*% Q1 %*% sqrt(mu22m) * pphi) %*% G2 %*% C2
    t22 <- t2[id == i]
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
  vcm <- ginv(M) %*% fdm %*% t(ginv(M))
  var_beta <- as.matrix(vcm[1: dim(xxxx)[2],1: dim(xxxx)[2]]  ) 
  sd_beta <- sqrt(diag(vcm)[1: dim(xxxx)[2]])
  var_F<-(vcm)[ (dim(xxxx)[2] + 1):( dim(xxxx)[2] + kk),(dim(xxxx)[2] + 1):( dim(xxxx)[2] + kk)]
  sd_F<- sqrt(diag(var_F))
  list(var_beta = var_beta, sd_beta = sd_beta,var_F = var_F, sd_F = sd_F)
}
