#' Title  Estimate Variance of coefficients estimation from the GEE approach with a sandwich formula
#'
#' @description Calculate the variance estimates using the sandwich formula.
#'
#' @param Time right censored data which is the follow up time.
#' @param Status the censoring indicator, normally 0 = event of interest happens, and 0 = censoring
#' @param X a matrix of covariates that may have effect on the failure times.
#' @param id a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.
#' @param beta covariate coefficients estimation from GEE method.
#' @param corstr a character string specifying the correlation structure. The following are permitted: \code{independence}, \code{exchangeable} and \code{AR1}.
#'
#' @export vargee

vargee <- function(Time, Status, X, id, beta, corstr) {
  Kn<- length(id)
  K <- length(unique(id))
  n <- as.vector(table(id))

  mu <- exp(X %*% beta )

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
  baseF1[baseF1==0] <- 1e-8
  S1 <- (Status - mu*baseF)

  Y1 <- Status/baseF1

  be <- beta
  g1 <- rep(1,Kn)
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


  z2 <- X
  c2 <- Status
  mu2 <- exp(z2 %*% be)
  baseF0 <- baseF
  ABC1 <- rep(0, K)
  VA1 <- matrix(0, dim(z2)[2], dim(z2)[2])
  for (v in 1:dim(z2)[2]) {
    for (w1 in 1:dim(z2)[2]) {
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

        IQ1 <- solve(QC)
        B2 <- matrix(0, n[i], n[i])
        z22 <- matrix(z2[id == i, ], nrow = n[i], )
        A2 <- t(z22[, v])
        c22 <- c2[id == i]
        g22 <- g1[id == i]
        mu22 <- mu2[id == i]
        baseF22 <- baseF0[id == i]
        BB1 <- (mu22^(1/2)) %*% ((t(mu22))^(-1/2)) * IQ1
        for (s in 1:n[i]) {
          for (l in 1:n[i]) {
            B2[s, l] <- (1/2) * (z22[s, w1] - z22[l, w1]) * BB1[s, l]
          }
        }
        C2 <- c22 - mu22 * baseF22
        D2 <- BB1
        E2 <- z22[, w1] * mu22
        G2 <- diag(g22 * baseF22)
        ABC1[i] <- A2 %*% (B2 %*%  C2 - D2 %*% G2 %*% E2)
      }
      VA1[v, w1] <- sum(ABC1) * (pphi^(-1))
      ABC1 <- rep(0, K)
    }
  }


  M <- -VA1

  fdv <- rep(0, dim(X)[2])
  fdm <- matrix(0,dim(X)[2], dim(X)[2] )
  for (i in 1:K) {
    z22 <- X[id == i, ]
    mu22 <- mu[id == i]
    g11 <- g1[id == i]
    mu22m <- diag(mu22, n[i], n[i])
    baseF01 <- baseF[id == i]

    c22 <- c2[id == i]
    C2 <- (c22) - mu22 * baseF01
    if(corstr=="independence"){
      Q1 <- diag(n[i])
    }else  if(corstr=="exchangeable"){
      Q1=matrix(rho,n[i],n[i])
      diag(Q1)=1
    }else if(corstr=="AR1"){
      exponent <- abs(matrix(1:n[i] - 1, nrow = n[i], ncol = n[i], byrow = TRUE) - (1:n[i] - 1))
      Q1=rho^exponent
    }
    fdv <- t(mu22m %*% z22) %*% ginv(sqrt(mu22m) %*% Q1 %*% sqrt(mu22m) * pphi) %*% C2
    fdm1 <- fdv %*% t(fdv)
    fdm <- fdm + fdm1
  }

  vcm <- ginv(M) %*% fdm %*% t(ginv(M))
  var_beta <- diag( as.matrix(vcm) )
  sd_beta <- sqrt(var_beta)
  list(var_beta = var_beta, sd_beta = sd_beta)
}
