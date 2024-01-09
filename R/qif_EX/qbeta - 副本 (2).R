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
    pphi <- sum(res^2)/(sum(n) - dim(X)[2] ) 
    
    repeat{
      mu <- exp(xxxx %*% qbeta )
      S1 <- newY1 - mu
      
      G=matrix(0,2*dim(xxxx)[2],1)
      C=matrix(0,2*dim(xxxx)[2],2*dim(xxxx)[2])
      G1=matrix(0,2*dim(xxxx)[2],dim(xxxx)[2])
      dotC <- array(0,c(2*dim(xxxx)[2],2*dim(xxxx)[2],dim(X)[2]))
      for (i in 1:K) {
        M1 <- diag(n[i])
        M2 <-matrix(1,n[i],n[i])
        diag(M2) <- 0
        
        W1=diag(Lambda[id==i])
        D1=diag(mu[id==i])%*%(as.matrix(xxxx[id==i,]))
        
        V1=MASS::ginv( sqrt(diag(mu[id==i])) )%*%M1%*%MASS::ginv( sqrt(diag(mu[id==i])) )
        V2=MASS::ginv( sqrt(diag(mu[id==i])) )%*%M2%*%MASS::ginv( sqrt(diag(mu[id==i])) )
        #   V1=( sqrt(diag(1/mu[id==i])) )%*%M1%*%( sqrt(diag(1/mu[id==i])) )
        #   V2=( sqrt(diag(1/mu[id==i])) )%*%M2%*%( sqrt(diag(1/mu[id==i])) )
        
        G_1=t(D1)%*%V1%*%W1%*%S1[id==i]
        G_1=rbind(G_1,t(D1)%*%V2%*%W1%*%S1[id==i])
        G=G+G_1
        C=C+G_1%*%t(G_1)
      }
      ####
      
      ABC1 <- rep(0, K)
      VA1 <- matrix(0, dim(X)[2], dim(X)[2])
      for (v in 1:dim(X)[2]) {
        for (w1 in 1:dim(X)[2]) {
          for (i in 1:K) {
            M1 <- diag(n[i])
            IQ1 <- M1
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
          VA1[v, w1] <- sum(ABC1)
          ABC1 <- rep(0, K)
        }
      }
      ABC1 <- rep(0, K)
      VA2 <- matrix(0, dim(X)[2], dim(X)[2])
      for (v in 1:dim(X)[2]) {
        for (w1 in 1:dim(X)[2]) {
          for (i in 1:K) {
            M2 <-matrix(1,n[i],n[i])
            diag(M2) <- 0
            IQ1 <- M2
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
          VA2[v, w1] <- sum(ABC1)
          ABC1 <- rep(0, K)
        }
      }
      G1 <- rbind(VA1,VA2)
      ##
      
 
      
      for (i in 1:K) {
        M1 <- diag(n[i])
        M2 <-	matrix(1,n[i],n[i])
        diag(M2) <- 0
        
        W1=diag(Lambda[id==i])
        D1=diag(mu[id==i])%*%(as.matrix(xxxx[id==i,]))
        
        V1=MASS::ginv( sqrt(diag(mu[id==i])) )%*%M1%*%MASS::ginv( sqrt(diag(mu[id==i])) )
        V2=MASS::ginv( sqrt(diag(mu[id==i])) )%*%M2%*%MASS::ginv( sqrt(diag(mu[id==i])) )
        #   V1=( sqrt(diag(1/mu[id==i])) )%*%M1%*%( sqrt(diag(1/mu[id==i])) )
        #   V2=( sqrt(diag(1/mu[id==i])) )%*%M2%*%( sqrt(diag(1/mu[id==i])) )
        
        G_1=t(D1)%*%V1%*%W1%*%S1[id==i]
        G_1=rbind(G_1,t(D1)%*%V2%*%W1%*%S1[id==i])
        for(mmm in 1:dim(X)[2]){
          VA1 <- matrix(0, dim(X)[2], 1)
          for (v in 1:dim(X)[2]) {
            
            IQ1 <- M1
            B2 <- matrix(0, n[i], n[i])
            z22 <- matrix(X[id == i, ], nrow = n[i], )
            A2 <- t(z22[, v])
            c22 <- Status[id == i]
            mu22 <- mu[id == i]
            Lambda22 <- Lambda[id == i]
            BB1 <- (mu22^(1/2)) %*% ((t(mu22))^(-1/2)) * IQ1
            for (s in 1:n[i]) {
              for (l in 1:n[i]) {
                B2[s, l] <- (1/2) * (z22[s, mmm] - z22[l, mmm]) * BB1[s, l]
              }
            }
            C2 <- c22/Lambda22 - mu22
            D2 <- BB1
            E2 <- z22[, mmm] * mu22
            G2 <- diag(Lambda22)  
            VA1[v, 1] <-VA2[v, 1]  +  A2 %*% (B2 %*% G2 %*% C2 - D2 %*% G2 %*% E2)
            
          }
          
          VA2 <- matrix(0, dim(X)[2], 1)
          
          for (v in 1:dim(X)[2]) {
            
            IQ1 <- M2
            B2 <- matrix(0, n[i], n[i])
            z22 <- matrix(X[id == i, ], nrow = n[i], )
            A2 <- t(z22[, v])
            c22 <- Status[id == i]
            mu22 <- mu[id == i]
            Lambda22 <- Lambda[id == i]
            BB1 <- (mu22^(1/2)) %*% ((t(mu22))^(-1/2)) * IQ1
            for (s in 1:n[i]) {
              for (l in 1:n[i]) {
                B2[s, l] <- (1/2) * (z22[s, mmm] - z22[l, mmm]) * BB1[s, l]
              }
            }
            C2 <- c22/Lambda22 - mu22
            D2 <- BB1
            E2 <- z22[, mmm] * mu22
            G2 <- diag(Lambda22)  
            VA2[v, 1] <- VA2[v, 1]+ A2 %*% (B2 %*% G2 %*% C2 - D2 %*% G2 %*% E2) 
          }
          VA1dot <- rbind(VA1,VA2)
          
          dotC[,,mmm] <- dotC[,,mmm] + VA1dot%*%t( G_1 )+ G_1 %*% t( VA1dot )
        }
      }
      
      
      
      
      ##
      
      ####
      # a <- C
      # storage.mode(a)<-"double"
      #  lda<-as.integer(nrow(a))
      # npp<-as.integer(ncol(a))
      #print("a:")
      #print(a)
      # z <- .Fortran("finv", a=a, lda, npp)
      # Cinv <- z$a
      #print("a inv:")
      #print(arcinv)
      Q1dott <- rep(0,dim(X)[2])
      for(iii in 1:dim(X)[2]){
        Q1dott[iii] <- t(G)%*%MASS::ginv(C)%*%dotC[,,iii]%*%MASS::ginv(C)%*%G 
      }
      Q1 <- 2*t(G1)%*%MASS::ginv(C)%*%G #-Q1dott  # Q1=t(G1)%*%MASS::ginv(C)%*%G
      
      Q2 <- 2*t(G1)%*%MASS::ginv(C)%*%G1 
      #Q1=t(G1)%*%Cinv%*%G
      # Q2=t(G1)%*%Cinv%*%G1
      
      # a <- Q2
      #  storage.mode(a)<-"double"
      #  lda<-as.integer(nrow(a))
      #  npp<-as.integer(ncol(a))
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
