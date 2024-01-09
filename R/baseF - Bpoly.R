baseF<- function(Time, Status, X, beta,N) {
  Kn <- length(Time)
  t2 <- Time
  t11 <- sort(Time)
  c11 <- Status[order(Time)]
  x111 <- as.matrix(X[order(Time), ])
  tt1 <- unique(t11[c11 == 1])
  kk <- length(table(t11[c11 == 1]))
  dd <- as.matrix(table(t11[c11 == 1]))
  tau <- max( Time[Status==1] ) - 0.0001  #max(Time) + 0.1 # max( Time[Status==1] ) + 0.1 # 

  b.kN <- function(t,k,N){  #k=1,...,(m+1)
  # t: observe surval time, a vector
  # k: must be 1 dimention
  y <- choose(N,k) * (t/tau)^k * (1 - t/tau)^(N-k)
  return(y)
  }

  D_b.kN = function(t,k,N){
  # result = D(expression((t/tau)^k * (1 - t/tau)^(N-k)), 't')  
  result = (t/tau)^(k - 1) * (k * (1/tau)) * (1 - t/tau)^(N - k) - (t/tau)^k * ((1 - t/tau)^((N - k) - 1) * ((N - k) * (1/tau)))
  y = choose(N,k) * result
  return(y)
  }
  f.t.star = function(t.star, gamma){
  ## F(t*.m)
  #  t is m * 1
  A.derive <- matrix(0, nrow = length(t.star), ncol = N+1)
  for(k in 0:N){
    A.derive[,k+1] <- D_b.kN(t.star,k,N) 
  }
  
exp_ga <- exp(gamma)
 PPHI <- cumsum(exp_ga)/sum(exp_ga)

  f.t = A.derive %*% PPHI 
  f.t = as.vector(f.t)
  return(f.t)
}


  F.t.star = function(t.star, gamma){
  ## F(t*.m)
  #  t is m * 1
  A <- matrix(0, nrow = length(t.star), ncol = N+1)
  for(k in 0:N){
    A[,k+1] <- b.kN(t.star,k,N) 
   }
exp_ga <- exp(gamma)
 PPHI <- cumsum(exp_ga)/sum(exp_ga)

  F.t = A %*% PPHI 
  F.t = as.vector(F.t)
  return(F.t)
}

 A <- matrix(0, nrow = length(Time), ncol = N+1)
  for(k in 0:N){
    A[,k+1] <- b.kN(Time,k,N) 
  }
  
  A.derive <- matrix(0, nrow = length(Time), ncol = N+1)
  for(k in 0:N){
    A.derive[,k+1] <- D_b.kN(Time,k,N) 
  }
  

Loglike<-function(gam){  #gam N+1 *1

    gamma0 <- gam
    exp_ga <- exp(gamma0)
    PPHI <- cumsum(exp_ga)/sum(exp_ga)

    F.Yi = A %*% PPHI 

    f.Yi = A.derive %*% PPHI  
   Like = Status * log(f.Yi) - (exp(X %*% beta) * F.Yi) 
    Like= sum(Like)
    return(Like)
}

Loglike1<-function(gam){  #gam N+1 *1

    gamma0 <- gam
    exp_ga <- exp(gamma0)
    PPHI <- cumsum(exp_ga)/sum(exp_ga)

    F.Yi = A %*% PPHI 

    f.Yi = A.derive %*% PPHI  
   Like = Status * log(f.Yi)+Status*(X %*% beta) - (exp(X %*% beta) * F.Yi) 
    Like= sum(Like)
    return(Like)
}

Eqtion = function(gam){
  gamma0 <- gam
  exp_ga <- exp(gamma0) 
  PPHI <- cumsum(exp_ga)/sum(exp_ga)
  DPPHI <- eql <- rep(0,(N+1))
  f.Yi = A.derive %*% PPHI  
  
  for(s in 1:(N+1)){
    for(j in 1:(N+1)){
      if(j >= s) {
        DPPHI[j] <- exp_ga[s]* (sum(exp_ga)- sum(exp_ga[1:j]))/((sum(exp_ga))^2)}else{
          DPPHI[j] <-  exp_ga[s]*  (- sum(exp_ga[1:j]) )/((sum(exp_ga))^2)
        }
    }
    DF = A %*% DPPHI
    
    Df = A.derive %*% DPPHI
    
    
    eql [s]<-sum( Status * (1/f.Yi) * Df - exp(X %*% beta) * DF)
    
    
  }
  #   temp.gamma = Status * sweep(A.derive, 1, f.Yi, '/') - sweep(A, 1, exp(X %*% beta), '*')
  #   f= apply(temp.gamma, 2, sum)
  
  
  return(eql)
}
  
 # C1<-matrix(0,N+1,N+1)
 # diag(C1)<-1

 # for(j1 in 1:N)
 # {
 #   C1[j1+1,j1]<--1
 # }
 # C1<-rbind(C1,c(rep(0,N),-1),c(rep(0,N),1))
 # C<-cbind(matrix(0,N+3,3),C1)
 # D<-c(rep(0,N+1),1.01,-0.99)


re_esti<-maxLik::maxLik(Loglike, start = rep(0.1,N+1), method = 'BFGS', grad = Eqtion)$estimate #
                        
                        re_esti<-as.double(re_esti)
     
                        gSS3 <-  F.t.star(Time, re_esti)
fhat <- f.t.star(Time, re_esti)
logL <- Loglike1(re_esti)
qic <- -2*logL + 2*(dim(X)[2]+length( re_esti) )
list(gSS3 = gSS3,fhat=fhat,logL=logL,qic=qic)
}