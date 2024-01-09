
vartso <- function(Time, Status, X, id, beta, gS, Lambda) {
  Kn<- length(id)
  K <- length(unique(id))
  n <- as.vector(table(id))
  xxxx<-xxx<-as.matrix(X)
  t2 <- Time
  c1 <- Status
  t11 <- sort(Time)
  c11 <- Status[order(Time)]
  tt1 <- unique(t11[c11 == 1])
  x111 <- as.matrix(xxx[order(t2), ])
  kk <- length(table(t11[c11 == 1]))
  dd <- as.matrix(table(t11[c11 == 1]))
  Lambda<-Lambda
  Y1=Status/Lambda
  be=beta
  mu=exp(xxxx%*%be)
  S1=Y1-mu
  gg1 <- rep(1,Kn) 

  
  BBC <- matrix(0, kk, dim( xxxx)[2])
  for (j in 1:dim( xxxx)[2]) {
    for (s in 1:(kk)) {
      BCm <- gS[s] * exp(xxxx[(c1 == 1) & (t2 == tt1[s]), ,drop = FALSE] %*% be)  
      BBC[s, j] <- sum(exp(xxxx[(c1 == 1) & (t2 == tt1[s]), ,drop = FALSE] %*% be) * (exp(-BCm) + BCm * exp(-BCm) - 1)/((1 - exp(-BCm))^2) * xxxx[(c1 == 1) & (t2 == tt1[s]), j]) + sum(gg1[t2 >= tt1[s]] * exp(xxxx[t2 >= tt1[s], ,drop = FALSE] %*% be) * xxxx[t2 >= tt1[s], j])
    }
  }
  
  CCC <- rep(0, (kk))
  for (s in 1:(kk)) {
    CCm <- gS[s] * exp(xxxx[(c1 == 1) & (t2 == tt1[s]), ,drop = FALSE] %*% be)
    CCC[s] <- sum(exp(2 * (xxxx[(c1 == 1) & (t2 == tt1[s]), ,drop = FALSE] %*% be) - CCm)/(1 - exp(-CCm))^2)
  }
  
  fun11=function(B1){
    fi2<-matrix(0,dim(xxx)[2],dim(xxx)[2])
    for(s in 1:dim(xxx)[2]){
      for(l in 1:dim(xxx)[2]){
        zegaB=exp( xxx%*%B1 )
        dataorderallB=data.frame(id,c1,t2,zegaB)
        zega1B<-dataorderallB$zegaB
        tze1B<-dataorderallB$t2
        zegasiB=rep(0,kk)
        for (i in 1:kk) {
          zegasiB[i]=sum(zega1B[tze1B>=tt1[i]])
        }
        
        
        zegaBx1=exp( xxx%*%B1 )*xxx[,s]
        dataorderallBx1=data.frame(id,c1,t2,zegaBx1)
        zega1Bx1<-dataorderallBx1$zegaBx1
        tze1Bx1<-dataorderallBx1$t2
        zegasiBx1=rep(0,kk)
        for (i in 1:kk) {
          zegasiBx1[i]=sum(zega1Bx1[tze1Bx1>=tt1[i]])
        }
        zegaBx2=exp( xxx%*%B1 )*xxx[,l]
        dataorderallBx2=data.frame(id,c1,t2,zegaBx2)
        zega1Bx2<-dataorderallBx2$zegaBx2
        tze1Bx2<-dataorderallBx2$t2
        zegasiBx2=rep(0,kk)
        for (i in 1:kk) {
          zegasiBx2[i]=sum(zega1Bx2[tze1Bx2>=tt1[i]])
        }
        
        
        zegaBxx=exp( xxx%*%B1 )*xxx[,s]*xxx[,l]
        dataorderallBxx=data.frame(id,c1,t2,zegaBxx)
        zega1Bxx<-dataorderallBxx$zegaBxx
        tze1Bxx<-dataorderallBxx$t2
        zegasiBxx=rep(0,kk)
        for (i in 1:kk) {
          zegasiBxx[i]=sum(zega1Bxx[tze1Bxx>=tt1[i]])
        }
        
        fun1=-zegasiBxx/zegasiB+(zegasiBx1*zegasiBx2)/zegasiB^2
        too <- t2[c1==1]
        a  <- 0
        for(o in 1:length(t2[c1==1])){
          b<- fun1[too[o]==tt1]
          a<-a+b
        }
        
        
        fi2[s,l]<-a
      }
    }
    
    return(fi2)
  }
  
  VA1 <- - fun11( beta)
  M22=VA1
  M23=matrix(0,dim(X)[2],kk)
  M32=t(t(BBC))
  M33=diag(CCC)
  
  M=rbind(cbind(M22,M23),cbind(M32,M33))
  
  
  Ui=function(ii,B1){
    fi1<-rep(0,dim(xxx)[2])
    for(s in 1:dim(xxx)[2]){
      zegaB=exp( xxx%*%B1 )
      dataorderallB=data.frame(id,c1,t2,zegaB)
      zega1B<-dataorderallB$zegaB
      tze1B<-dataorderallB$t2
      zegasiB=rep(0,kk)
      for (i in 1:kk) {
        zegasiB[i]=sum(zega1B[tze1B>=tt1[i]])
      }
      
      zegaBx=exp( xxx%*%B1 )*xxx[,s]
      dataorderallBx=data.frame(id,c1,t2,zegaBx)
      zega1Bx<-dataorderallBx$zegaBx
      tze1Bx<-dataorderallBx$t2
      zegasiBx=rep(0,kk)
      for (i in 1:kk) {
        zegasiBx[i]=sum(zega1Bx[tze1Bx>=tt1[i]])
      }
      if(all(c1[id==ii]==0)){
        a<-0}else{
          too <- t2[c1==1 &id==ii]
          a  <- 0
          
          for(o in 1:length(t2[c1==1&id==ii])){
            b<- (zegasiBx/zegasiB)[too[o]==tt1]
            a<-a+b
          }
        }
      id11<- id[order(t2)]
      
      #  fun1=-zegasiBx/zegasiB
      fun2=(c11*x111[,s])[id11==ii]
      fi1[s]<-sum(fun2)-a
    }
    return(fi1)
  }
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
   
    fdv[1:dim(xxxx)[2]] <-  Ui(i,beta)
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
  var_beta <- diag(vcm)[1: dim(xxxx)[2]]
  sd_beta <- sqrt(var_beta)
  var_F<-(vcm)[ (dim(xxxx)[2] + 1):( dim(xxxx)[2] + kk),(dim(xxxx)[2] + 1):( dim(xxxx)[2] + kk)]
  sd_F<- sqrt(diag(var_F)) 
  
  list(var_beta = var_beta, sd_beta = sd_beta,var_F = var_F, sd_F = sd_F)
}
