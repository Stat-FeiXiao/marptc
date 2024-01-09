baseF<- function(Time, Status, X, beta,alpha_esti=NULL) {

  Kn <- length(Time)
  t2 <- Time
  CTime <- Time[Status==1]
  OrdTime <- t11 <- sort(Time)
  c11 <- Status[order(Time)]
  ordX <- x111 <- as.matrix(X[order(Time), ])
  OrdCTime <- t11[c11 == 1]
  tt1 <- unique(OrdCTime)
  kk <- length(table(t11[c11 == 1]))
  dd <- as.matrix(table(t11[c11 == 1]))
  tau <- max(tt1)+0.001
  print(all(as.double(dd) ==1    ))
  F <- function(alp){
    gSS <- rep(0, kk-1)
    FF <- rep(0,Kn)
    gSS[1] <- 1 - exp( - alp[1] )
    for (i in 1:(kk - 2)) {
      gSS[i + 1] <- 1 - exp( - sum(alp[1:(i+1)] ) )
    }
    gS=c(gSS[1],gSS[2:(kk-1)]-gSS[1:(kk-2)], 1- gSS[kk-1])
    for (i in 1:(Kn)) {
      kk1 <- 1
      
      if (t2[i] < tt1[1]) {
        FF[i] <- 1e-08
      } else {
        if (t2[i] >= tt1[kk]) {
          FF[i] <-  1
        } else {
          repeat {
            if(t2[i]>=tt1[kk1]) kk1=kk1+1
            else break
          }
          {
            FF[i] <- gSS[kk1 - 1]
          }
        }
      }
    }
    ff = rep(gS,dd)
    list(FF=FF,ff=ff)  #FF = F(Time)  ; ff = f(  order(Time[cens==1] ) ) = f(ordTime[c11==1])
  }
  
  Loglike<-function(alpha)  
  {
    
    alpha0 <- alpha
    
    F.Yi = F( alpha0)$FF
    f.Yi = F( alpha0)$ff

 #   Like =sum( log(exp(ordX[c11==1,,drop=FALSE] %*% beta) * f.Yi) )- sum((exp(X %*% beta) * F.Yi) )
 Like =sum( log( f.Yi) )- sum((exp(X %*% beta) * F.Yi) )
    return(Like)
    
  }

  ##
  F1 <-  function(alpha){
    alp <- alpha
    gSS <- rep(0, kk-1)
    FF <- rep(0,Kn)
    gSS[1] <-   exp( - alp[1] )
    for (i in 1:(kk - 2)) {
      gSS[i + 1] <-   exp( - sum(alp[1:(i+1)] ) )
    }
    for (i in 1:(Kn)) {
      kk1 <- 1
      
      if (t2[i] < tt1[1]) {
        FF[i] <- 0
      } else {
        if (t2[i] >= tt1[kk]) {
          FF[i] <-  0
        } else {
          repeat {
            if(t2[i]>=tt1[kk1]) kk1=kk1+1
            else break
          }
          {
            FF[i] <- gSS[kk1 - 1]
          }
        }
      }
    }
    return(FF)
  }
  ordFfun <- function(alpha){
    alpha0 <- alpha
    F1 <-  F1( alpha0)
    ordF1 <- F1[order(Time)]
    return(ordF1)
  }
  ordf1fun <- function(alpha){
    alp <- alpha
    gSS <- rep(0, kk-2)
gSS [1] <- 0
    for (i in 2:(kk - 1)) {
      gSS[i] <-  -  exp( - sum(alp[1:(i-1)] ) ) + exp( - sum(alp[1:(i)] ) )
    }
    return(ordf=gSS)
  }


  ##
  Eqtion = function(alpha){
    ordF <- ordFfun(alpha)
    ff = F( alpha)$ff
    ordf1 <- ordf1fun(alpha)
     
    f = rep(0,kk-1)
    
   

#    for(p in 1:(kk-2)){        
 #     f[p] = sum( ordX[OrdTime ==max(tt1) &c11==1,,drop=FALSE]  %*% beta * 1/(ff[OrdCTime ==max(tt1)  ] ) * (-exp(-sum(alpha))) )      +
#      sum( ordX[OrdTime <max(tt1) &OrdTime>tt1[p]&c11==1,,drop=FALSE]  %*% beta * 1/(ff[OrdCTime <max(tt1) & OrdCTime>tt1[p]]) * ordf1 [(p+1):(kk-1)] )      +
#       sum( ordX[OrdTime==tt1[p]&c11==1,,drop=FALSE] %*% beta * 1/(ff[OrdCTime==tt1[p]]) *exp(-sum(alpha[1:p])) ) -  
#      sum( exp(ordX %*% beta) [ OrdTime>= tt1[p] ] *   ordF[ OrdTime>= tt1[p] ]    )
#    }
         
 #     f[kk-1] = sum( ordX[OrdTime ==max(tt1) &c11==1,,drop=FALSE]  %*% beta * 1/(ff[OrdCTime ==max(tt1)  ] ) * (-exp(-sum(alpha))) )   +   
#      sum( ordX[OrdTime==tt1[p]&c11==1,,drop=FALSE] %*% beta * 1/(ff[OrdCTime==tt1[p]]) *exp(-sum(alpha[1:(kk-1)])) )   -  
#     sum( exp(ordX %*% beta) [ OrdTime>= tt1[kk-1] ] *   ordF[ OrdTime>= tt1[kk-1] ]    )

    for(p in 1:(kk-2)){        
      f[p] = sum(  1/(ff[OrdCTime ==max(tt1)  ] ) * (-exp(-sum(alpha)))  )      +
     sum(  1/(ff[OrdCTime <max(tt1) & OrdCTime>tt1[p]]) * ordf1 [(p+1):(kk-1)] )      +
      sum(  1/(ff[OrdCTime==tt1[p]]) *exp(-sum(alpha[1:p])) ) -  
      sum( exp(ordX %*% beta) [ OrdTime>= tt1[p] ] *   ordF[ OrdTime>= tt1[p] ]    )
    }
         
      f[kk-1] = sum(  1/(ff[OrdCTime ==max(tt1)  ] ) * (-exp(-sum(alpha))) )   +   
      sum(  1/(ff[OrdCTime==tt1[kk-1]]) *exp(-sum(alpha[1:(kk-1)])) )   -  
     sum( exp(ordX %*% beta) [ OrdTime>= tt1[kk-1] ] *   ordF[ OrdTime>= tt1[kk-1] ]    )
    
    return(f)
  }



  C1<-matrix(0,kk-1,kk-1)
  diag(C1)<-1
  D <- rep(0,kk-1)
#print("##########")
#print( Eqtion(alpha_esti))
#print( Loglike(alpha_esti))
  re_esti<-maxLik::maxLik(Loglike, start =alpha_esti, method = 'BFGS', grad = Eqtion )$estimate #
  
  alpha_esti<-as.double(re_esti)
  print(re_esti)
  gSS3 <- F(alpha_esti)$FF
  
  list(gSS3 = gSS3,alpha_esti=alpha_esti)
}