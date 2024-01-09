#' Title Estimation of the baseline cumulative distribution function
#'
#' @param Time right censored data which is the follow up time.
#' @param Status the censoring indicator, 1 = event of interest happens, and 0 = censoring.
#' @param X a matrix of covariates that that may have effect on the failure times.
#' @param beta initial values of covariates coefficients.
#'
#' @return estimation of the baseline cumulative distribution function
#' @export
#'
baseF<- function(Time, Status, X, beta,alpha_esti) {
  Kn <- length(Time)
  t2 <- Time
  t11 <- sort(Time)
  c11 <- Status[order(Time)]
  X <- as.matrix(X[,-1])
  beta <- as.double(beta)[-1]
  x111 <- as.matrix(X[order(Time), ])
  tt1 <- unique(t11[c11 == 1])
  kk <- length(table(t11[c11 == 1]))
  dd <- as.matrix(table(t11[c11 == 1]))
  gSS <- rep(0, kk)
  gSS1 <- rep(1, kk)
  g11 <- rep(1,length(Time))
  beta <- as.matrix(beta)
  gSS[1] <- dd[1]/(sum(exp(x111[min((1:Kn)[t11 == tt1[1]]):Kn, ]%*%beta )))
  for (i in 1:(kk - 1)) {
    gSS[i + 1] <- gSS[i] + dd[i + 1]/(sum(exp(x111[min((1:Kn)[t11 == tt1[i + 1]]):Kn, ]%*%beta )))
  }
  gSS1=exp(-gSS)
  gS=c(gSS[1],gSS[2:kk]-gSS[1:(kk-1)])
  gSS=gSS/max(gSS)
  gss=seq(1,kk)/kk
  gs <- rep(0, Kn)
  gSS3 <- rep(0, Kn)
  for (i in 1:(Kn)) {
    kk1 <- 1

    if (t2[i] < tt1[1]) {
      gs[i]<- 1e-08
      gSS3[i] <- 1e-08
    } else {
      if (t2[i] >= tt1[kk]) {
        gs[i]<- 1
        gSS3[i] <-  1
      } else {
        repeat {
          if(t2[i]>=tt1[kk1]) kk1=kk1+1
          else break
        }
        {
          gs[i] <- gss[kk1 - 1]
          gSS3[i] <- gSS[kk1 - 1]
        }
      }
    }
  }
  lambda=1        #proper c.d.f 是参数为2的指数分布 
  cureThres <- tau <- max(Time[Status==1])
F_0 = function(t){

 # f =   1 - exp(-lambda*t)
   f=(1-exp(- lambda*t)) / (1-exp(- lambda* tau)) * ifelse(0<=t&t<=tau,1,0)+1* ifelse(0<=t&t<=tau,0,1)
  return(f)
}

gSS3 <- F_0(Time)
  gSS1 <- gSS1
  gSS3 <- gSS3
  gS <- gS


  list(gSS1=gSS1, gSS3 = gSS3, gS = gS,gs=gs)
}
