
baseF<- function(Time, Status, X, beta)  {
  Kn <- length(Time)
  t2 <- Time
  t11 <- sort(Time)
  c11 <- Status[order(Time)]
  x111 <- as.matrix(X[order(Time), ])
  tt1 <- unique(t11[c11 == 1])
  kk <- length(table(t11[c11 == 1]))
  dd <- as.matrix(table(t11[c11 == 1]))
  gSS <- rep(0, kk)
  gSS1 <- rep(1, kk)
  g11 <- rep(1,length(Time))
  gSS[1] <- dd[1]/(sum(g11[min((1:Kn)[t11 == tt1[1]]):Kn] * exp(x111[min((1:Kn)[t11 == tt1[1]]):Kn, ]%*%beta )))
  for (i in 1:(kk - 1)) {
    gSS[i + 1] <- gSS[i] + dd[i + 1]/(sum(g11[min((1:Kn)[t11 == tt1[i + 1]]):Kn] * exp(x111[min((1:Kn)[t11 == tt1[i + 1]]):Kn, ]%*%beta )))
  }
  gSS1=exp(-gSS)
  gS=c(gSS[1],gSS[2:kk]-gSS[1:(kk-1)])

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
        gSS3[i] <-gSS[kk]
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
  gSS1 <- gSS1
  gSS3 <- gSS3
  gS <- gS


  list(gSS1=gSS1, gSS3 = gSS3, gS = gS,gs=gs)
}
