#' Title Estimation of the initial covariates coefficients using Tsodikov method
#'
#' @param Time right censored data which is the follow up time.
#' @param Status the censoring indicator, 1 = event of interest happens, and 0 = censoring.
#' @param X a matrix of covariates that that may have effect on the failure times.
#' @param id a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.
#'
#' @return  the initial covariate coefficients.
#' @export
tbeta <- function(Time, Status, X,  id,stad) {

  c1 <- Status
  t2 <- Time
  K <- length(unique(id))
  n <- as.vector(table(id))
  Kn <- sum(n)
  xxx<-as.matrix(X)
  if(stad){
    stadX <-  scale(xxx) #std(X)
  }else{stadX <- xxx}
  xxx <- stadX
  t11 <- sort(t2)
  c11 <- Status[order(t2)]
  x111 <- as.matrix(xxx[order(t2), ])
  tt1 <- unique(t11[c11 == 1])
  kk <- length(table(t11[c11 == 1]))
  dd <- as.matrix(table(t11[c11 == 1]))
  betaest<-rep(0,dim(xxx)[2])
  fun1=function(B1){
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
      too <- t2[c1==1]
      a  <- 0
      for(o in 1:length(t2[c1==1])){
        b<- (zegasiBx/zegasiB)[too[o]==tt1]
        a<-a+b
      }

      #  fun1=-zegasiBx/zegasiB
      fun2=c11*x111[,s]
      fi1[s]<-sum(fun2)-a
    }
    return(fi1)
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

  repeat{
    oldbetaest= betaest
    betaest=betaest-MASS::ginv(fun11(betaest))%*%fun1(betaest)
    if( all(abs( betaest- oldbetaest)<1e-6) ) break
  }


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
  medvar<-matrix(0,dim(xxx)[2],dim(xxx)[2])
  for(ii in 1:K){
    medvar <- medvar + Ui(ii,betaest)%*%t(Ui(ii,betaest))
  }
  varrt<-(MASS::ginv(fun11(betaest))%*%medvar%*%MASS::ginv(fun11(betaest)))

  varrt.nv<- (-MASS::ginv(fun11(betaest)))
  if(stad){
    betaest2 <- betaest/attr(stadX,"scaled:scale")
  }else{betaest2 <- betaest}
list(tbeta = betaest,tbeta2 = betaest2,varrt=varrt,varrt.nv=varrt.nv)
}


tsocure <- function(formula, data, id, Var = TRUE,stad=FALSE) {
  call <- match.call()
  data <- data
  id <- id #这里输入的id可能不是按照顺序来的，而且id可能不是连续的即1???3???4???6，，，，
  uid <- sort(unique(id))
  newid <- rep(0, length(id))
  for (i in 1:length(id)) {
    j <- 1
    repeat {
      if (id[i] != uid[j])
        j <- j + 1 else {
          newid[i] <- j           #newid 是将id转换成连续的数，即newid是没有缺某个数的id???
          break                   #id=c( 3 ,1,5) => newid=c( 2,1,3)
        }
    }
  }
  data$id <- newid
  data1 <- data[data$id == 1, ]
  for (i in 2:length(uid)) {
    data1 <- rbind(data1, data[data$id == i, ])
  }
  data <- data1
  id <- data$id
  Kn <- length(id)
  K <- length(unique(id))
  n <- as.vector(table(id))
  mf <- stats::model.frame(formula, data)
  X <- stats::model.matrix(attr(mf, "terms"), mf)[,-1]   #[, -1]去掉截距???
  X <- as.matrix(X)
  colnames(X) <- colnames(model.matrix(attr(mf, "terms"), mf))[-1]
  beta_name <- colnames(X)
  beta_length <- ncol(X)
  Y <- model.extract(mf, "response")
  if (!inherits(Y, "Surv"))
    stop("Response must be a survival object")

  Time <- Y[, 1]
  Status <- Y[, 2]
  beta <- tbeta(Time, Status, X,id,stad)$tbeta
  beta2 <- tbeta(Time, Status, X,id,stad)$tbeta2
  fit <- list()
  class(fit) <- c("tsocure")

  fit$beta <- beta
  fit$beta2 <- beta2
  if(Var){
    varrt=tbeta(Time, Status, X,id,stad)$varrt
    varrt.nv=tbeta(Time, Status, X,id,stad)$varrt.nv
    if(stad){
    stadX <-  scale(X) #std(X)
    fit$beta_var <- diag(varrt) / (attr(stadX ,"scaled:scale") ^2)
    gS<- baseF(Time, Status, stadX ,  beta)$gS
    Lambda<- baseF(Time, Status, stadX ,  beta)$ gSS3
    var_F <- vartso(Time, Status, stadX , id, beta, gS, Lambda)$var_F
    fit$gamma=gamma_0_with_var(gS,var_F)$gamma_0 
    fit$var_gamma=gamma_0_with_var(gS,var_F)$var_gamma
    fit$gamma_sd <- sqrt(fit$var_gamma)
    fit$gamma_zvalue <- fit$gamma/fit$gamma_sd
    fit$gamma_pvalue <- (1 - stats::pnorm(abs(fit$gamma_zvalue))) * 2
    }else{fit$beta_var <- diag(varrt) 
    fit$beta_sd <- sqrt(fit$beta_var)
    fit$beta_zvalue <- beta2/fit$beta_sd
    fit$beta_pvalue <- (1 - stats::pnorm(abs(fit$beta_zvalue))) * 2
    fit$beta_var.nv <- diag(varrt.nv)
    stadX <- X #std(X)
    gS<- baseF(Time, Status, stadX ,  beta)$gS
    Lambda<- baseF(Time, Status, stadX ,  beta)$ gSS3
    var_F <- vartso(Time, Status, stadX , id, beta, gS, Lambda)$var_F
    fit$gamma=gamma_0_with_var(gS,var_F)$gamma_0 
    fit$var_gamma=gamma_0_with_var(gS,var_F)$var_gamma
    fit$gamma_sd <- sqrt(fit$var_gamma)
    fit$gamma_zvalue <- fit$gamma/fit$gamma_sd
    fit$gamma_pvalue <- (1 - stats::pnorm(abs(fit$gamma_zvalue))) * 2
    }
    
  }



  fit$num_of_clusters <- K
  fit$max_cluster_size <- max(n)

  fit$call <- call

  fit$beta_name <- beta_name
  fit$Time <- Time


  class(fit) = "tsocure"
  fit
}


print.tsocure <- function(x, ...) {
  if (!is.null(cl <- x$call)) {
    cat("Call:\n")
    dput(cl)
  }
  cat("\nPromotion time Cure Model:\n")
 

    bt <- array(c(x$gamma,x$beta2), c(1+length(x$beta2), 4))
    rownames(bt) <- c("Intercept",x$beta_name)
    colnames(bt) <- c("Estimate", "Std.Error", "Z value", "Pr(>|Z|)")
    bt[, 2] <- c(x$gamma_sd ,x$beta_sd)
    bt[, 3] <- c(x$gamma_zvalue ,x$beta_zvalue)
    bt[, 4] <- c(x$gamma_pvalue ,x$beta_pvalue)

  

  print(bt)
  cat("\n")
  cat("Number of clusters:", x$num_of_clusters)
  cat("       Maximum cluster size:", x$max_cluster_size)
  cat("\n")
  invisible(x)
}

