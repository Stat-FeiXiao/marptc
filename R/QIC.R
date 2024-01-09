qif <- function(Time,Status,Lambda,fhat,beta,N) {
 Kn <- length(Time)
    t2 <- Time
    t11 <- sort(Time)
    c11 <- Status[order(Time)]
    x111 <- as.matrix(X[order(Time), ])
    tt1 <- unique(t11[c11 == 1])

    kk <- length(table(t11[c11 == 1]))
    K <- length(unique(id))
    nt <- as.vector(table(id))
    dd <- as.matrix(table(t11[c11 == 1]))
    OM <- matrix(0,p,p)
    beta <- fitf $ beta
    Lambda<- baseF(Time, Status, X, beta)$gSS3
    fhat<- baseF(Time, Status,X, beta)$fhat
   
    xx11 <- as.matrix(x111[c11 == 1,])
    fun1<- Status log(exp(xx11%*%beta)*fhat)
    fun2<- -exp(X%*%beta)*Lambda
    loglikelihood<- sum(fun1)+sum(fun2)

   qic <- -2*loglikelihood+ 2*(dim(X)[2]+N+1)
   return(qic)
}