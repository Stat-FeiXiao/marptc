#' Title Fit the semiparametric proportional hazards cure model using GEE method improved by quadratic inference functions
#'
#' @param formula a formula expression, of the form \code{response ~ predictors}. The \code{response} is a \code{Surv} object with right censoring. The expression to the right of the "~" specifies the effect of covariates on the failure time.
#' @param data  a data frame in which to interpret the variables named in the \code{formula} and the \code{cureform}.
#' @param id a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.
#' @param Var If it is TRUE, the program returns Std.Error. By default, \code{Var = TRUE}.
#' @param esmax specifies the maximum iteration number. If the convergence criterion is not met, the iteration will be stopped after \code{esmax} iterations and the estimates will be based on the last iteration. The default \code{esmax = 100}.
#' @param eps tolerance for convergence. The default is \code{eps = 1e-6}. Iteration stops once the relative change in deviance is less than \code{eps}.
#'
#' @return An object of class \code{qifcure} is returned. It can be examined by \code{print.marcure()}.
#' @export
#'
#' @examples {
#' #Example. Fit the marginal semiparametric PHC model for the bmt data.
#'  data(bmt)
#'   qifbmt <- qifcure(Surv(T2, d3) ~ Z3,  data = bmt, id = bmt$Z9)
#'
#'}

qifcure <- function(formula, data, id, Var = TRUE, stad=FALSE,  esmax = 100, eps = 1e-06) {
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
  colnames(X) <- colnames(stats::model.matrix(attr(mf, "terms"), mf))[-1]
  beta_name <- colnames(X)
  beta_length <- ncol(X)
  Y <- model.extract(mf, "response")
  if (!inherits(Y, "Surv"))
    stop("Response must be a survival object")

  Time <- Y[, 1]
  Status <- Y[, 2]
  esmax <- esmax
  eps <- eps
  esfit <- esqif(Time, Status, X,stad, id,esmax, eps,formula, data)
  beta <- esfit$beta
  beta2 <- esfit$beta2
  hatF <- esfit$gSS3
  bgS <- esfit$gS
  fit <- list()
  class(fit) <- c("qifcure")

  fit$beta <- beta
  fit$beta2 <- beta2
  if (Var) {
    if(stad){
      stadX <-  std(X) #std(X)
      varfit <- varrqif(Time, Status,stadX,id, beta,gS=bgS,Lambda=hatF)
      fit$beta_var <- (varfit$var_beta)/ (attr(stadX ,"scale") ^2)
      stadX <-  scale(X) #std(X)
      gS<- baseF(Time, Status, stadX ,  beta)$gS
      Lambda<- baseF(Time, Status, stadX ,  beta)$ gSS3
      var_F <- vartso(Time, Status, stadX , id, beta, gS, Lambda)$var_F
      fit$gamma=gamma_0_with_var(gS,var_F)$gamma_0 
      fit$var_gamma=gamma_0_with_var(gS,var_F)$var_gamma
      fit$gamma_sd <- sqrt(fit$var_gamma)
      fit$gamma_zvalue <- fit$gamma/fit$gamma_sd
      fit$gamma_pvalue <- (1 - stats::pnorm(abs(fit$gamma_zvalue))) * 2
      fit$beta_sd <-  sqrt(fit$beta_var)
      fit$beta_zvalue <- beta2/fit$beta_sd
      fit$beta_pvalue <- (1 - stats::pnorm(abs(fit$beta_zvalue))) * 2
    }else{stadX <- X
  
    varfit <- varrqif(Time, Status,stadX,id, beta,gS=bgS,Lambda=hatF)
    fit$beta_var <- varfit$var_beta
    gS<- baseF(Time, Status, stadX ,  beta)$gS
    Lambda<- baseF(Time, Status, stadX ,  beta)$ gSS3
    var_F <- vartso(Time, Status, stadX , id, beta, gS, Lambda)$var_F
    fit$gamma=gamma_0_with_var(gS,var_F)$gamma_0 
    fit$var_gamma=gamma_0_with_var(gS,var_F)$var_gamma
    fit$gamma_sd <- sqrt(fit$var_gamma)
    fit$gamma_zvalue <- fit$gamma/fit$gamma_sd
    fit$gamma_pvalue <- (1 - stats::pnorm(abs(fit$gamma_zvalue))) * 2
    fit$beta_sd <-  sqrt(fit$beta_var)
    fit$beta_zvalue <- beta2/fit$beta_sd
    fit$beta_pvalue <- (1 - stats::pnorm(abs(fit$beta_zvalue))) * 2
    }  
 
      
   
  
  }


  fit$num_of_clusters <- K
  fit$max_cluster_size <- max(n)

  fit$call <- call

  fit$beta_name <- beta_name
  fit$Time <- Time
  fit$Var <- Var

  class(fit) = "qifcure"
  fit
}

print.qifcure <- function(x, ...) {
  if (!is.null(cl <- x$call)) {
    cat("Call:\n")
    dput(cl)
  }
  cat("\nPromotion time Cure Model:\n")
  if (x$Var) {
    bt <- array(c(x$gamma,x$beta2), c(1+length(x$beta2), 4))
    rownames(bt) <- c("Intercept",x$beta_name)
    colnames(bt) <- c("Estimate", "Std.Error", "Z value", "Pr(>|Z|)")
    bt[, 2] <- c(x$gamma_sd ,x$beta_sd)
    bt[, 3] <- c(x$gamma_zvalue ,x$beta_zvalue)
    bt[, 4] <- c(x$gamma_pvalue ,x$beta_pvalue)
  }
  else {
    bt <- array(x$beta2, c(length(x$beta), 1))
    rownames(bt) <- x$beta_name
    colnames(bt) <- "Estimate"
  }
  print(bt)
  cat("\n")
  cat("Number of clusters:", x$num_of_clusters)
  cat("       Maximum cluster size:", x$max_cluster_size)
  cat("\n")
  invisible(x)
}
