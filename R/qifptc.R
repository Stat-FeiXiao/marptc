#' Title  Fit the marginal semiparametric promotion time cure model under clustered survival time
#'
#' @description Fit the marginal promotion time cure (PTC) model based on the general estimation equations (GEE) and the quadratic inference functions (QIF).
#' We consider three common correlation structures (\code{independence}, \code{exchangeable} and \code{AR1}) in this funciton.
#' The variances can be consistently estimated by a sandwich variance estimator. We also present a bootstrap procedure for variance estimation.
#'
#' @param formula a formula expression, of the form \code{response ~ predictors}. The \code{response} is a \code{Surv} object with right censoring.
#' The expression to the right of the "~" specifies the effect of covariates on the failure time.
#' @param data  a data frame in which to interpret the variables named in the \code{formula}.
#' @param id a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.
#' @param Var if it is TRUE, the program returns Std.Error. By default, \code{Var = TRUE}.
#' @param Ibeta the initial value of the covariate coefficients.
#' @param boots if it is TRUE, the program returns Std.Error by the bootstrap method. By default, \code{boots = FALSE}
#' @param stad if it is TRUE, all the covariates in the \code{formula} are standardized. By default, \code{stdz = TRUE}.
#' @param nboot the number of bootstrap samples. The default is \code{nboot = 100}.
#' @param method using GEE or QIF method to estimate the covariate coefficients.
#' @param corstr a character string specifying the correlation structure. The following are permitted: \code{independence}, \code{exchangeable} and \code{AR1}.
#' @param itermax specifies the maximum iteration number. If the convergence criterion is not met, the iteration will be stopped after \code{esmax} iterations and
#' the estimates will be based on the last iteration. The default \code{itermax = 100}.
#' @param eps tolerance for convergence. The default is \code{eps = 1e-6}. Iteration stops once the relative change in deviance is less than \code{eps}.
#'
#' @return An object of class \code{qifptc} is returned. It can be examined by \code{print.qifptc()}.
#'
#' @export qifptc

qifptc <- function(formula, data, id, Var = TRUE,Ibeta=NULL, stad=TRUE,boots=FALSE,nboot=100, method = "GEE", corstr="independence",itermax = 100, eps = 1e-06) {

  call <- match.call()

 if (is.null(id))  stop("Id variable not found!")

  uid <- sort(unique(id))
  newid <- rep(0, length(id))
  for (i in 1:length(id)) {
    j <- 1
    repeat {
      if (id[i] != uid[j])
        j <- j + 1 else {
          newid[i] <- j
          break
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
  X <- stats::model.matrix(attr(mf, "terms"), mf)

  QR <- qr(X)
  if (QR$rank < ncol(X))  stop("rank-deficient model matrix!")

  beta_name <- colnames(X)

  Y <- model.extract(mf, "response")

  if (!inherits(Y, "Surv"))  stop("Response must be a survival object!")

  if (length(id) != nrow(Y))  stop("Id and response not same length!")

  Time <- Y[, 1]
  Status <- Y[, 2]

  CORSTRS <- c("independence", "exchangeable", "AR1")
  corstrv <- pmatch(corstr, CORSTRS, -1)
  if (corstrv == -1)  stop("invalid corstr!")

  Method <- c("GEE", "QIF")
  corstrv <- pmatch(method, Method, -1)
  if (corstrv == -1)  stop("invalid method!")

  if(is.null(Ibeta)){

    KMfit <- survival::survfit(survival::Surv(Time,Status)~1)
    cure.rate <- min(KMfit$surv)
    Ibeta <- c(ifelse(cure.rate==0,0,log(-log(cure.rate))),rep(0,ncol(X)-1))

  }else{
  if(length(Ibeta) != ncol(X)) {stop("Dimension of beta != ncol(X)!")}
  }

  fit <- list()

  betafit <- gqbeta(Time, Status, X,stad, Ibeta,id,itermax, eps,method,corstr)

  behat <- betafit$beta2

  if(method=="GEE")  fit$rho <- betafit$rho

  fit$beta <- behat

  if (Var) {

      if(method=="GEE"){
        varfit <- vargee(Time, Status, X, id, behat,corstr)
      }else if(method=="QIF"){
      varfit <- varqif(Time, Status, X, id, behat,corstr)
      }
    fit$beta_var <- varfit$var_beta
    fit$beta_sd <-  varfit$sd_beta
    fit$beta_zvalue <- fit$beta/fit$beta_sd
    fit$beta_pvalue <- (1 - stats::pnorm(abs(fit$beta_zvalue))) * 2

  }


  if (boots) {
    Bootsample <- nboot
    BMs <- matrix(0, Bootsample, ncol(X))
    BMnus <- rep(0, Bootsample)
    for (rrn in 1:Bootsample) {
      repeat{  bootid<- sample((1:K), replace = TRUE)
      bootdata <- data[id == bootid[1], ]
      bootdata$id <- rep(1, sum(id == bootid[1]))
      for (ll in 2:K) {
        bootdata1 <- data[id == bootid[ll], ]
        bootdata1$id <- rep(ll, sum(id == bootid[ll]))
        bootdata <- rbind(bootdata, bootdata1)
      }
      id_boot <- bootdata$id
      Kn_boot <- length(id_boot)
      K_boot <- length(unique(id_boot))
      n_boot <- as.vector(table(id_boot))
      mf <- model.frame(formula, bootdata)
      X <- model.matrix(attr(mf, "terms"), mf)
      X <- as.matrix(X)
      colnames(X) <- colnames(model.matrix(attr(mf, "terms"), mf))
      Y <- model.extract(mf, "response")
      if (!inherits(Y, "Surv"))
        stop("Response must be a survival object")
      Time <- Y[, 1]
      Status <- Y[, 2]
      tryboot <- try(
        withTimeout(
          gqbeta(Time, Status, X,stad, Ibeta,id = id_boot,itermax, eps,method,corstr)
          ,  timeout=300 , onTimeout="error")
        , silent = F)
      if(is(tryboot,"try-error") == FALSE)    break
      }
      esfitboot <- tryboot
      BMs[rrn, ] <- esfitboot$beta2
      if(method=="GEE") BMnus[rrn] <- esfitboot$rho
    }

    fit$boots_beta_var <- apply(BMs, 2, var)
    fit$boots_beta_sd <- apply(BMs, 2, sd)
    fit$boots_beta_zvalue <- fit$beta/fit$boots_beta_sd
    fit$boots_beta_pvalue <- (1 - pnorm(abs(fit$boots_beta_zvalue))) * 2
    if(method=="GEE"){
      fit$boots_rho_var<- var(BMnus)
      fit$boots_rho_sd <- sd(BMnus)
      fit$boots_rho_zvalue <- fit$rho/fit$sd_rho_boots
      fit$boots_rho_pvalue <- (1 - pnorm(abs(fit$boots_rho_zvalue))) * 2
    }
}

  fit$num_of_clusters <- K
  fit$max_cluster_size <- max(n)
  fit$call <- call
  fit$method <- method
  fit$beta_name <- beta_name
  fit$Time <- Time
  fit$Var <- Var
  fit$boots <- boots
  class(fit) = "qifptc"
  fit
}

#' Title Print qifptc object
#'
#' @param x an object of \code{qifptc}.
#' @param ... further arguments to be added in the \code{print.qifptc} function.
#'
#' @export print.qifptc

print.qifptc <- function(x, ...) {
  if (!is.null(cl <- x$call)) {
    cat("Call:\n")
    dput(cl)
  }
  cat("\nEstimation of Marginal Promotion Time Cure Model:\n")
  if (x$Var) {
    bt <- array(c(x$beta), c(length(x$beta), 4))
    rownames(bt) <- c(x$beta_name)
    colnames(bt) <- c("Estimate", "Std.Error", "Z value", "Pr(>|Z|)")
    if(x$boots){
    bt[, 2] <- c(x$boots_beta_sd)
    bt[, 3] <- c(x$boots_beta_zvalue)
    bt[, 4] <- c(x$boots_beta_pvalue)
  }else{
    bt[, 2] <- c(x$beta_sd)
    bt[, 3] <- c(x$beta_zvalue)
    bt[, 4] <- c(x$beta_pvalue)
 }

  }
  else {
    bt <- array(c(x$beta), c(length(x$beta), 1))
    rownames(bt) <-  c(x$beta_name)
    colnames(bt) <- "Estimate"
  }
  print(bt)
  cat("\n")
  cat("Estimation method:", x$method)
  cat("\n")
  cat("Number of clusters:", x$num_of_clusters)
  cat("        Maximum cluster size:", x$max_cluster_size)
  cat("\n")
  invisible(x)
}
