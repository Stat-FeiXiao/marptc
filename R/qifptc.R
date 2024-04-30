#' Title  Fit the marginal semiparametric promotion time cure model using QIF method
#'
#' @description Fit marginal promotion time cure (PTC) model based on the general estimation equations improved by quadratic inference functions (QIF) method. We consider three common correlation structures in this funciton.
#' The baseline cumulative distribution function in PTC model is appromximated by the Bernstein polynomial.
#'
#' @param formula a formula expression, of the form \code{response ~ predictors}. The \code{response} is a \code{Surv} object with right censoring. The expression to the right of the "~" specifies the effect of covariates on the failure time.
#' @param data  a data frame in which to interpret the variables named in the \code{formula}.
#' @param id a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.
#' @param Var if it is TRUE, the program returns Std.Error. By default, \code{Var = TRUE}.
#' @param Ibeta the initial value of the covariate coefficients.
#' @param boots if it is TRUE, the program returns Std.Error by the bootstrap method. By default, \code{boots = FALSE}
#' @param stad if it is TRUE, all the covariates in the \code{formula} are standardized. By default, \code{stdz = FALSE}.
#' @param nboot the number of bootstrap samples. The default is \code{nboot = 100}.
#' @param N the degree of Bernstein polynomial used to approximate the baseline distribution function in model.
#' @param tau the cure threshold. Individuals with survival time greater than \code{tau} are considered cured.
#' @param corstr a character string specifying the correlation structure. The following are permitted: \code{independence}, \code{exchangeable} and \code{AR1}.
#' @param itermax specifies the maximum iteration number. If the convergence criterion is not met, the iteration will be stopped after \code{esmax} iterations and the estimates will be based on the last iteration. The default \code{itermax = 100}.
#' @param eps tolerance for convergence. The default is \code{eps = 1e-6}. Iteration stops once the relative change in deviance is less than \code{eps}.
#'
#' @return An object of class \code{qifptc} is returned. It can be examined by \code{print.qifptc()}.
#'
#' @export qifptc
qifptc <- function(formula, data, id, Var = TRUE,Ibeta=NULL, stad=FALSE,boots=FALSE,nboot=100,N,tau=NULL,  corstr="independence",itermax = 100, eps = 1e-06) {
  call <- match.call()
  data <- data
 if (is.null(id)) {
        stop("Id variable not found")
    }
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
  if (QR$rank < ncol(X))
    stop("rank-deficient model matrix")
  beta_name <- colnames(X)
  beta_length <- ncol(X)
  Y <- model.extract(mf, "response")
  if (!inherits(Y, "Surv"))
    stop("Response must be a survival object")
   if (length(id) != dim(Y)[1])
        stop("Id and response not same length")
  CORSTRS <- c("independence", "exchangeable", "AR1")
  corstrv <- pmatch(corstr, CORSTRS, -1)
  if (corstrv == -1) stop("invalid corstr.")

  if(is.null(Ibeta)) Ibeta <- rep(0,dim(X)[2])
  Time <- Y[, 1]
  Status <- Y[, 2]
  if(is.null(tau)) tau <- max(Time[Status==1]) + 1e-4

  betafit <- qbeta(Time, Status, X,stad, Ibeta,id,itermax, eps,N,tau,corstr)
  beta1 <- betafit$beta1
  beta2 <- betafit$beta2
  Fhat <- betafit$Fhat
  gam <- betafit$gam

  fit <- list()
  fit$beta <-beta2

  if (Var) {

stadX <- X
    if(stad){
      for (i in 2:ncol(X)) {
        stadX[, i] <- (X[, i] - mean(X[, i]))/sd(X[, i])
      }
    }

      varfit <- varrqif(Time, Status,stadX,id, beta1, Fhat,tau,gam,N,corstr)
      fit$beta_var <- diag(varfit$var_beta)
      beta_var1 <- diag(varfit$var_beta)
      if(stad){
      AA1 <- c(1, -apply( X[, -1, drop = FALSE], 2, mean)/apply( X[, -1, drop = FALSE], 2, sd) )
      fit$beta_var <- c( (t(AA1))%*% (varfit$var_beta) %*% t(t(AA1))  ,(beta_var1[-1]) /((apply( X[,  -1, drop = FALSE], 2, sd)        )^2)       )             
    }
    fit$beta_sd <-  sqrt(fit$beta_var)
    fit$beta_zvalue <- fit$beta/fit$beta_sd
    fit$beta_pvalue <- (1 - stats::pnorm(abs(fit$beta_zvalue))) * 2
    fit$var_F <- varfit$var_F

  }


  if (boots) {
    Bootsample <- nboot
    BMs <- matrix(0, Bootsample, ncol(X))
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
      mf_boot <- model.frame(formula, bootdata)
      X_boot <- model.matrix(attr(mf_boot, "terms"), mf)
      X_boot <- as.matrix(X_boot)
      colnames(X_boot) <- colnames(model.matrix(attr(mf_boot, "terms"), mf_boot))
      Y_boot <- model.extract(mf_boot, "response")
      if (!inherits(Y_boot, "Surv"))
        stop("Response must be a survival object")
      Time_boot <- Y_boot[, 1]
      Status_boot <- Y_boot[, 2]
      tryboot <- try(
        withTimeout(
          qbeta(Time_boot, Status_boot, X_boot,stad,Ibeta,id = id_boot,itermax, eps,N,tau,corstr)
          ,  timeout=50 , onTimeout="error")
        , silent = F)
      if(is(tryboot,"try-error") == FALSE)    break
      }
      esfitboot <- tryboot
      BMs[rrn, ] <- esfitboot$beta2
    }
    var_beta_boots <- apply(BMs, 2, var)
    sd_beta_boots <- sqrt(var_beta_boots)

    fit$boots_beta_var <- var_beta_boots
    fit$boots_beta_sd <- sd_beta_boots
    fit$boots_beta_zvalue <- beta2/sd_beta_boots
    fit$boots_beta_pvalue <- (1 - pnorm(abs(fit$boots_beta_zvalue))) * 2
}

  fit$num_of_clusters <- K
  fit$max_cluster_size <- max(n)

  fit$call <- call

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
  cat("\nEstimation of Promotion Time Cure Model Based on QIF:\n")
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
  cat("Number of clusters:", x$num_of_clusters)
  cat("        Maximum cluster size:", x$max_cluster_size)
  cat("\n")
  invisible(x)
}
