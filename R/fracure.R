#' Title Semiparametric proportional hazards cure model Fit the semiparametric proportional hazards cure model using frailty method
#'
#' @param formula a formula expression, of the form \code{response ~ predictors}. The \code{response} is a \code{Surv} object with right censoring. The expression to the right of the "~" specifies the effect of covariates on the failure time.
#' @param data  a data frame in which to interpret the variables named in the \code{formula} and the \code{cureform}.
#' @param id a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.
#' @param Var If it is TRUE, the program returns Std.Error. By default, \code{Var = TRUE}.
#' @param emmax specifies the maximum iteration number. If the convergence criterion is not met, the iteration will be stopped after \code{esmax} iterations and the estimates will be based on the last iteration. The default \code{esmax = 100}.
#' @param eps tolerance for convergence. The default is \code{eps = 1e-6}. Iteration stops once the relative change in deviance is less than \code{eps}.
#'
#' @return An object of class \code{fracure} is returned. It can be examined by \code{print.marcure()}.
#' @export
#'
#' @examples {
#' #Example. Fit the marginal semiparametric PHC model for the bmt data.
#'  data(bmt)
#'   frabmt <- fracure(Surv(T2, d3) ~ Z3,  data = bmt, id = bmt$Z9)
#'
#'}

fracure <- function(formula, data, id, Var = TRUE,emmax = 100, eps = 1e-06) {
  call <- match.call()
  data <- data
  id <- id #这里输入的id可能不是按照顺序来的，而且id可能不是连续的即1，3，4，6，，，，
  uid <- sort(unique(id))
  newid <- rep(0, length(id))
  for (i in 1:length(id)) {
    j <- 1
    repeat {
      if (id[i] != uid[j])
        j <- j + 1 else {
          newid[i] <- j           #newid 是将id转换成连续的数，即newid是没有缺某个数的id列
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
  X <- stats::model.matrix(attr(mf, "terms"), mf)[,-1]   #[, -1]去掉截距项
  X <- as.matrix(X)
  colnames(X) <- colnames(stats::model.matrix(attr(mf, "terms"), mf))[-1]
  beta_name <- colnames(X)
  beta_length <- ncol(X)
  Y <- model.extract(mf, "response")
  if (!inherits(Y, "Surv"))
    stop("Response must be a survival object")

  Time <- Y[, 1]
  Status <- Y[, 2]
  emmax <- emmax
  eps <- eps
  emfit <- em(Time, Status, X,id,emmax , eps)
  beta <- emfit$beta
  stat <- emfit$stat
  gSS3 <- emfit$gSS3
  Lambda <- emfit$Lambda
  varrest<- emfit$varrest
  bgS <- emfit$gS
  fit <- list()
  class(fit) <- c("fracure")

  fit$beta <- beta
  if (Var) {
    varfit <- varrfra(Time, Status, X,id, beta,gS=bgS,Lambda, gSS3,stat,varrest)
    fit$beta_var <- varfit$var_beta
    fit$beta_sd <- varfit$sd_beta
    fit$beta_zvalue <- beta/fit$beta_sd
    fit$beta_pvalue <- (1 - stats::pnorm(abs(fit$beta_zvalue))) * 2
  }


  fit$num_of_clusters <- K
  fit$max_cluster_size <- max(n)

  fit$call <- call

  fit$beta_name <- beta_name
  fit$Time <- Time
  fit$Var <- Var

  class(fit) = "fracure"
  fit
}

print.fracure <- function(x, ...) {
  if (!is.null(cl <- x$call)) {
    cat("Call:\n")
    dput(cl)
  }
  cat("\nNon-Mixture Cure Model:\n")
  if (x$Var) {
    bt <- array(x$beta, c(length(x$beta), 4))
    rownames(bt) <- x$beta_name
    colnames(bt) <- c("Estimate", "Std.Error", "Z value", "Pr(>|Z|)")
    bt[, 2] <- x$beta_sd
    bt[, 3] <- x$beta_zvalue
    bt[, 4] <- x$beta_pvalue

  }
  else {
    bt <- array(x$beta, c(length(x$beta), 1))
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
