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
#' @param stad if it is TRUE, all the covariates in the \code{formula} are standardized. By default, \code{stad = TRUE}.
#' @param nboot the number of bootstrap samples. The default is \code{nboot = 100}.
#' @param method using GEE or QIF method to estimate the covariate coefficients.
#' @param corstr a character string specifying the correlation structure. The following are permitted: \code{independence}, \code{exchangeable} and \code{AR1}.
#' @param IC if \code{TRUE}, the function returns the QIC for the GEE method or the BIQIF for the QIF method. Default is \code{FALSE}.
#' @param itermax specifies the maximum iteration number. If the convergence criterion is not met, the iteration will be stopped after \code{itermax} iterations and
#' the estimates will be based on the last iteration. The default \code{itermax = 100}.
#' @param eps tolerance for convergence. The default is \code{eps = 1e-6}. Iteration stops once the relative change in deviance is less than \code{eps}.
#' @return An object of class \code{marptc} is returned. It can be examined by \code{print.marptc()}.
#'
#' @export marptc

marptc <- function(formula, data, id, Var = TRUE, Ibeta = NULL, stad = TRUE, boots = FALSE, 
                   nboot = 100, method = "GEE", corstr = "independence", IC = FALSE, itermax = 100, eps = 1e-06) {

  call <- match.call()

  if (is.null(id)) stop("Id variable not found!")

  data$id <- match(id, sort(unique(id)))
  data <- data[order(data$id), , drop = FALSE]
  id <- data$id
  
  K  <- length(unique(id))
  n  <- as.vector(table(id))   

  mf <- stats::model.frame(formula, data)
  X  <- stats::model.matrix(attr(mf, "terms"), mf)

  QR <- qr(X)
  if (QR$rank < ncol(X)) stop("Rank-deficient model matrix!")

  beta_name <- colnames(X)
  Y <- stats::model.extract(mf, "response")

  if (!inherits(Y, "Surv")) stop("Response must be a survival object!")
  if (length(id) != nrow(Y)) stop("Id and response not same length!")

  Time <- Y[, 1]
  Status <- Y[, 2]

  CORSTRS <- c("independence", "exchangeable", "AR1")
  if (!(corstr %in% CORSTRS)) stop("Invalid corstr! Choose 'independence', 'exchangeable', or 'AR1'.")
  if (!(method %in% c("GEE", "QIF"))) stop("Invalid method! Choose 'GEE' or 'QIF'.")

  if (is.null(Ibeta)) {
    KMfit <- survival::survfit(survival::Surv(Time, Status) ~ 1)
    cure.rate <- min(KMfit$surv)
    Ibeta <- c(ifelse(cure.rate == 0, 0, log(-log(cure.rate))), rep(0, ncol(X) - 1))
  } else {
    if (length(Ibeta) != ncol(X)) stop("Dimension of Ibeta does not match ncol(X)!")
  }

  fit <- list()

  betafit <- gqbeta(Time, Status, X, stad, Ibeta, id, IC, itermax, eps, method, corstr)
  
  behat <- betafit$beta2

  if(method=="GEE"){
    fit$rho <- betafit$rho
    fit$pphi <- betafit$pphi
    fit$V_Ind <- betafit$V_Ind
    fit$Q_beta <- betafit$Q_beta
  }else{
    fit$Quad <- betafit$Quad
  }

  fit$beta <- behat
  fit$convergence <- betafit$convergence
  fit$baseF <- betafit$baseF

  if (Var) {
    if (method == "GEE") {
      varfit <- vargee(Time, Status, X, id, behat, corstr)
    } else if (method == "QIF") {
      varfit <- varqif(Time, Status, X, id, behat, corstr)
    }
    fit$beta_var <- varfit$var_beta
    fit$beta_sd <-  varfit$sd_beta
    fit$beta_var_matrix <- varfit$var_beta_matrix
    fit$beta_zvalue <- fit$beta/fit$beta_sd
    fit$beta_pvalue <- (1 - stats::pnorm(abs(fit$beta_zvalue))) * 2
  }

  if(IC){
  if(method=="GEE"){
    fit$QIC <- -2*fit$Q_beta + 2*sum(diag(  as.matrix(fit$V_Ind%*%(fit$beta_var_matrix)) ) ) 
  }else if(method=="QIF"){
    fit$BIQIF <- as.numeric(fit$Quad) + ncol(X) * log(K)
  }
}

if (boots) {
    cluster_idx_list <- split(seq_len(nrow(X)), id)
    pbet <- ncol(X)
    
    n_cores <- parallel::detectCores()

    cores_to_use <- max(1, ifelse(is.na(n_cores), 1, n_cores - 3)) 
	cl <- parallel::makeCluster(cores_to_use)
    doParallel::registerDoParallel(cl)
    
    # Export internal solver explicitly to avoid environment scoping issues
    boot_res <- foreach::foreach(rrn = 1:nboot, .combine = "rbind", 
                                 .packages = c("stats", "survival", "R.utils", "doRNG"),
                                 .export = c("gqbeta","geq","qeq","getlambda")) %dorng% {
          

      current_beta <- rep(NA, pbet)
      current_rho <- NA
      
           repeat{
  
        bootid <- sample(K, replace = TRUE)
        boot_rows <- unlist(cluster_idx_list[bootid], use.names = FALSE)
        
        cluster_sizes <- sapply(cluster_idx_list[bootid], length)
        id_boot <- rep(seq_len(K), times = cluster_sizes)
        
        Time_boot <- Time[boot_rows]
        Status_boot <- Status[boot_rows]
        X_boot <- X[boot_rows, , drop = FALSE]
        
        tryboot <- tryCatch({
          R.utils::withTimeout({
      tryboot <-      gqbeta(Time = Time_boot, Status = Status_boot, X = X_boot, 
                   stad = stad, Ibeta = behat, id = id_boot, IC = IC, 
                   itermax = itermax, eps = eps, method = method, corstr = corstr)
          }, timeout = 300, onTimeout = "error")
        }, error = function(e) {
           list(beta2 = rep(NA, pbet))
        })
  
       if (!any(is.na(tryboot$beta2))) {
          current_beta <- as.numeric( tryboot$beta2)
          if (method == "GEE") current_rho <-  as.numeric(tryboot$rho)
        break  
       }}
      return(c(current_beta, current_rho))
    }

    parallel::stopCluster(cl) 
    foreach::registerDoSEQ() 
    gc()
    
 
    BMs <- boot_res[, 1:pbet, drop = FALSE]
    fit$boots_beta_var <- apply(BMs, 2, stats::var)
    fit$boots_beta_sd <- apply(BMs, 2, stats::sd)
    fit$boots_beta_zvalue <- fit$beta / fit$boots_beta_sd
    fit$boots_beta_pvalue <- (1 - stats::pnorm(abs(fit$boots_beta_zvalue))) * 2
    
    if (method == "GEE") {
       BMnus <- boot_res[, pbet + 1]
        fit$boots_rho_var <- stats::var(BMnus)
        fit$boots_rho_sd <- stats::sd(BMnus)
        fit$boots_rho_zvalue <- fit$rho / fit$boots_rho_sd
        fit$boots_rho_pvalue <- (1 - stats::pnorm(abs(fit$boots_rho_zvalue))) * 2
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
  class(fit) <- "marptc"
  return(fit)
}


#' Title Print marptc object
#'
#' @param x an object of \code{marptc}.
#' @param ... further arguments to be added in the \code{print.marptc} function.
#'
#' @export print.marptc

print.marptc <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  if (!is.null(cl <- x$call)) {
    cat("\nCall:\n")
    dput(cl)
    cat("\n")
  }
  
  cat(sprintf("Marginal Promotion Time Cure Model (Method: %s)\n", x$method))
  cat(sprintf("Number of clusters: %d | Maximum cluster size: %d\n", x$num_of_clusters, x$max_cluster_size))
  cat(rep("-", 50), "\n", sep = "")
  
  cat("\nCoefficients:\n")
  if (x$Var) {
    bt <- matrix(0, nrow = length(x$beta), ncol = 4)
    rownames(bt) <- x$beta_name
    colnames(bt) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    
    bt[, 1] <- x$beta
      bt[, 2] <- x$beta_sd
      bt[, 3] <- x$beta_zvalue
      bt[, 4] <- x$beta_pvalue
    stats::printCoefmat(bt, digits = digits, P.values = TRUE, has.Pvalue = TRUE)
  } else {
    bt <- matrix(x$beta, ncol = 1)
    rownames(bt) <- x$beta_name
    colnames(bt) <- "Estimate"
    print(bt, digits = digits)
  }
  
  if (x$method == "GEE" && !is.null(x$rho)) {
    cat("\nCorrelation Parameter (rho):\n")
    if (x$boots) {
      hatr <- matrix(c(x$rho, x$boots_rho_sd, x$boots_rho_zvalue, x$boots_rho_pvalue), nrow = 1)
      rownames(hatr) <- "rho"
      colnames(hatr) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
      stats::printCoefmat(hatr, digits = digits, P.values = TRUE, has.Pvalue = TRUE)
    } else {
      cat(sprintf("Estimate: %s\n", format(x$rho, digits = digits)))
    }
  }
  
  cat("\n")
  invisible(x)
}