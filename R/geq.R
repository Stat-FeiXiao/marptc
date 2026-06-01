#' Title General estimation equations for the covariate coefficients estimation
#'
#' @description Estimation the covariate coefficients using the GEE approach.
#'
#' @param Time right censored data which is the follow up time.
#' @param Status right censored data which is the follow up time.
#' @param X a matrix of covariates that that may have effect on the failure times.
#' @param beta the initial value of the covariate coefficients.
#' @param id a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.
#' @param IC if \code{TRUE}, the function returns the QIC for the GEE method or the BIQIF for the QIF method. Default is \code{FALSE}.
#' @param eps tolerance for convergence. Iteration stops once the relative change in deviance is less than \code{eps}.
#' @param corstr a character string specifying the correlation structure. The following are permitted: \code{independence}, \code{exchangeable} and \code{AR1}.
#' @param itermax specifies the maximum iteration number. If the convergence criterion is not met, the iteration will be stopped after \code{itermax} iterations and
#' the estimates will be based on the last iteration.
#'
#' @export geq

geq <- function(Time, Status, X, beta, id, IC =FALSE, eps = 1e-6, corstr, itermax) {
  
  K <- length(unique(id))
  n <- as.vector(table(id))
  gbeta <- beta1 <- beta
  p <- ncol(X)
  
  ord <- order(Time)
  Time_sort <- Time[ord]
  first_idx <- match(Time, Time_sort)
  
  idx1 <- which(Status == 1)
  Time1 <- Time[idx1]
  ord1 <- order(Time1)
  Time1_sort <- Time1[ord1]
  idx_interval <- findInterval(Time, Time1_sort)
  
  cluster_idx <- split(seq_along(id), id)
  SK1 <- 1
  
  repeat {
    mu <- exp(as.vector(X %*% gbeta)) 
    
    mu_sort <- mu[ord]
    cum_mu_pad <- c(0, cumsum(mu_sort))
    dRbet <- sum(mu) - cum_mu_pad[first_idx]
    
    # 极速预切割
    dRbet_d <- dRbet[Status == 1]
    Rbet.min <- min(dRbet_d)
    interval <- c(Rbet.min - sum(Status), Rbet.min - 1)
    lambet <- uniroot(getlambda, interval, tol = .Machine$double.eps^0.75, dRbet = dRbet_d)$root
    
    term <- 1 / (dRbet_d - lambet)
    term_sort <- term[ord1]
    cum_term_pad <- c(0, cumsum(term_sort))
    baseF <- cum_term_pad[idx_interval + 1]
    
    baseF1 <- baseF
    baseF1[baseF1 == 0] <- 1e-08
    newY1 <- Status / baseF1
    
    res <- (newY1 - mu) / sqrt(mu)
    pphi <- sum(res^2) / (length(Time) - p)
    res_list <- split(res, id)
    
    if (corstr == "independence") {
      rho <- 0
    } else if (corstr == "exchangeable") {
      rres <- sum(vapply(res_list, function(v) {
        if (length(v) == 1) return(v)
        (sum(v)^2 - sum(v^2)) / 2
      }, numeric(1)))
      rho <- (pphi^(-1)) * rres / (sum(n * (n - 1)) / 2 - p)
      
      # =======================================================================
      # 【核心修改 1】：动态限制 exchangeable 结构的 rho 下界，防止矩阵分母出现 0 或负数
      # =======================================================================
      max_n <- max(n)
      min_rho <- ifelse(max_n > 1, -1 / (max_n - 1) + 1e-4, -0.95)
      rho <- max(min(rho, 0.95), min_rho) 
      # =======================================================================
      
    } else if (corstr == "AR1") {
      rres <- sum(vapply(res_list, function(v) {
        if (length(v) == 1) return(v)
        sum(v[-length(v)] * v[-1])
      }, numeric(1)))
      rho <- (pphi^(-1)) * rres / (sum(n - 1) - p)
      rho <- max(min(rho, 0.95), -0.95) # AR1 不受影响，保留固定边界
    }
    
    if (corstr == "exchangeable") {
      inv_phi <- 1 / pphi
      a_exch <- inv_phi / (1 - rho)
      b_exch_vec <- inv_phi * (-rho / ((1 - rho) * (1 - rho + n * rho)))
    } else if (corstr == "AR1") {
      invQC_list <- vector("list", K)
      c_ar1 <- (1 / pphi) / (1 - rho^2)
      for (i in seq_len(K)) {
        ni <- n[i]
        if (ni == 1) {
          invQC_list[[i]] <- matrix(c_ar1, 1, 1)
        } else {
          invQC <- matrix(0, ni, ni)
          diag(invQC) <- c_ar1 * (1 + rho^2)
          invQC[1, 1] <- c_ar1
          invQC[ni, ni] <- c_ar1
          idx_off <- abs(row(invQC) - col(invQC)) == 1
          invQC[idx_off] <- -c_ar1 * rho
          invQC_list[[i]] <- invQC
        }
      }
    }
    
    S1 <- Status - mu * baseF
    G <- matrix(0, p, 1)
    G1 <- matrix(0, p, p)
    
    for (i in seq_len(K)) {
      idx <- cluster_idx[[i]]
      X_i <- X[idx, , drop = FALSE]
      mu_i <- mu[idx]
      sqrt_mu_i <- sqrt(mu_i)
      baseF_i <- baseF[idx]
      
      A <- t(X_i * sqrt_mu_i)                
      B_vec <- S1[idx] / sqrt_mu_i           
      C_mat <- X_i * (baseF_i * sqrt_mu_i)   
      
      if (corstr == "independence") {
        inv_phi <- 1 / pphi
        G_i <- inv_phi * (A %*% B_vec)
        G1_i <- -inv_phi * (A %*% C_mat)
      } else if (corstr == "exchangeable") {
        term1_G  <- A %*% B_vec
        term2_G  <- rowSums(A) * sum(B_vec)
        G_i <- a_exch * term1_G + b_exch_vec[i] * term2_G
        
        term1_G1 <- A %*% C_mat
        term2_G1 <- outer(rowSums(A), colSums(C_mat)) 
        G1_i <- - (a_exch * term1_G1 + b_exch_vec[i] * term2_G1)
      } else if (corstr == "AR1") {
        part1 <- A %*% invQC_list[[i]]
        G_i <- part1 %*% B_vec
        G1_i <- - (part1 %*% C_mat)
      }
      
      G <- G + G_i
      G1 <- G1 + G1_i
    }
    
    gbeta_update <- tryCatch(solve(G1, G), error = function(e) MASS::ginv(G1) %*% G)
    gbeta <- gbeta - gbeta_update
    
    if (any(abs(gbeta - beta1) > eps) && SK1 <= itermax) {
      beta1 <- gbeta
      SK1 <- SK1 + 1
    } else {
      break
    }
  }
  
  convergence <- (SK1 <= itermax)
if(IC){
  V_Ind <- crossprod(X, sweep(X, 1, as.numeric(baseF * mu), "*")) / pphi
  Q_beta <- sum(newY1 * as.vector(X %*% gbeta) - exp(as.vector(X %*% gbeta))) / pphi
}else{V_Ind <- Q_beta <-  NULL}

  return(list(beta = gbeta, rho = rho, pphi = pphi, baseF = baseF, V_Ind = V_Ind, Q_beta= Q_beta, convergence = convergence))
}