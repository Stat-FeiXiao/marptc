#' Title  Estimate Variance of coefficients estimation from the GEE approach with a sandwich formula
#'
#' @description Calculate the variance estimates using the sandwich formula.
#'
#' @param Time right censored data which is the follow up time.
#' @param Status the censoring indicator, normally 0 = event of interest happens, and 0 = censoring
#' @param X a matrix of covariates that may have effect on the failure times.
#' @param id a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.
#' @param beta covariate coefficients estimation from GEE method.
#' @param corstr a character string specifying the correlation structure. The following are permitted: \code{independence}, \code{exchangeable} and \code{AR1}.
#'
#' @export vargee

vargee <- function(Time, Status, X, id, beta, corstr) {
  
  K <- length(unique(id))
  n <- as.vector(table(id))
  p <- ncol(X)
  
  mu <- exp(as.vector(X %*% beta))
  
  ord <- order(Time)
  Time_sort <- Time[ord]
  first_idx <- match(Time, Time_sort)
  
  idx1 <- which(Status == 1)
  Time1 <- Time[idx1]
  ord1 <- order(Time1)
  Time1_sort <- Time1[ord1]
  idx_interval <- findInterval(Time, Time1_sort)
  
  mu_sort <- mu[ord]
  cum_mu_pad <- c(0, cumsum(mu_sort))
  dRbet <- sum(mu) - cum_mu_pad[first_idx]
  
  dRbet_d <- dRbet[Status == 1]
  Rbet.min <- min(dRbet_d)
  interval <- c(Rbet.min - sum(Status), Rbet.min - 1)
  lambet <- uniroot(getlambda, interval, tol = .Machine$double.eps^0.75, dRbet = dRbet_d)$root
  
  term <- 1 / (dRbet_d - lambet)
  term_sort <- term[ord1]
  cum_term_pad <- c(0, cumsum(term_sort))
  baseF <- cum_term_pad[idx_interval + 1]
  
  baseF1 <- baseF
  baseF1[baseF1 == 0] <- 1e-8
  
  newY1 <- Status / baseF1
  S1 <- Status - mu * baseF
  
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
    # 【核心修改 2】：动态限制 exchangeable 结构的 rho 下界，防止三明治方差计算崩溃
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
    rho <- max(min(rho, 0.95), -0.95) # AR1 维持原样
  }
  
  if (corstr == "exchangeable") {
    a_exch <- 1 / (1 - rho)
    b_exch_vec <- -rho / ((1 - rho) * (1 - rho + n * rho))
  } else if (corstr == "AR1") {
    invQC_list <- vector("list", K)
    c_ar1 <- 1 / (1 - rho^2)
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
  
  cluster_idx <- split(seq_along(id), id)
  VA1 <- matrix(0, p, p)
  fdm <- matrix(0, p, p)
  
  for (i in seq_len(K)) {
    idx <- cluster_idx[[i]]
    X_i <- X[idx, , drop = FALSE]
    mu_i <- mu[idx]
    sqrt_mu_i <- sqrt(mu_i)
    baseF_i <- baseF[idx]
    R_i <- S1[idx]
    
    R_scaled <- R_i / sqrt_mu_i
    RX_scaled <- X_i * R_scaled
    mu_base_scaled <- X_i * (baseF_i * sqrt_mu_i)
    
    if (corstr == "independence") {
      v <- R_i
      T2 <- X_i * R_i
      T3 <- mu_base_scaled * sqrt_mu_i
    } else if (corstr == "exchangeable") {
      S_R <- sum(R_scaled)
      v <- a_exch * R_i + b_exch_vec[i] * sqrt_mu_i * S_R
      
      S_RX <- colSums(RX_scaled)
      T2 <- a_exch * (X_i * R_i) + outer(sqrt_mu_i, b_exch_vec[i] * S_RX)
      
      S_mu_base <- colSums(mu_base_scaled)
      T3 <- a_exch * (mu_base_scaled * sqrt_mu_i) + outer(sqrt_mu_i, b_exch_vec[i] * S_mu_base)
    } else if (corstr == "AR1") {
      M_inv <- invQC_list[[i]]
      v <- sqrt_mu_i * as.vector(M_inv %*% R_scaled)
      T2 <- sqrt_mu_i * (M_inv %*% RX_scaled)
      T3 <- sqrt_mu_i * (M_inv %*% mu_base_scaled)
    }
    
    T1 <- X_i * v
    inner_term <- 0.5 * (T1 - T2) - T3
    
    VA1_i <- (1 / pphi) * (t(X_i) %*% inner_term)
    fdv_i <- (1 / pphi) * (t(X_i) %*% v)
    
    VA1 <- VA1 + VA1_i
    fdm <- fdm + fdv_i %*% t(fdv_i)
  }
  
  M <- -VA1
  inv_M <- tryCatch(solve(M), error = function(e) MASS::ginv(M))
  vcm <- inv_M %*% fdm %*% t(inv_M)
  
  var_beta <- diag(as.matrix(vcm))
  sd_beta <- sqrt(pmax(var_beta, 0)) 
  
  return(list(var_beta = var_beta, sd_beta = sd_beta, var_beta_matrix = vcm))
}