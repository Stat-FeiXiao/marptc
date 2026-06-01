#' Title Estimation the covariate coefficients using the QIF approach
#'
#' @description Estimation the covariate coefficients using the QIF approach.
#'
#' @param Time right censored data which is the follow up time.
#' @param Status right censored data which is the follow up time.
#' @param X a matrix of covariates that that may have effect on the failure times.
#' @param beta the initial value of the covariate coefficients.
#' @param id a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.
#' @param IC if \code{TRUE}, the function returns the QIC for the GEE method or the BIQIF for the QIF method. Default is \code{FALSE}.
#' @param eps tolerance for convergence. The default is \code{eps = 1e-6}. Iteration stops once the relative change in deviance is less than \code{eps}.
#' @param corstr a character string specifying the correlation structure. The following are permitted: \code{independence}, \code{exchangeable} and \code{AR1}.
#' @param itermax specifies the maximum iteration number. If the convergence criterion is not met, the iteration will be stopped after \code{itermax} iterations and
#' the estimates will be based on the last iteration.
#'
#' @export qeq

qeq <- function(Time, Status, X, beta, id, IC =FALSE, eps = 1e-6, corstr, itermax) {

  K <- length(unique(id))
  n <- as.vector(table(id))
  n_mean <- mean(n) 
  qbeta <- beta1 <- beta
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
  
  # 【优化4】：提前缓存切片数据，极大地加速 5000 次循环的速度
  X_list <- lapply(cluster_idx, function(idx) X[idx, , drop = FALSE])
  Status_list <- lapply(cluster_idx, function(idx) Status[idx])
  
  mu <- exp(as.vector(X %*% qbeta))
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

  dim_multiplier <- if (corstr == "independence") 1 else 2
  G_dim <- p * dim_multiplier

  S1_init <- Status - mu * baseF
  C_frozen <- matrix(0, G_dim, G_dim)
  
  for (i in seq_len(K)) {
    n_i <- n[i]
    w_i <- sqrt(n_mean / n_i) 
    
    X_i <- X_list[[i]]
    S1_i <- Status_list[[i]] - mu[cluster_idx[[i]]] * baseF[cluster_idx[[i]]]
    sqrt_mu_i <- sqrt(mu[cluster_idx[[i]]])
    
    g_1i <- t(X_i) %*% S1_i
    
    if (corstr == "independence") {
      g_i <- g_1i
    } else {
      A <- t(X_i * sqrt_mu_i)
      if (corstr == "exchangeable") {
        part2 <- rowSums(A) - A
        if (n_i > 1) part2 <- part2 / (n_i - 1)
      } else if (corstr == "AR1") {
        # 【优化1】：零矩阵内存消耗的 AR1 向量化计算
        if (n_i == 1) {
          part2 <- matrix(0, nrow = p, ncol = 1)
        } else {
          part2 <- matrix(0, nrow = p, ncol = n_i)
          part2[, 1:(n_i-1)] <- part2[, 1:(n_i-1)] + A[, 2:n_i, drop = FALSE]
          part2[, 2:n_i] <- part2[, 2:n_i] + A[, 1:(n_i-1), drop = FALSE]
        }
      }
      g_2i <- part2 %*% (S1_i / sqrt_mu_i)
      g_i <- rbind(g_1i, g_2i)
    }
    
    g_i_weighted <- g_i * w_i
    C_frozen <- C_frozen + g_i_weighted %*% t(g_i_weighted)
  }
  
  trace_C <- sum(diag(C_frozen))
  ridge_lambda <- max(1e-4, trace_C / (G_dim * K)) 
  C_frozen <- C_frozen + diag(ridge_lambda, G_dim)
  
  C_inv_frozen <- tryCatch(solve(C_frozen), error = function(e) MASS::ginv(C_frozen))

  SK1 <- 1
  repeat {
    mu <- exp(as.vector(X %*% qbeta))
    S1 <- Status - mu * baseF
    
    G <- matrix(0, G_dim, 1)
    G1 <- matrix(0, G_dim, p)
    
    for (i in seq_len(K)) {
      idx <- cluster_idx[[i]]
      X_i <- X_list[[i]]
      mu_i <- mu[idx]
      sqrt_mu_i <- sqrt(mu_i)
      baseF_i <- baseF[idx]
      S1_i <- S1[idx]
      
      n_i <- n[i]
      w_i <- sqrt(n_mean / n_i)
      
      g_1i <- t(X_i) %*% S1_i
      va_1i <- - t(X_i * (baseF_i * mu_i)) %*% X_i
      
      if (corstr == "independence") {
        g_i <- g_1i
        va_i <- va_1i
      } else {
        A <- t(X_i * sqrt_mu_i)
        
        if (corstr == "exchangeable") {
          part2 <- rowSums(A) - A
          if (n_i > 1) part2 <- part2 / (n_i - 1)
        } else if (corstr == "AR1") {
          if (n_i == 1) {
            part2 <- matrix(0, nrow = p, ncol = 1)
          } else {
            part2 <- matrix(0, nrow = p, ncol = n_i)
            part2[, 1:(n_i-1)] <- part2[, 1:(n_i-1)] + A[, 2:n_i, drop = FALSE]
            part2[, 2:n_i] <- part2[, 2:n_i] + A[, 1:(n_i-1), drop = FALSE]
          }
        }
        
        g_2i <- part2 %*% (S1_i / sqrt_mu_i)
        part2_scaled <- sweep(part2, 2, baseF_i * sqrt_mu_i, "*")
        va_2i <- - (part2_scaled %*% X_i)
        
        g_i <- rbind(g_1i, g_2i)
        va_i <- rbind(va_1i, va_2i)
      }
      
      G <- G + (g_i * w_i)
      G1 <- G1 + (va_i * w_i)
    }
    
    Q1 <- t(G1) %*% C_inv_frozen %*% G
    Q2 <- t(G1) %*% C_inv_frozen %*% G1
    
    # 【优化3】：强制对称并加载微小对角 Ridge 惩罚，保证 GEE 级别的稳定求逆
    Q2 <- (Q2 + t(Q2)) / 2 
    Q2 <- Q2 + diag(1e-8, p) 
    
    delta_beta <- tryCatch(solve(Q2, Q1), error = function(e) MASS::ginv(Q2) %*% Q1)
    
    # 【优化2】：牛顿步长衰减保护 (Step-Halving Line Search)
    step_size <- 1
    max_halving <- 8
    success <- FALSE
    
    for (step in 1:max_halving) {
      qbeta_try <- qbeta - step_size * delta_beta
      mu_try <- exp(as.vector(X %*% qbeta_try))
      
      # 防范指数爆炸（溢出）
      if (any(mu_try > 1e10) || any(is.infinite(mu_try))) {
        step_size <- step_size / 2
        next
      }
      
      mu_sort_try <- mu_try[ord]
      cum_mu_pad_try <- c(0, cumsum(mu_sort_try))
      dRbet_try <- sum(mu_try) - cum_mu_pad_try[first_idx]
      dRbet_d_try <- dRbet_try[Status == 1]
      Rbet.min_try <- min(dRbet_d_try)
      interval_try <- c(Rbet.min_try - sum(Status), Rbet.min_try - 1)
      
      lambet_try <- tryCatch({
        uniroot(getlambda, interval_try, tol = .Machine$double.eps^0.75, dRbet = dRbet_d_try)$root
      }, error = function(e) NULL)
      
      if (is.null(lambet_try)) {
        step_size <- step_size / 2
      } else {
        qbeta <- qbeta_try
        mu <- mu_try
        lambet <- lambet_try
        dRbet_d <- dRbet_d_try
        success <- TRUE
        break
      }
    }
    
    if (!success) break # 如果连续衰减仍然失败，安全退出并返回目前最优解
    
    term <- 1 / (dRbet_d - lambet)
    term_sort <- term[ord1]
    cum_term_pad <- c(0, cumsum(term_sort))
    baseF <- cum_term_pad[idx_interval + 1]
    
    if (any(abs(qbeta - beta1) > eps) && SK1 <= itermax) {
      beta1 <- qbeta
      SK1 <- SK1 + 1
    } else {
      break
    }
  }

  convergence <- (SK1 <= itermax) && success
if(IC){
  Quad <- t(G) %*% C_inv_frozen %*% G
}else{ Quad <-NULL}
  return(list(beta = qbeta, baseF = baseF, convergence = convergence, Quad = Quad))
}