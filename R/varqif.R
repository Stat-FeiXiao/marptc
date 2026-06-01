#' Title  Estimate Variance of coefficients estimation from the QIF approach with a sandwich formula
#'
#' @description Calculate the variance estimates using the sandwich formula.
#'
#' @param Time right censored data which is the follow up time.
#' @param Status the censoring indicator, normally 0 = event of interest happens, and 0 = censoring
#' @param X a matrix of covariates that may have effect on the failure times.
#' @param id a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.
#' @param beta covariate coefficients estimation from QIF method.
#' @param corstr a character string specifying the correlation structure. The following are permitted: \code{independence}, \code{exchangeable} and \code{AR1}.
#'
#' @export varqif

varqif <- function(Time, Status, X, id, beta, corstr) {
  
  K <- length(unique(id))
  n <- as.vector(table(id))
  p <- ncol(X)
  
  mu <- exp(as.vector(X %*% beta))
  
  ord <- order(Time)
  Time_sort <- Time[ord]
  first_idx <- match(Time, Time_sort)
  
  mu_sort <- mu[ord]
  cum_mu_pad <- c(0, cumsum(mu_sort))
  dRbet <- sum(mu) - cum_mu_pad[first_idx]
  
  dRbet_d <- dRbet[Status == 1]
  Rbet.min <- min(dRbet_d)
  interval <- c(Rbet.min - sum(Status), Rbet.min - 1)
  lambet <- uniroot(getlambda, interval, tol = .Machine$double.eps^0.75, dRbet = dRbet_d)$root
  
  idx1 <- which(Status == 1)
  Time1 <- Time[idx1]
  term <- 1 / (dRbet_d - lambet)
  ord1 <- order(Time1)
  Time1_sort <- Time1[ord1]
  term_sort <- term[ord1]
  cum_term_pad <- c(0, cumsum(term_sort))
  idx_interval <- findInterval(Time, Time1_sort)
  baseF <- cum_term_pad[idx_interval + 1]
  
  baseF1 <- baseF
  baseF1[baseF1 == 0] <- 1e-8
  S1 <- Status - mu * baseF
  
  dim_multiplier <- if (corstr == "independence") 1 else 2
  G_dim <- p * dim_multiplier
  
  VA1 <- matrix(0, G_dim, p)
  C <- matrix(0, G_dim, G_dim)
  
  cluster_idx <- split(seq_along(id), id)
  
  M2_list <- vector("list", K)
  if (corstr == "AR1") {
    for (i in seq_len(K)) {
      ni <- n[i]
      M2 <- matrix(0, ni, ni)
      if (ni > 1) {
        idx_off <- abs(row(M2) - col(M2)) == 1
        M2[idx_off] <- 1
      }
      M2_list[[i]] <- M2
    }
  }
  
  for (i in seq_len(K)) {
    idx <- cluster_idx[[i]]
    X_i <- X[idx, , drop = FALSE]
    mu_i <- mu[idx]
    sqrt_mu_i <- sqrt(mu_i)
    baseF_i <- baseF[idx]
    S1_i <- S1[idx]
    
    g_1i <- t(X_i) %*% S1_i
    va_1i <- - t(X_i * (baseF_i * mu_i)) %*% X_i
    
    if (corstr == "independence") {
      g_i <- g_1i
      va_i <- va_1i
    } else {
      A <- t(X_i * sqrt_mu_i)
      R_star <- S1_i / sqrt_mu_i
      C_star <- X_i * (baseF_i * sqrt_mu_i)
      
      if (corstr == "exchangeable") {
        g_2i <- rowSums(A) * sum(R_star) - A %*% R_star
        va_2i <- - (outer(rowSums(A), colSums(C_star)) - A %*% C_star)
      } else if (corstr == "AR1") {
        M2 <- M2_list[[i]]
        part2 <- A %*% M2
        g_2i <- part2 %*% R_star
        va_2i <- - (part2 %*% C_star)
      }
      
      g_i <- rbind(g_1i, g_2i)
      va_i <- rbind(va_1i, va_2i)
    }
    
    C <- C + g_i %*% t(g_i)
    VA1 <- VA1 + va_i
  }
  
  C_inv <- tryCatch(solve(C), error = function(e) MASS::ginv(C))
  M_neg <- t(VA1) %*% C_inv %*% VA1
  vcm <- tryCatch(solve(M_neg), error = function(e) MASS::ginv(M_neg))
  
  var_beta <- diag(as.matrix(vcm))
  sd_beta <- sqrt(pmax(var_beta, 0)) 
  
  return(list(var_beta = var_beta, sd_beta = sd_beta, var_beta_matrix = vcm))
}