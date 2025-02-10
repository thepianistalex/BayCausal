generate_sim_data <- function(seed, n){
  
  # For the non-random part of the data, use same seed all the time
  set.seed(0)
  Q <- 5
  S <- 2
  P <- 2
  
  mu <- runif(Q, -1, 1)
  
  B <- matrix(c(
    0, 0.5, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0.4, -0.7,
    0.3, 0, 0, 0, 0,
    0, 0, 0, 0, 0
  ), nrow = Q, byrow = TRUE)
  
  A <- matrix(c(runif(Q*S, -1, 1)), Q, S)
  
  L <- matrix(c(
    0, 0,
    0.5, 0,
    0.3, 0,
    0, -0.5,
    0, 0.4), ncol = 2, byrow = TRUE)
  sparsity_matrix <- ifelse(L != 0, 1, 0)
  pivot <- numeric(ncol(L))
  for (p in 1:ncol(L)) {
    pivot[p] <- ifelse(any(L[, p] != 0), which(L[, p] != 0)[1], NA)
  }
  
  # For random part, use the seed
  set.seed(seed)
  X <- matrix(rnorm(n*S), n, S)
  C <- matrix(rnorm(n*P), n, P)
  
  sigma_e <- rep(sqrt(1/8), Q)
  E <- matrix(NA, n, Q)
  for(q in 1:Q){
    E[,q] <- LaplacesDemon::rlaplace(n, 0, 2*sigma_e[q])
  }
  
  Y <- matrix(NA, n, Q)
  mix_mat <- solve(diag(Q)-B)
  for(i in 1:n){
    Y[i,] <- mix_mat %*% (mu + A%*%X[i,] + L%*%C[i,] + E[i,])
  }
  
  data <- list(Y = Y, X = X, mu = mu, B = B, A = A, L = L, C = C, sigma_e = sigma_e, Q = Q, S = S, P = P, n = n)

  return(data)
}



are_all_close <- function(v, w, abs_tol) {
  abs_diff <- abs(v - w)
  are_all_within_atol <- all(abs_diff < abs_tol)
  return(are_all_within_atol)
}



glvcausal_check_mu <- function(data, mh_setup_lst, init_lst, prior_lst, chain_setup_lst, verbose){
  
  set.seed(chain_setup_lst$seed)
  
  param_all_lst <- init_lst
  mcmc_lst <- list()
  mcmc_lst_count <- 0
  
  for(it in 2:chain_setup_lst$Nit){
    
    param_all_lst$mu <- update_mu_rcpp(param_all_lst$B, param_all_lst$A, param_all_lst$L, param_all_lst$C, param_all_lst$tau, param_all_lst$sigma2, data$Y, data$X)
    param_all_lst$tau <- update_tau_rcpp(param_all_lst$B, param_all_lst$A, param_all_lst$L, param_all_lst$C, param_all_lst$mu, param_all_lst$sigma2, data$Y, data$X)
    
    if(save_check(it, chain_setup_lst$burn, chain_setup_lst$thin)){
      mcmc_lst_count <- mcmc_lst_count + 1
      mcmc_lst[[mcmc_lst_count]] <- save_param_res(param_all_lst)
    }
    
  }
  
  return(mcmc_lst)
  
}



glvcausal_check_A <- function(data, mh_setup_lst, init_lst, prior_lst, chain_setup_lst, verbose){
  
  set.seed(chain_setup_lst$seed)
  
  param_all_lst <- init_lst
  mcmc_lst <- list()
  mcmc_lst_count <- 0
  
  for(it in 2:chain_setup_lst$Nit){
    
    param_all_lst$A <- update_A_rcpp(param_all_lst$gamma_alpha, param_all_lst$nu_alpha, param_all_lst$B, param_all_lst$L, param_all_lst$C, param_all_lst$mu, param_all_lst$tau, param_all_lst$sigma2, data$Y, data$X)
    param_all_lst$gamma_alpha <- update_gamma_alpha(param_all_lst$A, param_all_lst$nu_alpha, param_all_lst$rho_alpha, prior_lst)
    param_all_lst$nu_alpha <- update_nu_alpha(param_all_lst$A, param_all_lst$gamma_alpha, prior_lst)
    param_all_lst$rho_alpha <- update_rho_alpha(param_all_lst$gamma_alpha, prior_lst)
    
    if(save_check(it, chain_setup_lst$burn, chain_setup_lst$thin)){
      mcmc_lst_count <- mcmc_lst_count + 1
      mcmc_lst[[mcmc_lst_count]] <- save_param_res(param_all_lst)
    }
    
  }
  
  return(mcmc_lst)
  
}