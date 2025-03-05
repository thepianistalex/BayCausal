glvcausal_free_L <- function(data, mh_setup_lst, init_lst, prior_lst, chain_setup_lst, verbose){
  
  set.seed(chain_setup_lst$seed)
  
  param_all_lst <- init_lst
  data$X <- matrix(0, nrow = nrow(data$Y), ncol = 2)
  param_all_lst$A <- matrix(0, nrow = ncol(data$Y), ncol = 2)
  mcmc_all_lst <- list()
  mcmc_all_lst_count <- 0
  
  for(it in 2:chain_setup_lst$Nit){
    
    param_all_lst$mu <- update_mu_rcpp(param_all_lst$B, param_all_lst$A, param_all_lst$L_free, param_all_lst$C_free, param_all_lst$tau, param_all_lst$sigma2, data$Y, data$X)
    
    param_all_lst$B <- update_B(param_all_lst$B, param_all_lst$gamma_beta, param_all_lst$nu_beta, param_all_lst$A, param_all_lst$mu, param_all_lst$L_free, param_all_lst$C_free, param_all_lst$sigma2, mh_setup_lst$B_step, data)
    param_all_lst$gamma_beta <- update_gamma_beta(param_all_lst$B, param_all_lst$nu_beta, param_all_lst$rho_beta, prior_lst)
    param_all_lst$nu_beta <- update_nu_beta(param_all_lst$B, param_all_lst$gamma_beta, prior_lst)
    param_all_lst$rho_beta <- update_rho_beta(param_all_lst$gamma_beta, prior_lst)
    
    param_all_lst$L_free <- update_L_free_rcpp(param_all_lst$mu, param_all_lst$A, param_all_lst$B, param_all_lst$C_free, param_all_lst$tau, param_all_lst$sigma2, data$Y, data$X)
    param_all_lst$C_free <- update_C_rcpp(param_all_lst$mu, param_all_lst$A, param_all_lst$B, param_all_lst$L_free, param_all_lst$tau, param_all_lst$sigma2, data$Y, data$X)
    
    param_all_lst$tau <- update_tau_rcpp(param_all_lst$B, param_all_lst$A, param_all_lst$L_free, param_all_lst$C_free, param_all_lst$mu, param_all_lst$sigma2, data$Y, data$X)
    param_all_lst$sigma2 <- c(update_sigma2_rcpp(param_all_lst$B, param_all_lst$A, param_all_lst$L_free, param_all_lst$C_free, param_all_lst$mu, param_all_lst$tau, data$Y, data$X, prior_lst$a_sigma, prior_lst$b_sigma))
    
    if(save_check(it, chain_setup_lst$burn, chain_setup_lst$thin)){
      mcmc_all_lst_count <- mcmc_all_lst_count + 1
      mcmc_all_lst[[mcmc_all_lst_count]] <- param_all_lst
      mcmc_all_lst[[mcmc_all_lst_count]]$tau <- NULL
    }
    
    if(verbose){
      if(it %% 1000 == 0){
        cat("Iteration: ", it, "/", Nit, " completed for first step \n")
      }
    }
    
  }
  
  return(mcmc_all_lst)
  
}