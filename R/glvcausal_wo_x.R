#' Implement the Bayesian causal discovery algorithm 
#'
#' Implement the Bayesian causal discovery algorithm in the presence of cycles and unobserved confounders when no covariates are included
#'
#' @param data The observed data
#' @param mh_setup_lst A list containing the Metropolis-Hastings setups
#' @param init_lst A list containing the initial values of the parameters
#' @param prior_lst A list containing the prior hyper parameters
#' @param chain_setup_lst A list containing the MCMC chain setups
#' @param verbose A boolean indicating whether to print the progress of the MCMC chain
#'
#' @returns a list containing mcmc samples
#' @export
glvcausal_wo_x <- function(data, mh_setup_lst, init_lst, prior_lst, chain_setup_lst, verbose){

  set.seed(chain_setup_lst$seed)

  param_all_lst <- init_lst
  data$X <- matrix(0, nrow = nrow(data$Y), ncol = 2)
  param_all_lst$A <- matrix(0, nrow = ncol(data$Y), ncol = 2)
  mcmc_lst <- list()
  mcmc_lst_count <- 0

  for(it in 2:chain_setup_lst$Nit){

    param_all_lst$mu <- update_mu_rcpp(param_all_lst$B, param_all_lst$A, param_all_lst$L, param_all_lst$C, param_all_lst$tau, param_all_lst$sigma2, data$Y, data$X)

    param_all_lst$B <- update_B(param_all_lst$B, param_all_lst$gamma_beta, param_all_lst$nu_beta, param_all_lst$A, param_all_lst$mu, param_all_lst$L, param_all_lst$C, param_all_lst$sigma2, mh_setup_lst$B_step, data)
    param_all_lst$gamma_beta <- update_gamma_beta(param_all_lst$B, param_all_lst$nu_beta, param_all_lst$rho_beta, prior_lst)
    param_all_lst$nu_beta <- update_nu_beta(param_all_lst$B, param_all_lst$gamma_beta, prior_lst)
    param_all_lst$rho_beta <- update_rho_beta(param_all_lst$gamma_beta, prior_lst)

    param_all_lst$tau <- update_tau_rcpp(param_all_lst$B, param_all_lst$A, param_all_lst$L, param_all_lst$C, param_all_lst$mu, param_all_lst$sigma2, data$Y, data$X)
    param_all_lst$sigma2 <- c(update_sigma2_rcpp(param_all_lst$B, param_all_lst$A, param_all_lst$L, param_all_lst$C, param_all_lst$mu, param_all_lst$tau, data$Y, data$X, prior_lst$a_sigma, prior_lst$b_sigma))

    param_all_lst$C <- update_C_rcpp(param_all_lst$mu, param_all_lst$A, param_all_lst$B, param_all_lst$L, param_all_lst$tau, param_all_lst$sigma2, data$Y, data$X)

    param_all_lst$a[1] <- update_a1(param_all_lst$a[1], param_all_lst$a[2], param_all_lst$sparsity_matrix, param_all_lst$pivot, mh_setup_lst$a1_step, prior_lst)
    param_all_lst$a[2] <- update_a2(param_all_lst$a[1], param_all_lst$a[2], param_all_lst$sparsity_matrix, param_all_lst$pivot, mh_setup_lst$a2_step, prior_lst)
    param_all_lst$zeta <- update_zeta(param_all_lst$a[1], param_all_lst$a[2], param_all_lst$sparsity_matrix, param_all_lst$pivot, param_all_lst$P_star, prior_lst)
    param_all_lst$kappa <- update_kappa(param_all_lst$sparsity_matrix, param_all_lst$L, param_all_lst$sigma2, prior_lst)
    param_all_lst$L0 <- compute_L0(param_all_lst$kappa, param_all_lst$P_star, data)
    param_all_lst$L <- update_L_uglt_shrink(param_all_lst$mu, param_all_lst$A, param_all_lst$B, param_all_lst$C, param_all_lst$tau, param_all_lst$sigma2, param_all_lst$L0, param_all_lst$sparsity_matrix, param_all_lst$P_star, data)
    param_all_lst$sparsity_matrix <- update_sparsity_indicator(param_all_lst$sparsity_matrix, param_all_lst$pivot, param_all_lst$zeta, param_all_lst$C, param_all_lst$B, param_all_lst$A, param_all_lst$mu, param_all_lst$tau, param_all_lst$L0, param_all_lst$P_star, prior_lst, data)
    mcmc_tmp <- update_pivot(param_all_lst$sparsity_matrix, param_all_lst$pivot, param_all_lst$a[1], param_all_lst$a[2], param_all_lst$P_star, param_all_lst$C, param_all_lst$B, param_all_lst$A, param_all_lst$mu, param_all_lst$tau, param_all_lst$L0, prior_lst, mh_setup_lst, data)
    param_all_lst$pivot <- mcmc_tmp$pivot
    param_all_lst$sparsity_matrix <- mcmc_tmp$sparsity_matrix

    mcmc_tmp <- update_P_star(param_all_lst$mu, param_all_lst$A, param_all_lst$B, param_all_lst$L, param_all_lst$C, param_all_lst$sigma2, param_all_lst$a[1], param_all_lst$a[2], param_all_lst$sparsity_matrix, param_all_lst$pivot, param_all_lst$P_star, param_all_lst$kappa, mh_setup_lst, prior_lst, data)
    param_all_lst$P_star <- mcmc_tmp$mcmc_P_star
    param_all_lst$L <- mcmc_tmp$mcmc_L
    param_all_lst$C <- mcmc_tmp$mcmc_C
    param_all_lst$sigma2 <- mcmc_tmp$mcmc_sigma2
    param_all_lst$sparsity_matrix <- mcmc_tmp$mcmc_sparsity_matrix
    param_all_lst$pivot <- mcmc_tmp$mcmc_pivot

    if(save_check(it, chain_setup_lst$burn, chain_setup_lst$thin)){
      mcmc_lst_count <- mcmc_lst_count + 1
      mcmc_lst[[mcmc_lst_count]] <- save_param_res_wo_x(param_all_lst)
      }
    
    if(verbose){
      if(it %% 1000 == 0){
        cat("Iteration: ", it, "/", Nit, " completed \n")
      }
    }

  }
  
  return(mcmc_lst)

}
