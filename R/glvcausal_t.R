#' Implement the Bayesian causal discovery algorithm with Student-t errors
#'
#' Implement the Bayesian causal discovery algorithm in the presence of cycles and unobserved confounders
#' assuming Student-t distributed errors
#'
#' @param data The observed data
#' @param mh_setup_lst A list containing the Metropolis-Hastings setups
#' @param init_lst A list containing the initial values of the parameters
#' @param prior_lst A list containing the prior hyper parameters. Must contain `nu_t` (degrees of freedom).
#' @param chain_setup_lst A list containing the MCMC chain setups
#' @param verbose A boolean indicating whether to print the progress of the MCMC chain
#'
#' @returns a list containing mcmc samples
#' @export
glvcausal_t <- function(data, mh_setup_lst, init_lst, prior_lst, chain_setup_lst, verbose){

  set.seed(chain_setup_lst$seed)

  param_all_lst <- init_lst
  mcmc_lst <- list()
  mcmc_lst_count <- 0
  
  if(is.null(prior_lst$nu_t)) {
    # Default nu_t if not provided, though recommended to provide it.
    # A large nu_t -> normal distribution. A small nu_t -> heavy tails.
    # Let's error if missing as per requirement "assuming a nu is fixed" (implies user choice)
    stop("prior_lst$nu_t is required for glvcausal_t")
  }

  for(it in 2:chain_setup_lst$Nit){

    param_all_lst$mu <- update_mu_rcpp(param_all_lst$B, param_all_lst$A, param_all_lst$L, param_all_lst$C, param_all_lst$tau, param_all_lst$sigma2, data$Y, data$X)

    param_all_lst$A <- update_A_rcpp(param_all_lst$gamma_alpha, param_all_lst$nu_alpha, param_all_lst$B, param_all_lst$L, param_all_lst$C, param_all_lst$mu, param_all_lst$tau, param_all_lst$sigma2, data$Y, data$X)
    param_all_lst$gamma_alpha <- update_gamma_alpha(param_all_lst$A, param_all_lst$nu_alpha, param_all_lst$rho_alpha, prior_lst)
    param_all_lst$nu_alpha <- update_nu_alpha(param_all_lst$A, param_all_lst$gamma_alpha, prior_lst)
    param_all_lst$rho_alpha <- update_rho_alpha(param_all_lst$gamma_alpha, prior_lst)

    param_all_lst$B <- update_B_t(param_all_lst$B, param_all_lst$gamma_beta, param_all_lst$nu_beta, param_all_lst$A, param_all_lst$mu, param_all_lst$L, param_all_lst$C, param_all_lst$sigma2, mh_setup_lst$B_step, data, prior_lst$nu_t)
    param_all_lst$gamma_beta <- update_gamma_beta(param_all_lst$B, param_all_lst$nu_beta, param_all_lst$rho_beta, prior_lst)
    param_all_lst$nu_beta <- update_nu_beta(param_all_lst$B, param_all_lst$gamma_beta, prior_lst)
    param_all_lst$rho_beta <- update_rho_beta(param_all_lst$gamma_beta, prior_lst)

    # Use student t update for tau
    param_all_lst$tau <- update_tau_t_rcpp(param_all_lst$B, param_all_lst$A, param_all_lst$L, param_all_lst$C, param_all_lst$mu, param_all_lst$sigma2, data$Y, data$X, prior_lst$nu_t)
    
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

    mcmc_tmp <- update_P_star_t(param_all_lst$mu, param_all_lst$A, param_all_lst$B, param_all_lst$L, param_all_lst$C, param_all_lst$sigma2, param_all_lst$a[1], param_all_lst$a[2], param_all_lst$sparsity_matrix, param_all_lst$pivot, param_all_lst$P_star, param_all_lst$kappa, mh_setup_lst, prior_lst, data, prior_lst$nu_t)
    param_all_lst$P_star <- mcmc_tmp$mcmc_P_star
    param_all_lst$L <- mcmc_tmp$mcmc_L
    param_all_lst$C <- mcmc_tmp$mcmc_C
    param_all_lst$sigma2 <- mcmc_tmp$mcmc_sigma2
    param_all_lst$sparsity_matrix <- mcmc_tmp$mcmc_sparsity_matrix
    param_all_lst$pivot <- mcmc_tmp$mcmc_pivot

    if(save_check(it, chain_setup_lst$burn, chain_setup_lst$thin)){
      mcmc_lst_count <- mcmc_lst_count + 1
      mcmc_lst[[mcmc_lst_count]] <- save_param_res(param_all_lst)
      }
    
    if(verbose){
      if(it %% 1000 == 0){
        cat("Iteration: ", it, "/", chain_setup_lst$Nit, " completed \n")
      }
    }

  }
  
  return(mcmc_lst)

}
