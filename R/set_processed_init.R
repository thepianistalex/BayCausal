set_processed_init <- function(data, mh_setup_lst, init_lst, prior_lst, 
                               chain_setup_lst_1, chain_setup_lst_2, verbose){
  
  init_processed <- init_lst
  res <- glvcausal_free_L(data, mh_setup_lst, init_processed, prior_lst, chain_setup_lst_1, verbose)
  init_processed <- update_init_wo_UGLT(init_processed, res, prior_lst)
  init_processed <- glvcausal_only_UGLT(data, mh_setup_lst, init_processed, prior_lst, chain_setup_lst_2, verbose)
  
  return(init_processed)
  
}



set_processed_init_wo_x <- function(data, mh_setup_lst, init_lst, prior_lst, 
                               chain_setup_lst_1, chain_setup_lst_2, verbose){
  
  init_processed <- init_lst
  res <- glvcausal_free_L_wo_x(data, mh_setup_lst, init_processed, prior_lst, chain_setup_lst_1, verbose)
  init_processed <- update_init_wo_UGLT_wo_x(init_processed, res, prior_lst)
  init_processed <- glvcausal_only_UGLT_wo_x(data, mh_setup_lst, init_processed, prior_lst, chain_setup_lst_2, verbose)
  
  return(init_processed)
  
}



update_init_wo_UGLT <- function(init_lst, res, prior_lst){
  
  init_lst$mu <- extract_post_mean(res, "mu", init_lst$P_star)
  
  gamma_alpha_tmp <- extract_post_mean(res, "gamma_alpha", init_lst$P_star)
  gamma_alpha_tmp <- ifelse(gamma_alpha_tmp > 0.5, 1, prior_lst$nu_0)
  init_lst$gamma_alpha <- gamma_alpha_tmp
  init_lst$nu_alpha <- extract_post_mean(res, "nu_alpha", init_lst$P_star)
  init_lst$rho_alpha <- extract_post_mean(res, "rho_alpha", init_lst$P_star)
  A_tmp <- extract_post_mean(res, "A", init_lst$P_star) * init_lst$gamma_alpha
  init_lst$A <- A_tmp * gamma_alpha_tmp
  
  gamma_beta_tmp <- extract_post_mean(res, "gamma_beta", init_lst$P_star)
  gamma_beta_tmp <- ifelse(gamma_beta_tmp > 0.5, 1, prior_lst$nu_0)
  for(q in 1:ncol(gamma_beta_tmp)){
    gamma_beta_tmp[q,q] <- 0
  }
  init_lst$gamma_beta <- gamma_beta_tmp
  init_lst$nu_beta <- extract_post_mean(res, "nu_beta", init_lst$P_star)
  init_lst$rho_beta <- extract_post_mean(res, "rho_beta", init_lst$P_star)
  B_tmp <- extract_post_mean(res, "B", init_lst$P_star)
  init_lst$B <- B_tmp * gamma_beta_tmp
  
  init_lst$sigma2 <- extract_post_mean(res, "sigma2", init_lst$P_star)
  
  return(init_lst)
}



update_init_wo_UGLT_wo_x <- function(init_lst, res, prior_lst){
  
  init_lst$mu <- extract_post_mean(res, "mu", init_lst$P_star)
  
  gamma_beta_tmp <- extract_post_mean(res, "gamma_beta", init_lst$P_star)
  gamma_beta_tmp <- ifelse(gamma_beta_tmp > 0.5, 1, prior_lst$nu_0)
  for(q in 1:ncol(gamma_beta_tmp)){
    gamma_beta_tmp[q,q] <- 0
  }
  init_lst$gamma_beta <- gamma_beta_tmp
  init_lst$nu_beta <- extract_post_mean(res, "nu_beta", init_lst$P_star)
  init_lst$rho_beta <- extract_post_mean(res, "rho_beta", init_lst$P_star)
  B_tmp <- extract_post_mean(res, "B", init_lst$P_star)
  init_lst$B <- B_tmp * gamma_beta_tmp
  
  init_lst$sigma2 <- extract_post_mean(res, "sigma2", init_lst$P_star)
  
  return(init_lst)
}