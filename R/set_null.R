set_init_null <- function(){
  
  init_list <- list(A=NULL, gamma_alpha=NULL, nu_alpha=NULL, rho_alpha=NULL,
                    B=NULL, gamma_beta=NULL, nu_beta=NULL, rho_beta=NULL,
                    C=NULL, L=NULL, a=NULL, tau=NULL, mu=NULL, sigma2=NULL, 
                    sparsity_matrix=NULL, pivot=NULL, zeta=NULL, kappa=NULL)
  
  return(init_list)    
}



set_prior_null <- function(data){
  
  prior_lst <- list()
  
  prior_lst$a_nu <- NULL
  prior_lst$b_nu <- NULL
  prior_lst$a_rho <- NULL
  prior_lst$b_rho <- NULL
  prior_lst$a_sigma <- NULL
  prior_lst$b_sigma <- NULL
  prior_lst$nu_0 <- NULL
  prior_lst$a_a1 <- NULL
  prior_lst$H <- NULL
  prior_lst$b_a1 <- NULL
  prior_lst$a_a2 <- NULL
  prior_lst$b_a2 <- NULL
  prior_lst$a_kappa <- NULL
  prior_lst$b_kappa <- NULL
  
  return(prior_lst)
}



set_mh_null <- function(){
  
  mh_setup_lst <- list()
  
  mh_setup_lst$B_step <- NULL
  mh_setup_lst$a1_step <- NULL
  mh_setup_lst$a2_step <- NULL
  mh_setup_lst$pshift <- NULL
  mh_setup_lst$pswitch <- NULL
  mh_setup_lst$pa <- NULL
  mh_setup_lst$ps <- NULL
  
  return(mh_setup_lst)
}



set_chain_null <- function(data){
  
  chain_setup_lst <- list()
  
  chain_setup_lst$Nit <- NULL
  chain_setup_lst$burn <- NULL
  chain_setup_lst$thin <- NULL
  chain_setup_lst$seed <- NULL
  
  return(chain_setup_lst)
}
