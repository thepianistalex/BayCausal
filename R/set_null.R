set_init_null <- function(){
  
  init_list <- list(A=NA, gamma_alpha=NA, nu_alpha=NA, rho_alpha=NA,
                    B=NA, gamma_beta=NA, nu_beta=NA, rho_beta=NA,
                    C=NA, L=NA, a=NA, tau=NA, mu=NA, sigma2=NA, 
                    sparsity_matrix=NA, pivot=NA, zeta=NA, kappa=NA)
  
  return(init_list)    
}



set_prior_null <- function(data){
  
  prior_lst <- list()
  
  prior_lst$a_nu <- NA
  prior_lst$b_nu <- NA
  prior_lst$a_rho <- NA
  prior_lst$b_rho <- NA
  prior_lst$a_sigma <- NA
  prior_lst$b_sigma <- NA
  prior_lst$nu_0 <- NA
  prior_lst$a_a1 <- NA
  prior_lst$H <- NA
  prior_lst$b_a1 <- NA
  prior_lst$a_a2 <- NA
  prior_lst$b_a2 <- NA
  prior_lst$a_kappa <- NA
  prior_lst$b_kappa <- NA
  
  return(prior_lst)
}



set_mh_null <- function(){
  
  mh_setup_lst <- list()
  
  mh_setup_lst$B_step <- NA
  mh_setup_lst$a1_step <- NA
  mh_setup_lst$a2_step <- NA
  mh_setup_lst$pshift <- NA
  mh_setup_lst$pswitch <- NA
  mh_setup_lst$pa <- NA
  mh_setup_lst$ps <- NA
  
  return(mh_setup_lst)
}



#' Title
#'
#' @returns a list containing chain setup
#' @export
set_chain <- function(Nit, burn, thin, seed){
  
  chain_setup_lst <- list()
  
  chain_setup_lst$Nit <- Nit
  chain_setup_lst$burn <- burn
  chain_setup_lst$thin <- thin
  chain_setup_lst$seed <- seed
  
  return(chain_setup_lst)
}
