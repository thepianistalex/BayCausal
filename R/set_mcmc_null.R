set_mcmc_default <- function(){

  mcmc_setup_lst <- list()

  mcmc_setup_lst$B_step <- NULL
  mcmc_setup_lst$a1_step <- NULL
  mcmc_setup_lst$a2_step <- NULL
  mcmc_setup_lst$pshift <- NULL
  mcmc_setup_lst$pswitch <- NULL
  mcmc_setup_lst$pa <- NULL
  mcmc_setup_lst$ps <- NULL

  return(mcmc_setup_lst)
}
