set_mcmc_default <- function(){

  mcmc_setup_lst <- list()

  mcmc_setup_lst$B_step <- 0.02
  mcmc_setup_lst$a1_step <- 0.1
  mcmc_setup_lst$a2_step <- 0.1
  mcmc_setup_lst$pshift <- 1/3
  mcmc_setup_lst$pswitch <- 1/3
  mcmc_setup_lst$pa <- 0.5
  mcmc_setup_lst$ps <- 0.5

  return(mcmc_setup_lst)
}
