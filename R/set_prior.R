set_prior <- function(data){

  prior_lst <- list()

  prior_lst$a_nu <- 1
  prior_lst$b_nu <- 1
  prior_lst$a_rho <- 1
  prior_lst$b_rho <- 1
  prior_lst$a_sigma <- 1
  prior_lst$b_sigma <- 1
  prior_lst$nu_0 <- 2.5e-4
  prior_lst$a_a1 <- 6
  Eq <- 1
  prior_lst$H <- ncol(data$Y) - 1
  prior_lst$b_a1 <- a_a1*(H-Eq)/(H*Eq)
  prior_lst$a_a2 <- a_a1
  prior_lst$b_a2 <- a_a1
  prior_lst$a_kappa <- 1
  prior_lst$b_kappa <- 1

  return(prior_lst)
}
