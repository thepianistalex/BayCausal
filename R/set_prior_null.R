set_prior_default <- function(data){

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
