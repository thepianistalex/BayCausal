update_gamma_alpha <- function(A, nu_alpha, rho_alpha, prior_lst){

  Q <- nrow(A)
  S <- ncol(A)

  gamma_alpha_update <- matrix(NA, Q, S)
  nu_0 <- prior_lst$nu_0

  for (q in 1:Q){
    for (s in 1:S){
      # prob = pr(gamma_alpha_lqs=1|~) / pr(gamma_alpha_lqs=nu_0|~)
      prob <- (sqrt(nu_0)*rho_alpha)/(1-rho_alpha)*exp((1-nu_0)*A[q,s]^2/(2*nu_0*nu_alpha[q,s]))

      if (prob == Inf){
        gamma_alpha_update[q,s] <- 1
      }else{
        gamma_alpha_update[q,s] <- sample(c(1,nu_0), 1, prob=c(prob, 1))
      }

    }
  }

  return(gamma_alpha_update)
}



update_nu_alpha <- function(A, gamma_alpha, prior_lst, data){

  Q <- nrow(A)
  S <- ncol(A)

  nu_alpha_update <- matrix(NA, Q, S)
  a_nu <- prior_lst$a_nu
  b_nu <- prior_lst$b_nu

  for (q in 1:Q){
    for (s in 1:S){
      a_nu_star <- a_nu + 1/2
      b_nu_star <- b_nu + A[q,s]^2/(2*gamma_alpha[q,s])
      nu_alpha_update[q,s] <- rinvgamma(1, a_nu_star, b_nu_star)
    }
  }

  return(nu_alpha_update)
}



update_rho_alpha <- function(gamma_alpha, prior_lst){

  a_rho <- prior_lst$a_rho
  b_rho <- prior_lst$b_rho
  nu_0 <- prior_lst$nu_0

  a_rho_star <- a_rho + sum(gamma_alpha == 1)
  b_rho_star <- b_rho + sum(gamma_alpha == nu_0)
  rho_alpha_update <- rbeta(1, a_rho_star, b_rho_star)

  return(rho_alpha_update)
}
