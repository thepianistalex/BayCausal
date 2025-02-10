update_B <- function(B, gamma_beta, nu_beta, A, mu, L, C, sigma2, mh_step_size, data){

  Y <- data$Y
  X <- data$X
  Q <- ncol(Y)

  B_update <- B

  for (q in 1:Q){
    for (p in 1:Q){
      if (p == q){
        B_update[q,p] <- 0
      }else{
        # propose the new beta_qp
        B_new <- B_update
        B_new[q,p] <- rnorm(1, B_update[q,p], mh_step_size)

        # jacobian matrix norm constant
        if (max(Mod(eigen(B_new)$values)) < 1){
          # calculate the acceptance ratio
          ratio <- compute_loglik_B_laplace_faster_rcpp(q, B_new, A, L, C, mu, sigma2, Y, X) -
            compute_loglik_B_laplace_faster_rcpp(q, B_update, A, L, C, mu, sigma2, Y, X) +
            dnorm(x=B_new[q,p], 0, sqrt(gamma_beta[q,p]*nu_beta[q,p]), log=T) -
            dnorm(x=B_update[q,p], 0, sqrt(gamma_beta[q,p]*nu_beta[q,p]), log=T)
          # accept or reject
          test <- log(runif(1))
          if (test<ratio){ B_update <- B_new }
        }
      }
    }
  }

  return(B_update)
}



update_gamma_beta <- function(B, nu_beta, rho_beta, prior_lst){

  Q <- ncol(B)

  gamma_beta_update <- array(NA, dim=c(Q, Q))
  nu_0 <- prior_lst$nu_0

  for (q in 1:Q){
    for (p in 1:Q){

      # Assign 0 to the diagonal elements, so rho_beta's update will not be affected
      if(q == p){

        gamma_beta_update[q,p] <- 0

      } else{

        # prob = pr(gamma_beta_lqp=1|~) / pr(gamma_beta_lqp=nu_0|~)
        prob <- ((sqrt(nu_0)*rho_beta)/(1-rho_beta)) * exp((1-nu_0)*B[q,p]^2/(2*nu_0*nu_beta[q,p]))
        if (prob == Inf){
          gamma_beta_update[q,p] <- 1
        }else{
          gamma_beta_update[q,p] <- sample(c(1,nu_0), 1, prob=c(prob, 1))
        }
      }

    }
  }

  return(gamma_beta_update)
}



update_nu_beta <- function(B, gamma_beta, prior_lst){

  Q <- ncol(B)

  nu_beta_update <- matrix(NA, Q, Q)
  a_nu <- prior_lst$a_nu
  b_nu <- prior_lst$b_nu

  for (q in 1:Q){
    for (p in 1:Q){

      # No need to update for the diagonal elements
      if(q == p){

        nu_beta_update[q,p] <- 0

      } else{

        a_nu_star <- a_nu + 1/2
        b_nu_star <- b_nu + B[q,p]^2/(2*gamma_beta[q,p])
        nu_beta_update[q,p] <- rinvgamma(1, a_nu_star, b_nu_star)

      }

    }
  }

  return(nu_beta_update)
}



update_rho_beta <- function(gamma_beta, prior_lst){

  a_rho <- prior_lst$a_rho
  b_rho <- prior_lst$b_rho
  nu_0 <- prior_lst$nu_0

  a_rho_star <- a_rho + sum(gamma_beta == 1)
  b_rho_star <- b_rho + sum(gamma_beta == nu_0)
  rho_beta_update <- rbeta(1, a_rho_star, b_rho_star)

  return(rho_beta_update)
}
