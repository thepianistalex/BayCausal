#' Title
#'
#' @param seed random seed
#' @param P_star  number of latent factors
#' @param data data list
#' @param prior_lst a list containing prior hyperparameters
#'
#' @returns a list containing initial values
#' @export
set_init_default <- function(seed, P_star, data, prior_lst){
  
  set.seed(seed)
  
  Y <- data$Y
  X <- data$X
  n <- nrow(Y)
  Q <- ncol(Y)
  S <- ncol(X)
  nu_0 <- prior_lst$nu_0
  
  A_init <- matrix(NA, Q, S)
  gamma_alpha_init <- matrix(NA, Q, S)
  nu_alpha_init <- matrix(NA, Q, S)
  rho_alpha_init <- 0.5
  
  B_init <- matrix(NA, Q, Q)
  gamma_beta_init <- matrix(NA, Q, Q)
  nu_beta_init <- matrix(NA, Q, Q)
  rho_beta_init <- 0.5
  
  tau_init <- matrix(NA, n, Q)
  mu_init <- rep(0,Q)
  sigma2_init <- rep(1,Q)
  
  L_init <- matrix(NA, Q, P_star)
  C_init <- matrix(rnorm(n*P_star,0,1), n, P_star)
  sparsity_matrix_init <- matrix(1, Q, P_star)
  pivot_init <- 1:P_star
  
  a_init <- c(1,2)
  
  kappa_init <- 1
  
  # initialize alpha
  for (q in 1:Q){
    for (s in 1:S){
      gamma_alpha_init[q,s] <- sample(c(nu_0,1), 1, prob=c(0.5,0.5))
      nu_alpha_init[q,s] <- 1
      A_init[q,s] <- rnorm(1, 0, sqrt(gamma_alpha_init[q,s]*nu_alpha_init[q,s]))
    }
  }
  
  # initialize beta
  for (q in 1:Q){
    for (p in 1:Q){
      gamma_beta_init[q,p] <- sample(c(nu_0,1), 1, prob=c(0.5,0.5))
      nu_beta_init[q,p] <- 1
      B_init[q,p] <- rnorm(1, 0, sqrt(gamma_beta_init[q,p]*nu_beta_init[q,p]))
    }
  }
  
  # initialize sparsity matrix
  if(P_star > 1){
    for(p in 2:P_star){
      for(q in 1:(p-1)){
        sparsity_matrix_init[q,p] <- 0
      }
    }
  }
  
  # initialize L
  L_init <- sparsity_matrix_init
  
  # initialize zeta
  zeta_init <- update_zeta(a_init[1], a_init[2], sparsity_matrix_init, pivot_init, P_star, prior_lst)
  
  # initialize tau
  for (i in 1:n){
    for (q in 1:Q){
      Y_tilde_iq <- Y[i,q] - t(Y[i,])%*%B_init[q,] - t(X[i,])%*%A_init[q,] - t(C_init[i,])%*%L_init[q,]
      tau_init[i,q] <- rinvgaussian(1, mu=sqrt(sigma2_init[q])/(2*abs(Y_tilde_iq)), lambda=1/4)
    }
  }
  
  init_list <- list(A=A_init, gamma_alpha=gamma_alpha_init, nu_alpha=nu_alpha_init, rho_alpha=rho_alpha_init,
                    B=B_init, gamma_beta=gamma_beta_init, nu_beta=nu_beta_init, rho_beta=rho_beta_init,
                    C=C_init, L=L_init, a=a_init, tau=tau_init, mu=mu_init, sigma2=sigma2_init, P_star=P_star,
                    sparsity_matrix=sparsity_matrix_init, pivot=pivot_init, zeta=zeta_init, kappa=kappa_init)
  
  return(init_list)    
}



#' Title
#'
#' @returns a list containing chain setup
#' @export
set_mh_default <- function(){
  
  mh_setup_lst <- list()
  
  mh_setup_lst$B_step <- 0.02
  mh_setup_lst$a1_step <- 0.1
  mh_setup_lst$a2_step <- 0.1
  mh_setup_lst$pshift <- 1/3
  mh_setup_lst$pswitch <- 1/3
  mh_setup_lst$pa <- 0.5
  mh_setup_lst$ps <- 0.5
  
  return(mh_setup_lst)
}


#' Title
#'
#' @returns a list containing chain setup
#' @export
set_prior_default <- function(){
  
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
  prior_lst$b_a1 <- prior_lst$a_a1*(prior_lst$H-Eq)/(prior_lst$H*Eq)
  prior_lst$a_a2 <- prior_lst$a_a1
  prior_lst$b_a2 <- prior_lst$a_a1
  prior_lst$a_kappa <- 1
  prior_lst$b_kappa <- 1
  
  return(prior_lst)
}
