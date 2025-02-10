update_sparsity_indicator <- function(sparsity_matrix, pivot, zeta, C, B, A, mu, tau, L_0, P_star, prior_lst, data){

  Y <- data$Y
  X <- data$X
  Q <- ncol(Y)
  a_sigma <- prior_lst$a_sigma
  b_sigma <- prior_lst$b_sigma

  sparsity_matrix_update <- sparsity_matrix

  for(p in 1:P_star){

    if(pivot[p] == Q){
      next
    }

    for(q in (pivot[p]+1) : Q){
      # log prior
      log_prior_0 <- log(1 - zeta[p])
      log_prior_1 <- log(zeta[p])

      # loglik
      sparsity_matrix_0 <- sparsity_matrix_update
      sparsity_matrix_0[q,p] <- 0
      loglik_0 <- compute_marginal_loglik_rcpp(sparsity_matrix_0, q, C, tau, L_0, Y, B, X, A, mu, a_sigma, b_sigma)

      sparsity_matrix_1 <- sparsity_matrix_update
      sparsity_matrix_1[q,p] <- 1
      loglik_1 <- compute_marginal_loglik_rcpp(sparsity_matrix_1, q, C, tau, L_0, Y, B, X, A, mu, a_sigma, b_sigma)

      # Compute log odds ratio
      log_O <- log_prior_1 + loglik_1 - log_prior_0 - loglik_0
      p1 <- 1 / (1 + exp(-log_O))
      sparsity_matrix_update[q,p] <- rbinom(1, 1, p1)
    }
  }

  return(sparsity_matrix_update)
}
