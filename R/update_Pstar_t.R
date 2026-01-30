update_P_star_t <- function(mu, A, B, L, C, sigma2, a1, a2, sparsity_matrix, pivot, P_star, kappa, mh_setup_lst, prior_lst, data, nu_t){

  # check if split/merge is applicable
  psplit <- compute_psplit(pivot, sparsity_matrix, mh_setup_lst$ps, prior_lst$H)

  if(is.na(psplit)){
    # no move is possible, return the current values
    res_lst <- list(mcmc_P_star = P_star, mcmc_L = L, mcmc_C = C, mcmc_sigma2 = sigma2,
                    mcmc_sparsity_matrix = sparsity_matrix, mcmc_pivot = pivot)
    return(res_lst)
  } else{
    move <- ifelse(rbinom(1, 1, psplit) == 1, "split", "merge")
  }

  if(move == "split"){
    res_lst <- update_P_star_split_t(mu, A, B, L, C, sigma2, a1, a2, sparsity_matrix, pivot, P_star, kappa, mh_setup_lst, prior_lst, data, nu_t)
  } else{
    res_lst <- update_P_star_merge_t(mu, A, B, L, C, sigma2, a1, a2, sparsity_matrix, pivot, P_star, kappa, mh_setup_lst, prior_lst, data, nu_t)
  }

  return(res_lst)
}



update_P_star_split_t <- function(mu, A, B, L, C, sigma2, a1, a2, sparsity_matrix, pivot, P_star, kappa, mh_setup_lst, prior_lst, data, nu_t){

  Y <- data$Y
  X <- data$X
  n <- nrow(Y)
  Q <- ncol(Y)
  H <- prior_lst$H
  a_sigma <- prior_lst$a_sigma
  b_sigma <- prior_lst$b_sigma
  ps <- mh_setup_lst$ps

  # Sample a pivot for the new spurious column
  available_pivots <- (1:Q)[!(1:Q) %in% pivot]
  if(length(available_pivots) == 1){
    pivot_propose <- available_pivots
  } else{
    pivot_propose <- sample(available_pivots, 1)
  }

  # Sample U for dimension matching
  U <- runif(1,-1,1)

  # Specify the new sparsity matrix and pivot
  pivot_new <- c(pivot, pivot_propose)
  column_new <- rep(0, Q)
  column_new[pivot_propose] <- 1
  sparsity_matrix_new <- cbind(sparsity_matrix, column_new)

  # New sigma2
  # Variance conservation for Student-t: Kt = nu/(nu-2)
  # Kt * sigma2_old = Kt * sigma2_new + L^2
  # sigma2_new = (1-U^2) * sigma2_old (unchanged proposal structure)
  # => L^2 = Kt * U^2 * sigma2_old
  sigma2_propose <- (1-U^2)*sigma2[pivot_propose]
  sigma2_new <- sigma2
  sigma2_new[pivot_propose] <- sigma2_propose

  # New factor
  C_propose <- rnorm(n, 0, 1)
  C_new <- cbind(C, C_propose)

  # New loading
  Kt <- nu_t / (nu_t - 2.0)
  l_propose <- U*sqrt(Kt*sigma2[pivot_propose])
  L_propose <- rep(0, Q)
  L_propose[pivot_propose] <- l_propose
  L_new <- cbind(L, L_propose)

  # log(Acceptance ratio) = log(Prior ratio) + log(Likelihood ratio) + log(Proposal ratio) + log(Jacobian)
  # log prior ratio, if using prior to propose factor, then cancelled out with the proposal ratio's factor component
  prior_ratio_sigma2 <- dgamma(1/sigma2_propose, a_sigma, b_sigma, log = TRUE) - dgamma(1/sigma2[pivot_propose], a_sigma, b_sigma, log = TRUE)
  prior_ratio_L <- sum(ifelse(L_new[pivot_propose,] == 0, 0, dnorm(L_new[pivot_propose,], 0, sqrt(kappa*sigma2_propose), log = TRUE))) -
    sum(ifelse(L[pivot_propose,] == 0, 0, dnorm(L[pivot_propose,], 0, sqrt(kappa*sigma2[pivot_propose]), log = TRUE)))
  r <- sum(colSums(sparsity_matrix) > 1)
  r_sp <- sum(colSums(sparsity_matrix) == 1)
  aH <- a1 * a2 / H
  prior_ratio_delta <- log(aH) + log(Q-r-r_sp) - log(a2-1+Q-r-r_sp)
  log_prior_ratio <- prior_ratio_sigma2 + prior_ratio_L + prior_ratio_delta

  # log likelihood ratio - USING T LIKELIHOOD
  loglik_old <- compute_loglik_rjmcmc_t_rcpp(pivot_propose, B, A, L, C, mu, sigma2, Y, X, nu_t)
  loglik_new <- compute_loglik_rjmcmc_t_rcpp(pivot_propose, B, A, L_new, C_new, mu, sigma2_new, Y, X, nu_t)
  log_lik_ratio <- loglik_new - loglik_old

  # log proposal ratio
  q_split <- ps / (H - r - r_sp)
  q_merge <- ps / (r_sp + 1)
  log_proposal_ratio <- log(2) + log(q_merge) - log(q_split)

  # Jacobian
  # J = sqrt(Kt * sigma2_old)
  log_jacobian <- log(sqrt(Kt)) + log(sqrt(sigma2[pivot_propose]))

  # Acceptance ratio
  log_accept_prob <- log_prior_ratio + log_lik_ratio + log_proposal_ratio + log_jacobian

  if(log(runif(1)) < log_accept_prob){
    return(list(mcmc_P_star = P_star + 1, mcmc_L = L_new, mcmc_C = C_new, mcmc_sigma2 = sigma2_new,
                mcmc_sparsity_matrix = sparsity_matrix_new, mcmc_pivot = pivot_new))
  } else{
    return(list(mcmc_P_star = P_star, mcmc_L = L, mcmc_C = C, mcmc_sigma2 = sigma2,
                mcmc_sparsity_matrix = sparsity_matrix, mcmc_pivot = pivot))
  }
}



## To be updated
update_P_star_merge_t <- function(mu, A, B, L, C, sigma2, a1, a2, sparsity_matrix, pivot, P_star, kappa, mh_setup_lst, prior_lst, data, nu_t){

  Y <- data$Y
  X <- data$X
  n <- nrow(Y)
  Q <- ncol(Y)
  H <- prior_lst$H
  a_sigma <- prior_lst$a_sigma
  b_sigma <- prior_lst$b_sigma
  ps <- mh_setup_lst$ps

  # Sample a spurious column to merge
  sp_columns <- which(colSums(sparsity_matrix) == 1)
  if(length(sp_columns) == 1){
    column_merge <- sp_columns
  } else{
    column_merge <- sample(sp_columns, 1)
  }
  pivot_merge <- pivot[column_merge]

  # Specify the new sparsity matrix and pivot
  pivot_new <- pivot[-column_merge]
  sparsity_matrix_new <- sparsity_matrix[,-column_merge, drop = FALSE]

  # New sigma2
  Kt <- nu_t / (nu_t - 2.0)
  sigma2_new <- sigma2
  sigma2_propose <- sigma2[pivot_merge] + (L[pivot_merge,column_merge]^2)/Kt
  sigma2_new[pivot_merge] <- sigma2_propose

  # New factor
  C_new <- C[,-column_merge, drop = FALSE]

  # New loading
  L_new <- L[,-column_merge, drop = FALSE]

  # log(Acceptance ratio) = log(Prior ratio) + log(Likelihood ratio) + log(Proposal ratio) + log(Jacobian)
  # log prior ratio, if using prior to propose factor, then cancelled out with the proposal ratio's factor component
  prior_ratio_sigma2 <- dgamma(1/sigma2_propose, a_sigma, b_sigma, log = TRUE) - dgamma(1/sigma2[pivot_merge], a_sigma, b_sigma, log = TRUE)
  prior_ratio_L <- sum(ifelse(L_new[pivot_merge,] == 0, 0, dnorm(L_new[pivot_merge,], 0, sqrt(kappa*sigma2_propose), log = TRUE))) -
    sum(ifelse(L[pivot_merge,] == 0, 0, dnorm(L[pivot_merge,], 0, sqrt(kappa*sigma2[pivot_merge]), log = TRUE)))
  r <- sum(colSums(sparsity_matrix) > 1)
  r_sp <- sum(colSums(sparsity_matrix) == 1)
  aH <- a1 * a2 / H
  prior_ratio_delta <- log(a2+Q-r-r_sp) - log(aH) - log(Q-r-r_sp+1)
  log_prior_ratio <- prior_ratio_sigma2 + prior_ratio_L + prior_ratio_delta

  # log likelihood ratio - USING T LIKELIHOOD
  loglik_old <- compute_loglik_rjmcmc_t_rcpp(pivot_merge, B, A, L, C, mu, sigma2, Y, X, nu_t)
  loglik_new <- compute_loglik_rjmcmc_t_rcpp(pivot_merge, B, A, L_new, C_new, mu, sigma2_new, Y, X, nu_t)
  log_lik_ratio <- loglik_new - loglik_old

  # log proposal ratio
  q_split <- ps / (H - r - r_sp + 1)
  q_merge <- ps / (r_sp)
  log_proposal_ratio <- log(q_split) - log(2) - log(q_merge)

  # Jacobian
  log_jacobian <- -(log(sqrt(Kt)) + log(sqrt(sigma2_new[pivot_merge])))

  # Acceptance ratio
  log_accept_prob <- log_prior_ratio + log_lik_ratio + log_proposal_ratio + log_jacobian

  if(log(runif(1)) < log_accept_prob){
    return(list(mcmc_P_star = P_star - 1, mcmc_L = L_new, mcmc_C = C_new, mcmc_sigma2 = sigma2_new,
                mcmc_sparsity_matrix = sparsity_matrix_new, mcmc_pivot = pivot_new))
  } else{
    return(list(mcmc_P_star = P_star, mcmc_L = L, mcmc_C = C, mcmc_sigma2 = sigma2,
                mcmc_sparsity_matrix = sparsity_matrix, mcmc_pivot = pivot))
  }
}
