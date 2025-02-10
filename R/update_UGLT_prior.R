update_a1 <- function(a1_current, a2, sparsity_matrix, pivot, mh_step_size, prior_lst) {

  a1_proposed <- rnorm(1, a1_current, mh_step_size)
  if(a1_proposed < 0){
    return(a1_current)
  }
  log_acceptance_ratio <- log_posterior_a1(a1_proposed, a2, sparsity_matrix, pivot, prior_lst) -
    log_posterior_a1(a1_current, a2, sparsity_matrix, pivot, prior_lst)

  if (log(runif(1)) < log_acceptance_ratio) {
    return(a1_proposed)
  } else {
    return(a1_current)
  }
}



update_a2 <- function(a1, a2_current, sparsity_matrix, pivot, mh_step_size, prior_lst) {

  a2_proposed <- rnorm(1, a2_current, mh_step_size)
  if(a2_proposed < 0){
    return(a2_current)
  }

  log_acceptance_ratio <- log_posterior_a2(a1, a2_proposed, sparsity_matrix, pivot, prior_lst) -
    log_posterior_a2(a1, a2_current, sparsity_matrix, pivot, prior_lst)

  if (log(runif(1)) < log_acceptance_ratio) {
    return(a2_proposed)
  } else {
    return(a2_current)
  }
}



log_posterior_a1 <- function(a1, a2, sparsity_matrix, pivot, prior_lst) {

  a_a1 <- prior_lst$a_a1
  b_a1 <- prior_lst$b_a1
  H <- prior_lst$H
  Q <- nrow(sparsity_matrix)

  # Compute values from sparsity matrix
  d_p <- colSums(sparsity_matrix)
  active_columns <- which(d_p > 1)
  sp_columns <- which(d_p == 1)
  r <- length(active_columns)
  r_sp <- length(sp_columns)

  # Compute prior and common factors
  log_prior <- dgamma(a1, shape = a_a1, rate = b_a1, log = TRUE) -
    H * lbeta(a1*a2/H, a2)

  # Compute 0 columns
  log_zero <- (H-r-r_sp) * lbeta(a1*a2/H, a2 + Q - r - r_sp)

  # Compute active term
  log_active <- 0
  for(active_column in active_columns) {
    log_active <- log_active + lbeta(a1*a2/H + d_p[active_column] - 1,
                                     a2 + Q - pivot[active_column] - d_p[active_column] + 1)
  }

  # Compute spurious term
  log_spurious <- 0
  for(p in 1:r_sp) {
    log_spurious <- log_spurious + lbeta(a1*a2/H + 1, a2 + Q - r - p)
  }

  # Putting it all together
  log_post <- log_prior + log_zero + log_active + log_spurious

  return(log_post)
}



log_posterior_a2 <- function(a1, a2, sparsity_matrix, pivot, prior_lst) {

  a_a2 <- prior_lst$a_a2
  b_a2 <- prior_lst$b_a2
  H <- prior_lst$H
  Q <- nrow(sparsity_matrix)

  # Compute values from sparsity matrix
  d_p <- colSums(sparsity_matrix)
  active_columns <- which(d_p > 1)
  sp_columns <- which(d_p == 1)
  r <- length(active_columns)
  r_sp <- length(sp_columns)

  # Compute prior and common factors
  log_prior <- dgamma(a2, shape = a_a2, rate = b_a2, log = TRUE) -
    H * lbeta(a1*a2/H, a2)

  # Compute 0 columns
  log_zero <- (H-r-r_sp) * lbeta(a1*a2/H, a2 + Q - r - r_sp)

  # Compute active term
  log_active <- 0
  for(active_column in active_columns) {
    log_active <- log_active + lbeta(a1*a2/H + d_p[active_column] - 1,
                                     a2 + Q - pivot[active_column] - d_p[active_column] + 1)
  }

  # Compute spurious term
  log_spurious <- 0
  for(p in 1:r_sp) {
    log_spurious <- log_spurious + lbeta(a1*a2/H + 1, a2 + Q - r - p)
  }

  # Putting it all together
  log_post <- log_prior + log_zero + log_active + log_spurious

  return(log_post)
}



update_zeta <- function(a1, a2, sparsity_matrix, pivot, P_star, prior_lst){

  zeta_update <- rep(NA, P_star)
  H <- prior_lst$H
  Q <- nrow(sparsity_matrix)

  # Compute values from sparsity matrix
  d_p <- colSums(sparsity_matrix)
  active_columns <- which(d_p > 1)
  sp_columns <- which(d_p == 1)

  # Actually they are the same
  for(active_column in active_columns) {
    zeta_update[active_column] <- rbeta(1, a1*a2/H + d_p[active_column] - 1,
                                        a2 + Q - pivot[active_column] - d_p[active_column] + 1)
  }

  for(sp_column in sp_columns) {
    zeta_update[sp_column] <- rbeta(1, a1*a2/H,
                                    a2 + Q - pivot[sp_column])
  }

  return(zeta_update)
}



update_kappa <- function(sparsity_matrix, L, sigma2, prior_lst){

  a_kappa <- prior_lst$a_kappa
  b_kappa <- prior_lst$b_kappa

  sum_sq <- sum(sparsity_matrix * L^2 / sigma2)

  a_new <- a_kappa + sum(sparsity_matrix)/2
  b_new <- b_kappa + sum_sq/2

  kappa_update <- 1/rgamma(1, a_new, b_new)
  return(kappa_update)
}
