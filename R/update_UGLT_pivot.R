update_pivot <- function(sparsity_matrix, pivot, a1, a2, P_star, C, B, A, mu, tau, L_0, prior_lst, mcmc_setup_lst, data){

  pshift <- mcmc_setup_lst$pshift
  pswitch <- mcmc_setup_lst$pswitch

  pivot_update <- pivot
  sparsity_matrix_update <- sparsity_matrix

  # update columns in a random order
  for(p in sample(1:P_star)){
    d_p <- colSums(sparsity_matrix_update)

    move <- sample(c("shift", "switch", "add/del"), size = 1, prob = c(pshift, pswitch, 1 - pshift - pswitch))

    # when changing a pivot, the sparsity matrix should also be updated
    if(move == "shift"){
      res_lst <- update_pivot_shift(p, d_p, sparsity_matrix_update, pivot_update, a1, a2, C, B, A, mu, tau, L_0, prior_lst, data)
    } else if(move == "switch"){
      res_lst <- update_pivot_switch(p, d_p, sparsity_matrix_update, pivot_update, a1, a2, P_star, C, B, A, mu, tau, L_0, prior_lst, data)
    } else{
      res_lst <- update_pivot_add_del(p, d_p, sparsity_matrix_update, pivot_update, a1, a2, C, B, A, mu, tau, L_0, prior_lst, mcmc_setup_lst, data)
    }

    pivot_update <- res_lst$pivot_update
    sparsity_matrix_update <- res_lst$sparsity_matrix_update
  }

  return(list(sparsity_matrix = sparsity_matrix_update, pivot = pivot_update))
}



update_pivot_shift <- function(p, d_p, sparsity_matrix, pivot, a1, a2, C, B, A, mu, tau, L_0, prior_lst, data){

  Y <- data$Y
  X <- data$X
  Q <- ncol(Y)
  H <- prior_lst$H
  a_sigma <- prior_lst$a_sigma
  b_sigma <- prior_lst$b_sigma

  pivot_update <- pivot
  sparsity_matrix_update <- sparsity_matrix

  pivot_star <- which(sparsity_matrix[,p] != 0)[2]

  if(is.na(pivot_star)){
    # If we are having a spurious column
    possible_pivot <- (1:Q)[!(1:Q) %in% pivot[-p]]
  } else{
    # If there is no room for shift, just return
    if (pivot_star <= 2){
      return(list(sparsity_matrix_update = sparsity_matrix_update, pivot_update = pivot_update))
    }
    possible_pivot <- (1:(pivot_star-1))[!(1:(pivot_star-1)) %in% pivot[-p]]
  }

  if (length(possible_pivot) < 2){
    # if only one possible pivot is available no shift move is performed since accept prob is one (the new proposed pivot would be it self)
    return(list(sparsity_matrix_update = sparsity_matrix_update, pivot_update = pivot_update))
  }

  # If shift is possible (more than one possible pivot)
  pivot_new <- sample(possible_pivot, 1)
  if(pivot_new == pivot[p]){
    # the old pivot is selected, directly return
    return(list(sparsity_matrix_update = sparsity_matrix_update, pivot_update = pivot_update))
  }

  # Compute the accept probability
  # Log prior
  log_prior_old <- lbeta(a1*a2/H + d_p[p] - 1, a2 + Q - pivot[p] - d_p[p] + 1)
  log_prior_new <- lbeta(a1*a2/H + d_p[p] - 1, a2 + Q - pivot_new - d_p[p] + 1)

  # log lik
  loglik_old <- compute_marginal_loglik_rcpp(sparsity_matrix, pivot[p], C, tau, L_0, Y, B, X, A, mu, a_sigma, b_sigma) +
    compute_marginal_loglik_rcpp(sparsity_matrix, pivot_new, C, tau, L_0, Y, B, X, A, mu, a_sigma, b_sigma)

  sparsity_matrix_new <- sparsity_matrix
  sparsity_matrix_new[pivot[p],p] <- 0
  sparsity_matrix_new[pivot_new,p] <- 1
  loglik_new <- compute_marginal_loglik_rcpp(sparsity_matrix_new, pivot[p], C, tau, L_0, Y, B, X, A, mu, a_sigma, b_sigma) +
    compute_marginal_loglik_rcpp(sparsity_matrix_new, pivot_new, C, tau, L_0, Y, B, X, A, mu, a_sigma, b_sigma)

  # MH step
  log_ratio <- log_prior_new + loglik_new - log_prior_old - loglik_old
  if(log(runif(1)) < log_ratio){
    pivot_update[p] <- pivot_new
    sparsity_matrix_update <- sparsity_matrix_new
  }

  return(list(sparsity_matrix_update = sparsity_matrix_update, pivot_update = pivot_update))
}



update_pivot_switch <- function(p, d_p, sparsity_matrix, pivot, a1, a2, P_star, C, B, A, mu, tau, L_0, prior_lst, data){

  Y <- data$Y
  X <- data$X
  Q <- ncol(Y)
  H <- prior_lst$H
  a_sigma <- prior_lst$a_sigma
  b_sigma <- prior_lst$b_sigma

  pivot_update <- pivot
  sparsity_matrix_update <- sparsity_matrix

  # check if switch is applicable
  if(P_star == 1){
    return(list(sparsity_matrix_update = sparsity_matrix_update,
                pivot_update = pivot_update))
  }

  # If switch is possible
  # sample a column to switch with
  possible_columns <- (1:P_star)[!(1:P_star) %in% p]
  if(length(possible_columns) == 1){
    column_2_switch <- possible_columns
  } else{
    column_2_switch <- sample(possible_columns, 1)
  }

  # Construct the new sparsity matrix and pivot
  tmp_lst <- switch_2_columns_pivots(sparsity_matrix, pivot, p, column_2_switch)
  sparsity_matrix_new <- tmp_lst$sparsity_matrix_new
  pivot_new <- tmp_lst$pivot_new
  diff_row_indices <- tmp_lst$diff_row_indices
  d_p_new <- colSums(sparsity_matrix_new)

  # Compute the accept probability
  # Log prior
  log_prior_old <- lbeta(a1*a2/H + d_p[p] - 1, a2 + Q - pivot[p] - d_p[p] + 1) +
    lbeta(a1*a2/H + d_p[column_2_switch] - 1, a2 + Q - pivot[column_2_switch] - d_p[column_2_switch] + 1)
  log_prior_new <- lbeta(a1*a2/H + d_p_new[p] - 1, a2 + Q - pivot_new[p] - d_p_new[p] + 1) +
    lbeta(a1*a2/H + d_p_new[column_2_switch] - 1, a2 + Q - pivot_new[column_2_switch] - d_p_new[column_2_switch] + 1)

  # log lik
  loglik_old <- 0
  loglik_new <- 0
  for(diff_row_index in diff_row_indices){
    loglik_old <- loglik_old + compute_marginal_loglik_rcpp(sparsity_matrix, diff_row_index, C, tau, L_0, Y, B, X, A, mu, a_sigma, b_sigma)
    loglik_new <- loglik_new + compute_marginal_loglik_rcpp(sparsity_matrix_new, diff_row_index, C, tau, L_0, Y, B, X, A, mu, a_sigma, b_sigma)
  }

  # MH step
  log_ratio <- log_prior_new + loglik_new - log_prior_old - loglik_old
  if(log(runif(1)) < log_ratio){
    pivot_update <- pivot_new
    sparsity_matrix_update <- sparsity_matrix_new
  }

  return(list(sparsity_matrix_update = sparsity_matrix_update, pivot_update = pivot_update))
}



update_pivot_add_del <- function(p, d_p, sparsity_matrix, pivot, a1, a2, C, B, A, mu, tau, L_0, prior_lst, mcmc_setup_lst, data){

  Y <- data$Y
  X <- data$X
  Q <- ncol(Y)
  H <- prior_lst$H
  pa <- mcmc_setup_lst$pa
  a_sigma <- prior_lst$a_sigma
  b_sigma <- prior_lst$b_sigma

  pivot_update <- pivot
  sparsity_matrix_update <- sparsity_matrix

  # Check if add/delete a pivot is applicable
  tmp_lst <- compute_padd(p, pivot, sparsity_matrix, pa)
  padd <- tmp_lst$padd
  addable_pivot <- tmp_lst$addable_pivot

  # Decide add/delete/stay
  if(is.na(padd)){
    return(list(sparsity_matrix_update = sparsity_matrix_update,
                pivot_update = pivot_update))
  } else{
    add <- ifelse(rbinom(1, 1, padd) == 1, TRUE, FALSE)
  }

  # MH step if applicable
  if(add){
    if(length(addable_pivot) == 1){
      pivot_2_add <- addable_pivot
    } else{
      pivot_2_add <- sample(addable_pivot, 1)
    }

    # log prior
    log_prior_0 <- lbeta(a1*a2/H + d_p[p] - 1, a2 + Q - pivot[p] - d_p[p] + 1)
    log_prior_1 <- lbeta(a1*a2/H + d_p[p], a2 + Q - pivot_2_add - d_p[p])

    # log lik
    sparsity_matrix_0 <- sparsity_matrix
    sparsity_matrix_0[pivot_2_add,p] <- 0
    loglik_0 <- compute_marginal_loglik_rcpp(sparsity_matrix_0, pivot_2_add, C, tau, L_0, Y, B, X, A, mu, a_sigma, b_sigma)

    sparsity_matrix_1 <- sparsity_matrix
    sparsity_matrix_1[pivot_2_add,p] <- 1
    loglik_1 <- compute_marginal_loglik_rcpp(sparsity_matrix_1, pivot_2_add, C, tau, L_0, Y, B, X, A, mu, a_sigma, b_sigma)

    # log proposal
    ## pivot_tmp created to compute padd for the new pivot proposal prob
    pivot_tmp <- pivot
    pivot_tmp[p] <- pivot_2_add
    padd_new <- compute_padd(p, pivot_tmp, sparsity_matrix_1, pa)$padd
    log_proposal_new_2_old <- log(1 - padd_new)
    log_proposal_old_2_new <- log(padd) - log(length(addable_pivot))

    # MH step
    log_ratio <- log_prior_1 + loglik_1 + log_proposal_new_2_old - log_prior_0 - loglik_0 - log_proposal_old_2_new

    if(log(runif(1)) < log_ratio){
      pivot_update[p] <- pivot_2_add
      sparsity_matrix_update <- sparsity_matrix_1
    }

  } else{
    pivot_after_del <- which(sparsity_matrix[,p] != 0)[2]

    # log prior
    log_prior_0 <- lbeta(a1*a2/H + d_p[p] - 2, a2 + Q - pivot_after_del - d_p[p] + 2)
    log_prior_1 <- lbeta(a1*a2/H + d_p[p] - 1, a2 + Q - pivot[p] - d_p[p] + 1)

    # log lik
    # 0 is to set the pivot to be deleted to 0
    sparsity_matrix_0 <- sparsity_matrix
    sparsity_matrix_0[pivot[p],p] <- 0
    loglik_0 <- compute_marginal_loglik_rcpp(sparsity_matrix_0, pivot[p], C, tau, L_0, Y, B, X, A, mu, a_sigma, b_sigma)

    sparsity_matrix_1 <- sparsity_matrix
    sparsity_matrix_1[pivot[p],p] <- 1
    loglik_1 <- compute_marginal_loglik_rcpp(sparsity_matrix_1, pivot[p], C, tau, L_0, Y, B, X, A, mu, a_sigma, b_sigma)

    # log proposal
    ## pivot_tmp created to compute padd for the new pivot proposal prob
    pivot_tmp <- pivot
    pivot_tmp[p] <- pivot_after_del
    tmp_lst <- compute_padd(p, pivot_tmp, sparsity_matrix_0, pa)
    log_proposal_new_2_old <- log(tmp_lst$padd) - log(length(tmp_lst$addable_pivot))
    log_proposal_old_2_new <- log(1 - padd)

    # MH step
    log_ratio <- log_prior_0 + loglik_0 + log_proposal_new_2_old - log_prior_1 - loglik_1 - log_proposal_old_2_new

    if(log(runif(1)) < log_ratio){
      pivot_update[p] <- pivot_after_del
      sparsity_matrix_update <- sparsity_matrix_0
    }
  }

  return(list(sparsity_matrix_update = sparsity_matrix_update,
              pivot_update = pivot_update))
}



