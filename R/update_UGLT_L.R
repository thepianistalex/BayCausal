update_L_uglt_shrink <- function(mu, A, B, C, tau, sigma2, L_0, sparsity_matrix, P_star, data){

  Q <- ncol(B)
  Y <- data$Y
  X <- data$X

  L_update <- sparsity_matrix

  for(q in 1:Q){

    non_zero_p <- which(sparsity_matrix[q, ] == 1)

    if(length(non_zero_p) == 0){
      next
    }

    Y_q_tilda <- Y[, q] - mu[q] - Y %*% B[q,] - X %*% A[q,]
    C_q <- as.matrix(C[, non_zero_p], ncol = length(non_zero_p))

    if (length(non_zero_p) == 1) {
      L_0_q <- matrix(L_0[q, non_zero_p], nrow = 1, ncol = 1)
    } else {
      L_0_q <- diag(L_0[q, non_zero_p])
    }

    V_q_n <- chol2inv(chol(t(C_q)%*%(tau[,q]*C_q) + solve(L_0_q)))
    b_q <- t(C_q)%*%(tau[,q]*Y_q_tilda)

    L_update[q, non_zero_p] <- as.vector(rmvn_rcpp(1, V_q_n%*%b_q, sigma2[q]*V_q_n))

  }

  return(L_update)
}
