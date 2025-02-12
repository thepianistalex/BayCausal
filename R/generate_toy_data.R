#' Title
#'
#' @returns toy data list
#' @export
generate_toy_data <- function(){
  
  set.seed(0)
  Q <- 3
  S <- 2
  P <- 1
  n <- 3000
  
  mu <- runif(Q, -1, 1)
  
  B <- matrix(c(
    0, 0.5, 0,
    0, 0, 0.7,
    0, 0, 0
  ), nrow = Q, byrow = TRUE)
  
  A <- matrix(c(runif(Q*S, -1, 1)), Q, S)
  
  L <- matrix(c(
    0,
    0.5,
    0.3), ncol = 1)
  sparsity_matrix <- ifelse(L != 0, 1, 0)
  pivot <- numeric(ncol(L))
  for (p in 1:ncol(L)) {
    pivot[p] <- ifelse(any(L[, p] != 0), which(L[, p] != 0)[1], NA)
  }
  
  X <- matrix(rnorm(n*S), n, S)
  C <- matrix(rnorm(n*P), n, P)
  
  sigma_e <- rep(sqrt(1/8), Q)
  E <- matrix(NA, n, Q)
  for(q in 1:Q){
    E[,q] <- LaplacesDemon::rlaplace(n, 0, 2*sigma_e[q])
  }
  
  Y <- matrix(NA, n, Q)
  mix_mat <- solve(diag(Q)-B)
  for(i in 1:n){
    Y[i,] <- mix_mat %*% (mu + A%*%X[i,] + L%*%C[i,] + E[i,])
  }
  
  data <- list(Y = Y, X = X, mu = mu, B = B, A = A, L = L, C = C, sigma_e = sigma_e, Q = Q, S = S, P = P, n = n)

  return(data)
}