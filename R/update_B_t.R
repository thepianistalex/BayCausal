#' Update B using Metropolis-Hastings with Student-t errors
#'
#' @param B current B matrix
#' @param gamma_beta inclusion indicators
#' @param nu_beta variances
#' @param A coefficient matrix for X
#' @param mu mean vector
#' @param L factor loadings
#' @param C factor scores
#' @param sigma2 noise variances (Normal scale parameter squared)
#' @param mh_step_size step size for MH
#' @param data data object
#' @param nu_t degrees of freedom for Student-t
#'
#' @export
update_B_t <- function(B, gamma_beta, nu_beta, A, mu, L, C, sigma2, mh_step_size, data, nu_t){

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
          # calculate the acceptance ratio using Student-t likelihood
          ratio <- compute_loglik_B_t_rcpp(q, B_new, A, L, C, mu, sigma2, Y, X, nu_t) -
            compute_loglik_B_t_rcpp(q, B_update, A, L, C, mu, sigma2, Y, X, nu_t) +
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
