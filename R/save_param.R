save_check <- function(it, burn, thin){
  if((it > burn) & ((it-burn) %% thin == 0)){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

save_param_res <- function(param_all_lst){
  
  param_res_lst <- list()
  param_res_lst$mu <- param_all_lst$mu
  param_res_lst$A <- param_all_lst$A
  param_res_lst$B <- param_all_lst$B
  param_res_lst$gamma_beta <- param_all_lst$gamma_beta
  param_res_lst$L <- param_all_lst$L
  param_res_lst$sigma2 <- param_all_lst$sigma2
  param_res_lst$P_star <- param_all_lst$P_star
  param_res_lst$sparsity_matrix <- param_all_lst$sparsity_matrix
  param_res_lst$pivot <- param_all_lst$pivot
  
  return(param_res_lst)
}