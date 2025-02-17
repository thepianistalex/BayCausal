#' Title
#'
#' @param res 
#' @param param 
#' @param P_star 
#'
#' @returns posterior mean of the parameter
#' @export
extract_post_mean <- function(res, param, P_star) {
  filtered_samples <- Filter(function(x) x[["P_star"]] == P_star, res)
  
  if (length(filtered_samples) == 0) {
    stop("No samples with matching P_star found.")
  }
  
  param_samples <- lapply(filtered_samples, function(x) x[[param]])
  post_mean_param <- Reduce("+", param_samples) / length(param_samples)
  
  return(post_mean_param)
}



#' Title
#'
#' @param res 
#'
#' @returns posterior samples of P_star
#' @export
extract_P_star <- function(res) {
  
  P_star_samples <- unlist(lapply(res, function(x) x[["P_star"]]))
  
  return(P_star_samples)
}



#' Title
#'
#' @param res 
#' @param order_flag 
#' @param cutoff 
#'
#' @returns post processed samples
#' @export
post_process_L <- function(res, order_flag, cutoff = Inf) {
  
  processed_res <- lapply(res, function(sample) {
    
    if (is.finite(cutoff)) {
      sample$L[abs(sample$L) < cutoff] <- 0
      sample$sparsity_matrix <- ifelse(sample$L == 0, 0, 1)
      sample$pivot <- apply(sample$L, 2, function(col_j) {
        nz <- which(col_j != 0)
        if (length(nz) == 0) {
          NA_integer_  
        } else {
          nz[1] 
        }
      })
    }
    
    if (order_flag) {
      ord <- order(sample$pivot, decreasing = FALSE)
      sample$L <- sample$L[, ord, drop = FALSE]
      sample$sparsity_matrix <- sample$sparsity_matrix[, ord, drop = FALSE]
      sample$pivot <- sample$pivot[ord]
    }
    
    
    return(sample)
  })
  
  return(processed_res)
}



#' Title
#'
#' @param v 
#'
#' @returns mode of the vector
#' @export
get_mode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}