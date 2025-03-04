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
#' @param param 
#'
#' @returns posterior mean of the parameter
#' @export
extract_post_mean_all <- function(res, param) {
  
  param_samples <- lapply(res, function(x) x[[param]])
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
#' @param res A list containing MCMC samples
#' @param order_flag A logical indicating whether to order the columns of L
#' @param pos_sign A logical indicating whether to enforce positive sign on the first non-zero element of each column of L
#' @param cutoff A numeric value indicating the threshold below which the elements of L are set to zero
#'
#' @returns post processed samples
#' @export
post_process_L <- function(res, order_flag, pos_sign, cutoff = -Inf) {
  
  processed_res <- lapply(res, function(sample) {
    
    if (is.finite(cutoff)) {
      sample$L[abs(sample$L) < cutoff] <- 0
      sample$sparsity_matrix <- ifelse(sample$L == 0, 0, 1)
      sample$pivot <- apply(sample$L, 2, function(col_j) {
        nz <- which(col_j != 0)
        if (length(nz) == 0) NA_integer_ else nz[1]
      })
      nonzero_cols <- !is.na(sample$pivot)
      sample$P_star <- sum(nonzero_cols)
      sample$L <- sample$L[, nonzero_cols, drop = FALSE]
      sample$pivot <- sample$pivot[nonzero_cols]
      sample$sparsity_matrix <- sample$sparsity_matrix[, nonzero_cols, drop = FALSE]
    }
    
    if (pos_sign) {
      for (j in 1:ncol(sample$L)) {
        if (sample$L[sample$pivot[j], j] < 0) {
          sample$L[, j] <- -sample$L[, j]
        }
      }
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
#' @param L_mean The posterior mean of L
#' @param order_flag A logical indicating whether to order the columns of L
#' @param pos_sign A logical indicating whether to enforce positive sign on the first non-zero element of each column of L
#' @param cutoff A numeric value indicating the threshold below which the elements of L are set to zero
#'
#' @returns post processed L posterior mean
#' @export
post_process_L_mean <- function(L_mean, order_flag, pos_sign, cutoff = -Inf) {
  
  processed_L_mean <- L_mean
  
  if (is.finite(cutoff)) {
    processed_L_mean[abs(processed_L_mean) < cutoff] <- 0
  }
  
  pivot <- apply(processed_L_mean, 2, function(col_j) {
    nz <- which(col_j != 0)
    if (length(nz) == 0) {
      NA_integer_  
    } else {
      nz[1] 
    }
  })
  nonzero_cols <- !is.na(pivot)
  processed_L_mean <- processed_L_mean[, nonzero_cols, drop = FALSE]
  pivot <- pivot[nonzero_cols]
  
  if (pos_sign) {
    for (j in 1:ncol(processed_L_mean)) {
      if (processed_L_mean[pivot[j], j] < 0) {
        processed_L_mean[, j] <- -processed_L_mean[, j]
      }
    }
  }
  
  if (order_flag) {
    ord <- order(pivot, decreasing = FALSE)
    processed_L_mean <- processed_L_mean[, ord, drop = FALSE]
  }
  
  return(processed_L_mean)
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