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

