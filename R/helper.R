compute_padd <- function(p, pivot, sparsity_matrix, pa){
  # Check add
  if(pivot[p] == 1){
    add_applicable <- FALSE
  } else{
    addable_pivot <- setdiff(1:(pivot[p]-1), pivot[-p])
    if(length(addable_pivot) > 0){
      add_applicable <- TRUE
    } else{
      add_applicable <- FALSE
    }
  }

  # check delete
  second_non_zero <- which(sparsity_matrix[,p] != 0)[2]
  if(is.na(second_non_zero)){
    del_applicable <- FALSE
  } else{
    if(second_non_zero %in% pivot){
      del_applicable <- FALSE
    } else{
      del_applicable <- TRUE
    }
  }

  # Decide padd
  if(add_applicable & del_applicable){
    padd <- pa
  } else if(add_applicable & !del_applicable){
    padd <- 1
  } else if(!add_applicable & del_applicable){
    padd <- 0
  } else{
    padd <- NA
  }

  if(add_applicable){
    return(list(padd = padd, addable_pivot = addable_pivot))
  } else{
    return(list(padd = padd, addable_pivot = NULL))
  }

}



switch_2_columns_pivots <- function(sparsity_matrix, pivot, p, column_2_switch){

  pivot_new <- pivot
  sparsity_matrix_new <- sparsity_matrix

  pivot_2_switch <- pivot[column_2_switch]
  involved_row_index <- min(pivot[p], pivot_2_switch):max(pivot[p], pivot_2_switch)
  diff_row_indices <- sparsity_matrix[involved_row_index, p] != sparsity_matrix[involved_row_index, column_2_switch]
  diff_row_indices <- involved_row_index[diff_row_indices]
  sparsity_matrix_new <- sparsity_matrix
  sparsity_matrix_new[diff_row_indices, p] <- sparsity_matrix[diff_row_indices, column_2_switch]
  sparsity_matrix_new[diff_row_indices, column_2_switch] <- sparsity_matrix[diff_row_indices, p]

  pivot_new[p] <- pivot_2_switch
  pivot_new[column_2_switch] <- pivot[p]

  return(list(sparsity_matrix_new = sparsity_matrix_new, pivot_new = pivot_new, diff_row_indices = diff_row_indices))
}


compute_L0 <- function(kappa, P_star, data){
  L_0 <- matrix(kappa, ncol(data$Y), P_star)
  return(L_0)
}



compute_psplit <- function(pivot, sparsity_matrix, ps, H){
  # if psplit is 0, then only a merge move is possible
  # if psplit is 1, then only a split move is possible
  # if psplit is ps, then both split and merge moves are possible
  # if psplit is NA, then no move is possible

  n_active_columns <- sum(colSums(sparsity_matrix) > 1)
  n_sp_columns <- sum(colSums(sparsity_matrix) == 1)
  n_zero_columns <- H - n_active_columns - n_sp_columns

  # Check split
  if(n_zero_columns == 0){
    split_applicable <- FALSE
  } else{
    split_applicable <- TRUE
  }

  # Check merge
  ## See if this is still valid MCMC by adding (n_active_columns == 0 & n_sp_columns == 1)
  if(n_sp_columns == 0  | (n_active_columns == 0 & n_sp_columns == 1)){
    merge_applicable <- FALSE
  } else{
    merge_applicable <- TRUE
  }

  # Decide psplit
  if(split_applicable & merge_applicable){
    psplit <- ps
  } else if(split_applicable & !merge_applicable){
    psplit <- 1
  } else if(!split_applicable & merge_applicable){
    psplit <- 0
  } else{
    psplit <- NA
  }

  return(psplit)

}
