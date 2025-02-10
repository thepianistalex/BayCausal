save_check <- function(it, burn, thin){
  if((it > burn) & ((it-burn) %% thin == 0)){
    return(TRUE)
  } else {
    return(FALSE)
  }
}