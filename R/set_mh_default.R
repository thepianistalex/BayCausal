set_mh_default <- function(){

  mh_setup_lst <- list()

  mh_setup_lst$B_step <- 0.02
  mh_setup_lst$a1_step <- 0.1
  mh_setup_lst$a2_step <- 0.1
  mh_setup_lst$pshift <- 1/3
  mh_setup_lst$pswitch <- 1/3
  mh_setup_lst$pa <- 0.5
  mh_setup_lst$ps <- 0.5

  return(mh_setup_lst)
}
