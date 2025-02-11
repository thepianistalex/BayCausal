test_that("check B", {
  
  Nit <- 5000
  burn <- 3000
  thin <- 10
  seed <- 0
  n <- 5000
  
  data <- generate_sim_data(seed, n)
  
  prior_lst <- set_prior_default(data)
  mh_setup_lst <- set_mh_default()
  
  init_lst <- set_init_default(seed, data$P, data, prior_lst)
  init_lst$B <- data$B
  init_lst$A <- data$A
  init_lst$L <- data$L
  init_lst$C <- data$C
  init_lst_sigma2 <- data$sigma_e^2
  
  chain_setup_lst <- set_chain_null()
  chain_setup_lst$Nit <- Nit
  chain_setup_lst$burn <- burn
  chain_setup_lst$thin <- thin
  chain_setup_lst$seed <- seed
  
  
  res <- glvcausal_check_B(data, mh_setup_lst, init_lst, prior_lst, chain_setup_lst, FALSE)
  post_mean <- extract_post_mean(res, "B", data$P)
  
  expect_true(are_all_close(post_mean, data$B, abs_tol = 0.1))
  
  }
  )