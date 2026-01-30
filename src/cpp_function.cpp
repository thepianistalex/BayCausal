//[[Rcpp::depends(RcppArmadillo)]]
# include <RcppArmadillo.h>


using namespace std;
using namespace arma;
using namespace Rcpp;



/* 
 helper functions
 */

// [[Rcpp::export]]
double rinvgamma_rcpp(const double &a, const double &b){
  
  return(1/R::rgamma(a, 1/b));
  
}



// [[Rcpp::export]]
arma::mat rmvn_rcpp(const int &n, arma::vec &mean, arma::mat &sigma){
  
  int k = sigma.n_cols; 
  mat z = randn(n, k);
  mat out = repmat(mean,1,n).t() + z*chol(sigma);
  return(out);
  
}



// [[Rcpp::export]]
double rinvgaussian_rcpp(const double &mu, const double &lambda){
  
  double out;
  double nu = R::rnorm(0, 1);
  double y = pow(nu, 2);
  double x = mu + ((pow(mu,2)*y)/(2*lambda)) - (mu/(2*lambda))*
    sqrt(4*mu*lambda*y + pow(mu,2)*pow(y,2));
  double z = R::runif(0, 1);
  
  if (z > (mu/(mu + x))){ out = mu*mu/x; 
  }else{ out = x; }
  return(out);
  
}



// [[Rcpp::export]]
arma::vec compute_Y_tilde_q(const int &q, const arma::mat &B, const arma::mat &A, 
                            const arma::mat &L, const arma::mat &C, 
                            const arma::vec &mu, const int &Q, const int &S,
                            const int &P, const arma::mat &Y, const arma::mat &X,
                            const bool &minus_mu, const bool &minus_BY, const bool &minus_AX, const bool &minus_LC) {
  
  vec Y_tilde_q = Y.col(q);  
  
  if (minus_mu) {
    Y_tilde_q -= mu(q);
  }
  
  if (minus_BY) {
    vec B_q = B.row(q).t();   // shape: Q x 1
    // Y is n x Q, B_q is Q x 1 => (Y * B_q) is n x 1
    Y_tilde_q -= Y * B_q;
  }
  
  if (minus_AX) {
    vec A_q = A.row(q).t();   // shape: S x 1
    Y_tilde_q -= X * A_q;
  }
  
  if (minus_LC) {
    vec L_q = L.row(q).t();   // shape: P x 1
    // C is n x P, L_q is P x 1 => (C * L_q) is n x 1
    Y_tilde_q -= C * L_q;
  }
  
  return Y_tilde_q;
}



/* 
 likelihood functions
 */

// [[Rcpp::export]]
double compute_loglik_B_rcpp(const int &q, const arma::mat &B, const arma::mat &A, 
                             const arma::mat &L, const arma::mat &C, const arma::vec &mu, 
                             const arma::vec &sigma2, const arma::mat &Y, const arma::mat &X){
  
  double loglik = 0;
  double n = Y.n_rows;
  double Q = Y.n_cols;
  double S = X.n_cols;
  double P = L.n_cols;
  double jacob_factor;
  int q_cpp = q - 1;
  
  // calculate the jacobian matrix 
  mat jacob; jacob.eye(Q,Q); jacob -= B;
  jacob_factor = abs(det(jacob));
  loglik += n*log(jacob_factor);
  vec Y_tilde_q = compute_Y_tilde_q(q_cpp, B, A, L, C, mu, Q, S, P, Y, X, true, true, true, true);
  
  for (int i=0; i<n; i++){
    double b  = 2.0 * sqrt(sigma2(q_cpp));
    double yy = abs(Y_tilde_q(i)) / b;
    double log_density = -log(2.0 * b) - yy;
    loglik += log_density;
  }
  
  return(loglik);
}



// [[Rcpp::export]]
double compute_loglik_B_t_rcpp(const int &q, const arma::mat &B, const arma::mat &A, 
                             const arma::mat &L, const arma::mat &C, const arma::vec &mu, 
                             const arma::vec &sigma2, const arma::mat &Y, const arma::mat &X,
                             const double &nu_t){
  
  double loglik = 0;
  double n = Y.n_rows;
  double Q = Y.n_cols;
  double S = X.n_cols;
  double P = L.n_cols;
  double jacob_factor;
  int q_cpp = q - 1;
  
  // calculate the jacobian matrix 
  mat jacob; jacob.eye(Q,Q); jacob -= B;
  jacob_factor = abs(det(jacob));
  loglik += n*log(jacob_factor);
  vec Y_tilde_q = compute_Y_tilde_q(q_cpp, B, A, L, C, mu, Q, S, P, Y, X, true, true, true, true);
  
  double sigma = sqrt(sigma2(q_cpp));
  double const_term = R::lgammafn((nu_t + 1.0) / 2.0) - R::lgammafn(nu_t / 2.0) - 0.5 * log(nu_t * M_PI) - log(sigma);

  for (int i=0; i<n; i++){
    double resid = Y_tilde_q(i);
    double z = resid / sigma;
    double log_density = const_term - ((nu_t + 1.0) / 2.0) * log(1.0 + (z * z) / nu_t);
    loglik += log_density;
  }
  
  return(loglik);
}



// [[Rcpp::export]]
double compute_marginal_loglik_rcpp(const arma::mat &sparsity_matrix, const int &q, const arma::mat &C, const arma::mat &tau,
                                    const arma::mat &L_0, const arma::mat &Y, const arma::mat &B, const arma::mat &X, const arma::mat &A, const arma::vec &mu,
                                    const double &a_sigma, const double &b_sigma) {
  
  int q_cpp = q - 1;
  int n = Y.n_rows;
  double log_margin_lik;
  
  vec Y_q_tilda = Y.col(q_cpp) - mu(q_cpp) - Y * B.row(q_cpp).t() - X * A.row(q_cpp).t();
  uvec non_zero_p = find(sparsity_matrix.row(q_cpp) == 1);
  
  if (non_zero_p.n_elem != 0) {
    
    mat C_q = C.cols(non_zero_p);
    rowvec L_0_row = L_0.row(q_cpp);
    vec L_0_row_elems = conv_to<vec>::from(L_0_row.elem(non_zero_p));
    mat L_0_q_inv = diagmat(1.0 / L_0_row_elems);
    
    mat M = C_q.t() * (tau.col(q_cpp) % C_q.each_col()) + L_0_q_inv;
    mat V_q_n = inv_sympd(M);
    vec b_q = C_q.t() * (tau.col(q_cpp) % Y_q_tilda);
    
    double a_sigma_n = a_sigma + n / 2.0;
    double term1 = dot(Y_q_tilda, tau.col(q_cpp) % Y_q_tilda);
    double term2 = dot(b_q, V_q_n * b_q);
    double b_sigma_n = b_sigma + 0.5 * (term1 - term2);
    
    double sign = 0.0, log_det_V_q_n = 0.0;
    log_det(log_det_V_q_n, sign, V_q_n);
    
    double sum_log_L_0_q = sum(log(L_0_row_elems));
    
    log_margin_lik = 0.5 * (log_det_V_q_n - sum_log_L_0_q) + a_sigma * log(b_sigma) -
      a_sigma_n * log(b_sigma_n) + lgamma(a_sigma_n) - lgamma(a_sigma);
    
  } else {
    // When non_zero_p is empty
    double a_sigma_n = a_sigma + n / 2.0;
    double term1 = dot(Y_q_tilda, tau.col(q_cpp) % Y_q_tilda);
    double b_sigma_n = b_sigma + 0.5 * term1;
    
    log_margin_lik = a_sigma * log(b_sigma) - a_sigma_n * log(b_sigma_n) +
      lgamma(a_sigma_n) - lgamma(a_sigma);
  }
  
  return log_margin_lik;
}



// [[Rcpp::export]]
double compute_loglik_rjmcmc_cpp(int &q, const arma::mat &B, const arma::mat &A, 
                                 const arma::mat &L, const arma::mat &C, 
                                 const arma::vec &mu, const arma::vec &sigma2, 
                                 const arma::mat &Y, const arma::mat &X
) {
  // make the index 0-based
  int q_cpp = q - 1;
  int n = Y.n_rows;
  arma::vec Y_tilde_q = Y.col(q_cpp) 
    - mu(q_cpp)
    - Y * B.row(q_cpp).t() 
    - X * A.row(q_cpp).t() 
    - C * L.row(q_cpp).t();
    
    double scale = 2.0 * sqrt(sigma2(q_cpp));
    double loglik = - (double)n * log(2.0 * scale) 
      - arma::sum(arma::abs(Y_tilde_q)) / scale;
    
    return loglik;
}

// [[Rcpp::export]]
double compute_loglik_rjmcmc_t_rcpp(int &q, const arma::mat &B, const arma::mat &A, 
                                 const arma::mat &L, const arma::mat &C, 
                                 const arma::vec &mu, const arma::vec &sigma2, 
                                 const arma::mat &Y, const arma::mat &X,
                                 const double &nu_t
) {
  // make the index 0-based
  int q_cpp = q - 1;
  int n = Y.n_rows;
  arma::vec Y_tilde_q = Y.col(q_cpp) 
    - mu(q_cpp)
    - Y * B.row(q_cpp).t() 
    - X * A.row(q_cpp).t() 
    - C * L.row(q_cpp).t();
    
  double sigma = sqrt(sigma2(q_cpp));
  double const_term = R::lgammafn((nu_t + 1.0) / 2.0) - R::lgammafn(nu_t / 2.0) - 0.5 * log(nu_t * M_PI) - log(sigma);
  
  double loglik = 0;
  for (int i=0; i<n; i++){
    double resid = Y_tilde_q(i);
    double z = resid / sigma;
    loglik += const_term - ((nu_t + 1.0) / 2.0) * log(1.0 + (z * z) / nu_t);
  }
    
  return loglik;
}



/* 
 update functions
 */

// [[Rcpp::export]]
arma::mat update_mu_rcpp(const arma::mat &B, const arma::mat &A, const arma::mat &L, 
                         const arma::mat &C,
                         const arma::mat &tau, const arma::vec &sigma2,
                         const arma::mat &Y, const arma::mat &X){
  
  double n = Y.n_rows;
  double Q = Y.n_cols;
  double S = X.n_cols;
  double P = L.n_cols;
  double Y_tilde_iq, tau_Ytilde_sum, sigma2_tilde_sum, mu_n, V_n;
  
  vec mu_update(Q);
  
  for (int q=0; q<Q; q++){
    
    tau_Ytilde_sum = 0; 
    sigma2_tilde_sum = 0;
    vec mu_placeholder;
    vec Y_tilde_q = compute_Y_tilde_q(q, B, A, L, C, mu_placeholder, Q, S, P, Y, X, false, true, true, true);
    
    for (int i=0; i<n; i++){
      
      Y_tilde_iq = Y_tilde_q(i);
      
      tau_Ytilde_sum += Y_tilde_iq*tau(i,q);
      sigma2_tilde_sum += tau(i,q)/sigma2(q);
    }
    
    // Gaussian full conditional distribution for mu
    V_n = 1/(sigma2_tilde_sum + 1/100);
    mu_n = V_n*tau_Ytilde_sum/sigma2(q);
    mu_update(q) = R::rnorm(mu_n, sqrt(V_n));
  }
  
  return(mu_update);
}



// [[Rcpp::export]]
arma::mat update_A_rcpp(const arma::mat &gamma_alpha, const arma::mat &nu_alpha,
                        const arma::mat &B, const arma::mat &L, const arma::mat &C,
                        const arma::vec &mu, const arma::mat &tau, const arma::vec &sigma2,
                        const arma::mat &Y, const arma::mat &X){
  
  double n = Y.n_rows;
  double Q = Y.n_cols;
  double S = X.n_cols;
  double P = L.n_cols;
  double Y_tilde_iq;
  
  rowvec X_i;
  vec tau_X_Ytilde_sum (S); 
  vec mu_n(S); 
  mat V_n(S,S);
  mat tau_X2_sum(S,S);
  mat V_alpha(S,S);
  mat A_update(Q, S);
  mat A_placeholder;
  
  for (int q=0; q<Q; q++){
    
    tau_X2_sum.fill(0); 
    tau_X_Ytilde_sum.fill(0);
    // need to plug in sth for the place of A for the function to work
    vec Y_tilde_q = compute_Y_tilde_q(q, B, A_placeholder, L, C, mu, Q, S, P, Y, X, true, true, false, true);
    
    for (int i=0; i<n; i++){
      
      X_i = X.submat(i, 0, i, S-1);
      tau_X2_sum += X_i.t()*X_i*tau(i,q);
      
      Y_tilde_iq = Y_tilde_q(i);
      tau_X_Ytilde_sum += conv_to<vec>::from(X_i*Y_tilde_iq*tau(i,q));
    }
    
    // Gaussian full conditional distribution for A
    V_alpha = diagmat(1 / (gamma_alpha.row(q).t() % nu_alpha.row(q).t()));
    V_n = inv_sympd(tau_X2_sum/sigma2(q) + V_alpha);
    mu_n = V_n*(tau_X_Ytilde_sum /sigma2(q));
    A_update.row(q) = conv_to<rowvec>::from(rmvn_rcpp(1, mu_n, V_n));
  }
  
  return(A_update);
}



// [[Rcpp::export]]
arma::mat update_L_free_rcpp(const arma::vec &mu, const arma::mat &A, const arma::mat &B,
                             const arma::mat &C, const arma::mat &tau, const arma::vec &sigma2,
                             const arma::mat &Y, const arma::mat &X) {
  
  int Q = mu.n_elem;
  int P_star = C.n_cols;
  mat L_update(Q, P_star, fill::zeros);
  
  for (int q = 0; q < Q; q++) {
    mat D_L_q_inv = 0.001 * eye(P_star, P_star);
    
    mat L_placeholder;
    vec Y_tilde_q = compute_Y_tilde_q(q, B, A, L_placeholder, C, mu, Q, X.n_cols, P_star, Y, X, true, true, true, false);
    
    vec weight = tau.col(q) / sigma2(q);
    
    // In R: t(C) %*% ((tau[,q]/sigma2[q])*C)
    // In Armadillo, multiply each column of C by the weight vector (elementwise multiplication on rows)
    mat weighted_C = C.each_col() % weight;
    mat temp = D_L_q_inv + C.t() * weighted_C;
    
    mat V_n = inv_sympd(temp);
    vec mu_n = V_n * ( C.t() * (Y_tilde_q % weight) );
    
    L_update.row(q) = conv_to<rowvec>::from(rmvn_rcpp(1, mu_n, V_n));
  }
  
  return L_update;
}



// [[Rcpp::export]]
arma::mat update_L_uni_rcpp(const arma::vec &mu, const arma::mat &A, const arma::mat &B,
                            const arma::mat &C, const arma::mat &tau, const arma::vec &sigma2,
                            const arma::mat &Y, const arma::mat &X) {
  // Same interface as update_L_free_rcpp, but constrain L to be diagonal:
  // only L(q,q) is updated (when that column exists); all off-diagonals are zero.
  int Q = mu.n_elem;
  int P_star = C.n_cols;
  int S = X.n_cols;
  int diag_len = std::min(Q, P_star);
  arma::mat L_update(Q, P_star, arma::fill::zeros);

  // prior precision for each diagonal element (matching free-L prior scale 0.001 I)
  const double prior_prec = 0.0000001;

  for (int q = 0; q < diag_len; ++q) {
    // Build Y_tilde without the LC term (so L contribution is excluded)
    arma::mat L_placeholder; // unused in compute_Y_tilde_q
    arma::vec Y_tilde_q = compute_Y_tilde_q(q, B, A, L_placeholder, C, mu, Q, S, P_star, Y, X,
                                            true,  // minus_mu
                                            true,  // minus_BY
                                            true,  // minus_AX
                                            false  // minus_LC (exclude)
                                            );

    arma::vec weight = tau.col(q) / sigma2(q);
    arma::vec C_q = C.col(q);

    // Scalar Normal posterior for L(q,q)
    double V_n_inv = prior_prec + arma::dot(C_q, C_q % weight);
    double V_n = 1.0 / V_n_inv;
    double mu_n = V_n * arma::dot(C_q, Y_tilde_q % weight);

    // Draw from N(mu_n, V_n)
    double sample = R::rnorm(mu_n, std::sqrt(V_n));
    L_update(q, q) = sample;
  }

  return L_update;
}

// [[Rcpp::export]]
arma::mat update_C_rcpp(const arma::vec &mu, const arma::mat &A, const arma::mat &B, 
                        const arma::mat &L, const arma::mat &tau, 
                        const arma::vec &sigma2,
                        const arma::mat &Y, const arma::mat &X) {
  
  int P_star = L.n_cols;
  int n = Y.n_rows;
  vec Y_tilde_i, mu_n;
  mat Sigma_i_inv, V_n;
  mat C_update(n, P_star, fill::none);
  
  for (int i = 0; i < n; i++) {
    
    Sigma_i_inv = diagmat(tau.row(i).t() / sigma2);
    
    Y_tilde_i = Y.row(i).t() - mu - B * Y.row(i).t() - A * X.row(i).t();
    
    V_n = inv_sympd(diagmat(ones(P_star)) + L.t() * Sigma_i_inv * L);
    mu_n = V_n * L.t() * Sigma_i_inv * Y_tilde_i;
    
    C_update.row(i) = conv_to<rowvec>::from(rmvn_rcpp(1, mu_n, V_n));
  }
  
  return C_update;
}



// [[Rcpp::export]]
arma::mat update_tau_rcpp(const arma::mat &B, const arma::mat &A, const arma::mat &L, 
                                const arma::mat &C, const arma::vec &mu, const arma::vec &sigma2,
                                const arma::mat &Y, const arma::mat &X){

  double n = Y.n_rows;
  double Q = Y.n_cols;
  double S = X.n_cols;
  double P = L.n_cols;
  double eps = 1e-3;
  double lambda_iq = 1.0/4;
  double tau_update_iq;
  mat tau_update(n, Q); tau_update.fill(datum::nan);

  for (int q=0; q<Q; q++){
    vec Y_tilde_q = compute_Y_tilde_q(q, B, A, L, C, mu, Q, S, P, Y, X, true, true, true, true);
    vec mu_q = sqrt(sigma2(q)) / (2.0 * arma::abs(Y_tilde_q));
    
    for (int i=0; i<n; i++){
      tau_update_iq = rinvgaussian_rcpp(mu_q(i), lambda_iq);
      tau_update(i,q) = max(tau_update_iq, eps);
    }
    
  }

  return(tau_update);
}



// [[Rcpp::export]]
arma::vec update_sigma2_rcpp(const arma::mat &B, const arma::mat &A, const arma::mat &L, const arma::mat &C, const arma::vec &mu, 
                             const arma::mat &tau,
                             const arma::mat &Y, const arma::mat &X,
                             const double &a_sigma, const double &b_sigma){
  
  double n = Y.n_rows;
  double Q = Y.n_cols;
  double S = X.n_cols;
  double P = L.n_cols;
  double tau_Ytilde_sum, Y_tilde_iq;
  double a_sigma_star, b_sigma_star;
  vec sigma2_update(Q);
  
  for (int q=0; q<Q; q++){
    
    tau_Ytilde_sum = 0;
    vec Y_tilde_q = compute_Y_tilde_q(q, B, A, L, C, mu, Q, S, P, Y, X, true, true, true, true);
    
    for (int i=0; i<n; i++){
      Y_tilde_iq = Y_tilde_q(i);
      tau_Ytilde_sum += pow(Y_tilde_iq,2)*tau(i,q);
    }
    // Inverse-Gamma posterior distribution
    a_sigma_star = a_sigma + n/2;
    b_sigma_star = b_sigma + tau_Ytilde_sum/2;
    sigma2_update(q) = rinvgamma_rcpp(a_sigma_star, b_sigma_star);
  }
  
  return(sigma2_update);
}



// [[Rcpp::export]]
arma::mat update_tau_t_rcpp(const arma::mat &B, const arma::mat &A, const arma::mat &L, 
                                const arma::mat &C, const arma::vec &mu, const arma::vec &sigma2,
                                const arma::mat &Y, const arma::mat &X, const double &nu_t){

  double n = Y.n_rows;
  double Q = Y.n_cols;
  double S = X.n_cols;
  double P = L.n_cols;

  mat tau_update(n, Q, fill::none);

  for (int q=0; q<Q; q++){
    // For tau update we need residuals. All terms included.
    vec Y_tilde_q = compute_Y_tilde_q(q, B, A, L, C, mu, Q, S, P, Y, X, true, true, true, true);
    
    // lambda_iq ~ Gamma((nu+1)/2, rate = (y-mu)^2/(2*sigma2) + nu/2)
    // R::rgamma uses scale = 1/rate
    double shape = (nu_t + 1.0) / 2.0;
    
    for (int i=0; i<n; i++){
       double resid_sq = pow(Y_tilde_q(i), 2);
       double rate = resid_sq / (2.0 * sigma2(q)) + nu_t / 2.0;
       double scale = 1.0 / rate;
       tau_update(i,q) = R::rgamma(shape, scale); 
    }
  }

  return(tau_update);
}