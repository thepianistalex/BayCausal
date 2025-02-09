//[[Rcpp::depends(RcppArmadillo)]]
# include <RcppArmadillo.h>


using namespace std;
using namespace arma;
using namespace Rcpp;



/* 
 helper functions
 */

// [[Rcpp::export]]
double rinvgamma_rcpp(const double& a, const double& b){
  
  return(1/R::rgamma(a, 1/b));
  
}



// [[Rcpp::export]]
arma::mat rmvn_rcpp(const int& n, arma::vec &mean, arma::mat &sigma){
  
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