//[[Rcpp::depends(RcppArmadillo)]]
# include <RcppArmadillo.h>


using namespace std;
using namespace arma;
using namespace Rcpp;



// [[Rcpp::export]]
double rinvgamma_rcpp(const double a, const double b){
  return(1/R::rgamma(a, 1/b));
}

// [[Rcpp::export]]
arma::mat rmvn_rcpp(const int n, arma::vec& mean, arma::mat& sigma){

  int k = sigma.n_cols;
  arma::mat z = arma::randn(n, k);
  arma::mat out = arma::repmat(mean,1,n).t() + z*arma::chol(sigma);
  return(out);
}
