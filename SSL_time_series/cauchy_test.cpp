#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

double ln_norm(const arma::vec& x, const arma::vec& mu, const double& sigma_sq){
  return -0.5*(x.size()*log(sigma_sq) + arma::dot((x-mu),(x-mu))/sigma_sq);
}

// [[Rcpp::export]]

arma::vec cauchy_test(const arma::vec& z){
  arma::vec sigma_store(10000);
  arma::vec mu(z.size(), fill::zeros);
  double sigma = 0.3;
  double sigma_prop;
  double mh;
  for(int i=0; i<10000; i++){
    sigma_prop = sigma + R::rnorm(0,0.1);
     mh = exp(R::dcauchy(sigma_prop,0,1,1) + ln_norm(z,mu,sigma_prop*sigma_prop) - R::dcauchy(sigma,0,1,1) - ln_norm(z,mu,sigma*sigma));
     if(R::runif(0,1) <= mh){
      sigma = sigma_prop;
    }
    sigma_store(i) = fabs(sigma);
  }
  return sigma_store;
}

//

/*** R
z = rnorm(1000, 0, 0.3)
sigma = cauchy_test(z)
plot(sigma, type='l')
abline(h=sd(z), col="red")
*/
