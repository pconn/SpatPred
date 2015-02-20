#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


arma::vec armaNorm(int n){
  NumericVector x = rnorm(n,0,1);
  arma::vec out(x.begin(), x.size(), false);
  return out;
}

arma::vec armaU(int n){
  NumericVector x = runif(n,0,1);
  arma::vec out(x.begin(), x.size(), false);
  return out;
}
arma::vec GCN(const arma::mat& V, const arma::vec& v){
  arma::mat out = solve(V,v) + solve(chol(V), armaNorm(v.n_elem));
  return out;
}

arma::vec ln_Pois(const arma::vec& k, const arma::vec& ln_lambda){
  return k%ln_lambda - exp(ln_lambda);
}

double ln_norm(const arma::vec& x, const double& sigma_sq){
  return -0.5*(x.size()*log(sigma_sq) + arma::dot(x,x)/sigma_sq);
}

double ln_log_hc(const double& x, const double& scale){
  return x - log(1 + exp(2*x)/(scale*scale));
}


// [[Rcpp::export]]
List PCtimeMCMC(const arma::vec& y, const arma::mat& K, const arma::vec& sampleIndex,
const double& beta_mean, const double& beta_prec, const double& phi_scale, const double& sigma_scale,
const int& block, const int& burn, const int& iter) {
  //matrices for fitting
  int na = K.n_cols;
  int ntot = K.n_rows;
  int ns = y.n_elem;
  int nns = ntot-ns;
  arma::uvec smp = find(sampleIndex==1);
  arma::uvec not_smp = find(sampleIndex==0);
  arma::vec X(ntot, fill::ones);
  
  // Storage matrices
  arma::vec beta_store(iter, fill::zeros);
  arma::mat alpha_store(iter, na, fill::zeros);
  arma::vec phi_store(iter, fill::zeros);
  arma::vec sigma_store(iter, fill::zeros);
  arma::mat pred_store(iter, ntot);
  
  // Z items
  arma::vec z_prop(ns);
  arma::vec U(ns);
  arma::vec MHR_z(ns);
  arma::mat jump_z(iter+burn, ns, fill::zeros);
  arma::uvec jump;
  arma::vec tune_z(ns);
  tune_z.fill(2.4*2.4);
  arma::vec zeta_curr(ns);
  arma::vec zeta_prop(ns);
  
  // Z tuning items
  arma::vec pv_z(ns);
  arma::mat z_smp_store(iter+burn, ns);
  arma::vec r_z(ns);
  arma::uvec i_uvec(1);
  arma::mat tune_store(burn+iter, ns);
  
  // Beta items
  double V_beta_inv;
  double v_beta;
  
  // alpha items
  arma::mat I_alpha(na, na, fill::eye);
  arma::mat V_alpha_inv(na, na);
  arma::vec v_alpha(na);
  
  // phi items
//  arma::vec Kbar(ntot);
//  double V_phi_inv;
//  double v_phi;
  double xi_phi = log(5);
  double phi = exp(xi_phi);
  double phi2_inv = exp(-2*xi_phi);
  double xi_phi_prop;
  double MHR_xi_phi;
  arma::vec xi_phi_store(iter+burn);
  arma::vec jump_xi_phi(iter+burn, fill::zeros);
  double r_xi_phi;
  double tune_xi_phi = 2.4*2.4;
  double pv_xi_phi = 1;

  
  // sigma items
  double xi = 0;
  double sigma = exp(xi);
  double sigma2_inv = 1/(sigma*sigma);
  double xi_prop;
  double sigma_prop;
  double MHR_xi;
  arma::vec xi_store(iter+burn);
  arma::vec jump_xi(iter+burn, fill::zeros);
  double r_xi;
  double tune_xi = 2.4*2.4;
  double pv_xi = 1;
  
  // Initial values
  arma::vec z_curr = log(y+1);
  double beta = log(mean(y));
  arma::vec z(ntot);
  z.elem(not_smp) = beta*ones<vec>(nns);
  arma::vec alpha(na, fill::zeros);
  arma::vec mu_z = beta + K*alpha;
  arma::vec mu_z_smp = mu_z(smp);
//  double sigma2_inv = 100; //(ns-1)/as_scalar((z - X_smp*beta).t()*(z - X_smp*beta));
  pv_z = 1/(y+1) + sigma*sigma;
  
  
  
  //Begin MCMC
  for(int i=0; i<iter+burn; i++){
    
    //update z
    mu_z_smp = mu_z(smp);
    z.elem(not_smp) = mu_z(not_smp) + armaNorm(nns)/sqrt(sigma2_inv);
    z_prop = z_curr + sqrt(tune_z)%sqrt(pv_z)%armaNorm(ns);
    MHR_z = exp( y%z_prop-exp(z_prop) - 0.5*(z_prop-mu_z_smp)%(z_prop-mu_z_smp)*sigma2_inv   
    -(y%z_curr-exp(z_curr))-(-0.5*(z_curr-mu_z_smp)%(z_curr-mu_z_smp)*sigma2_inv) );
    jump = find(armaU(ns)<=MHR_z);
    z_curr.elem(jump) = z_prop.elem(jump);
    i_uvec(0)=i;
    jump_z(i_uvec,jump) = ones<rowvec>(jump.n_elem);
    z.elem(smp) = z_curr;  
    z_smp_store.row(i) = z_curr.t();
    
    // adapt z MH tuning parameter
    if(i>0 & i%block==0){
      r_z = mean(jump_z.submat(i-block, 0, i, ns-1)).t();
      tune_z = exp(log(tune_z) + pow(i/block,-0.5)*(r_z-0.234));
      pv_z = pv_z + pow(i/block,-0.5)*(var(z_smp_store.submat(i-block, 0, i, ns-1)).t() - pv_z);
    }
    tune_store.row(i) = (sqrt(tune_z)%sqrt(pv_z)).t();
    
    // update beta
    V_beta_inv = sigma2_inv*as_scalar(X.t()*X) + beta_prec;
    v_beta = sigma2_inv * as_scalar(X.t()*(z - K*alpha) + beta_prec*beta_mean);
    beta = as<double>(rnorm(1, v_beta/V_beta_inv, 1/sqrt(V_beta_inv)));
    if(i>=burn){beta_store(i-burn) = beta;}
    
    // update alpha
    V_alpha_inv = sigma2_inv*K.t()*K + phi2_inv*I_alpha; 
    v_alpha = sigma2_inv*K.t()*(z-X*beta);
    alpha = GCN(V_alpha_inv, v_alpha);
    if(i>=burn){alpha_store.row(i-burn) = phi*alpha.t();}
    
    // update phi
    xi_phi_prop = xi_phi + R::rnorm(0,sqrt(tune_xi_phi)*sqrt(pv_xi_phi));
    MHR_xi_phi = exp( R::dcauchy(exp(xi_phi_prop), 0, phi_scale, 1) + ln_norm(alpha, exp(2*xi_phi_prop)) - R::dcauchy(exp(xi_phi), 0, phi_scale, 1) - ln_norm(alpha, exp(2*xi_phi)));
    if(R::runif(0,1) <= MHR_xi_phi){
      xi_phi = xi_phi_prop;
      phi = exp(xi_phi);
      phi2_inv = exp(-2*xi_phi);
      jump_xi_phi(i) = 1;
    }
    xi_phi_store(i) = xi_phi;
    if(i>=burn) phi_store(i-burn) = phi;
    
    // adapt sigma MH tuning parameter
    if(i>0 & i%block==0){
      r_xi_phi = mean(jump_xi_phi.subvec(i-block, i));
      tune_xi_phi = exp(log(tune_xi_phi) + pow(i/block,-0.5)*(r_xi_phi-0.234));
      pv_xi_phi = pv_xi_phi + pow(i/block,-0.5)*(var(xi_phi_store.subvec(i-block, i)) - pv_xi_phi);
    }

//    Kbar = K*alpha;
//    V_phi_inv = sigma2_inv*as_scalar(Kbar.t()*Kbar) + 1/(phi_scale*phi_scale);
//    v_phi = sigma2_inv*dot(Kbar, z-X*beta);
//    phi = as<double>(rnorm(1, v_phi/V_phi_inv, 1/sqrt(V_phi_inv)));
//    if(i>=burn) phi_store(i-burn) = sqrt(phi*phi); 
    
    // update xi = log(sigma)
    mu_z = beta + K*alpha;
    xi_prop = xi + R::rnorm(0,sqrt(tune_xi)*sqrt(pv_xi));
    sigma_prop = exp(xi_prop);
    MHR_xi = exp( R::dcauchy(sigma_prop, 0, sigma_scale, 1) + ln_norm(z-mu_z, sigma_prop*sigma_prop) - R::dcauchy(sigma, 0, sigma_scale, 1) - ln_norm(z-mu_z, sigma*sigma));
    if(R::runif(0,1) <= MHR_xi){
      xi = xi_prop;
      sigma = exp(xi);
      sigma2_inv = 1/(sigma*sigma);
      jump_xi(i) = 1;
    }
    xi_store(i) = xi;
    if(i>=burn) sigma_store(i-burn) = sigma;
    
    // adapt sigma MH tuning parameter
    if(i>0 & i%block==0){
      r_xi = mean(jump_xi.subvec(i-block, i));
      tune_xi = exp(log(tune_xi) + pow(i/block,-0.5)*(r_xi-0.234));
      pv_xi = pv_xi + pow(i/block,-0.5)*(var(xi_store.subvec(i-block, i)) - pv_xi);
    }

//    // update sigma2_inv
//    mu_z = beta + K*alpha;
//    sigma2_inv = as<double>(rgamma(1, ntot/2 + a_sigma, as_scalar((z-mu_z).t()*(z-mu_z))/2 + b_sigma));
//    if(i>=burn) sigma_store(i-burn) = 1/sqrt(sigma2_inv);
    
    // make prediction
    if(i>=burn){
      pred_store.row(i-burn) = exp(mu_z).t();
    }
    
    
  }
  
  return Rcpp::List::create(
    Rcpp::Named("beta") = beta_store, 
    Rcpp::Named("alpha") = alpha_store,
    Rcpp::Named("phi")=phi_store,
    Rcpp::Named("sigma")=sigma_store,
    Rcpp::Named("pred")=pred_store,
    Rcpp::Named("jump_z")=jump_z,
    Rcpp::Named("jump_xi")=jump_xi,
    Rcpp::Named("tune")=tune_store
    );  
}

