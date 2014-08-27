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

arma::vec ln_norm(const arma::vec& x, const arma::vec& mu, const arma::vec& sigma2_inv){
  return -0.5*(x-mu)%(x-mu)%sigma2_inv;
}


// [[Rcpp::export]]
List PCtimeMCMC(const arma::vec& y, const arma::mat& K_lf, const arma::mat& K_hf, const arma::vec& sampleIndex,
const double& beta_mean, const double& beta_prec, const double& phi_lf_scale, const double& phi_hf_scale, const double& a_sigma,
const double & b_sigma, const int& block, const int& burn, const int& iter, const bool& sample_sigma,
const bool& sample_hf) {
  //matrices for fitting
  int na_lf = K_lf.n_cols;
  int na_hf = K_hf.n_cols;
  int ntot = K_lf.n_rows;
  int ns = y.n_elem;
  int nns = ntot-ns;
  arma::uvec smp = find(sampleIndex==1);
  arma::uvec not_smp = find(sampleIndex==0);
  arma::vec X(ntot, fill::ones);
  
  // Storage matrices
  arma::vec beta_store(iter, fill::zeros);
  arma::mat alpha_lf_store(iter, na_lf, fill::zeros);
  arma::mat alpha_hf_store(iter, na_hf, fill::zeros);
  arma::vec phi_lf_store(iter, fill::zeros);
  arma::vec phi_hf_store(iter, fill::zeros);
  arma::vec sigma_store(iter, fill::zeros);
  arma::mat pred_store(iter, ntot);
  
  // Z items
  arma::vec z_prop(ns);
  arma::vec U(ns);
  arma::vec MHR(ns);
  arma::mat jump_idx(iter+burn, ns, fill::zeros);
  arma::uvec jump;
  arma::vec tune(ns);
  tune.fill(2.4*2.4);
  arma::vec zeta_curr(ns);
  arma::vec zeta_prop(ns);
  
  // Z tuning items
  arma::vec sigma_0(ns);
  arma::mat z_smp_store(iter+burn, ns);
  arma::vec r_hat(ns);
  arma::uvec i_uvec(1);
  
  // Beta items
  double V_beta_inv;
  double v_beta;
  
  // alpha items
  arma::mat Sigma_alpha_lf(na_lf, na_lf, fill::eye);
  arma::mat Sigma_alpha_hf(na_hf, na_hf, fill::eye);
  arma::mat V_alpha_lf_inv(na_lf, na_lf);
  arma::vec v_alpha_lf(na_lf);
  arma::mat V_alpha_hf_inv(na_hf, na_hf);
  arma::vec v_alpha_hf(na_hf);
  
  // phi items
  arma::vec Kbar(ntot);
  double V_phi_lf_inv;
  double v_phi_lf;
  double V_phi_hf_inv;
  double v_phi_hf;
  
  // Initial values
  arma::vec z_curr = log(y+1);
  double beta = log(mean(y));
  arma::vec z(ntot);
  z.elem(not_smp) = beta*ones<vec>(nns);
  double sigma2_inv = 1000000; //(ns-1)/as_scalar((z - X_smp*beta).t()*(z - X_smp*beta));
  double phi_lf = 0;
  double phi_hf = 0;
  arma::vec alpha_lf(na_lf, fill::zeros);
  arma::vec alpha_hf(na_hf, fill::zeros);
  arma::vec mu_z = beta + phi_lf*K_lf*alpha_lf + phi_hf*K_hf*alpha_hf;
  arma::vec mu_z_smp = mu_z(smp);
  sigma_0 = 1/sqrt(exp(z_curr) + sigma2_inv);
  
  
  
  arma::vec eigenvals;
  arma::mat eigenvecs;
  
  
  arma::mat tune_store(burn+iter, ns);
  
  
  //Begin MCMC
  for(int i=0; i<iter+burn; i++){
    
    //update z
    mu_z_smp = mu_z(smp);
    z.elem(not_smp) = mu_z(not_smp) + armaNorm(nns)/sqrt(sigma2_inv);
    z_prop = z_curr + sqrt(tune)%sigma_0%armaNorm(ns);
    MHR = exp( y%z_prop-exp(z_prop) - 0.5*(z_prop-mu_z_smp)%(z_prop-mu_z_smp)*sigma2_inv   
    -(y%z_curr-exp(z_curr))-(-0.5*(z_curr-mu_z_smp)%(z_curr-mu_z_smp)*sigma2_inv) );
    jump = find(armaU(ns)<=MHR);
    z_curr.elem(jump) = z_prop.elem(jump);
    i_uvec(0)=i;
    jump_idx(i_uvec,jump) = ones<rowvec>(jump.n_elem);
    z.elem(smp) = z_curr;  
    z_smp_store.row(i) = z_curr.t();
    
    // adapt z MH tuning parameter
    if(i>0 & i%block==0){
      r_hat = mean(jump_idx.submat(i-block, 0, i, ns-1)).t();
      tune = exp(log(tune) + 2*pow(i/block,-0.8)*(r_hat-0.234));
      sigma_0 = sigma_0 + pow(i/block,-0.8)*(stddev(z_smp_store.submat(i-block, 0, i, ns-1)).t() - sigma_0);
    }
    tune_store.row(i) = (sqrt(tune)%sigma_0).t();
    
    // update beta
    V_beta_inv = sigma2_inv*as_scalar(X.t()*X) + beta_prec;
    v_beta = sigma2_inv * as_scalar(X.t()*(z - phi_lf*K_lf*alpha_lf - phi_hf*K_hf*alpha_hf)) + beta_prec*beta_mean;
    beta = as<double>(rnorm(1, v_beta/V_beta_inv, 1/sqrt(V_beta_inv)));
    if(i>=burn){beta_store(i-burn) = beta;}
    
    // update alpha_lf
    V_alpha_lf_inv = phi_lf*phi_lf*sigma2_inv*K_lf.t()*K_lf + Sigma_alpha_lf; 
    v_alpha_lf = phi_lf*sigma2_inv*K_lf.t()*(z-X*beta-phi_hf*K_hf*alpha_hf);
    alpha_lf = GCN(V_alpha_lf_inv, v_alpha_lf);
    if(i>=burn){alpha_lf_store.row(i-burn) = phi_lf*alpha_lf.t();}
    
    // update alpha_hf
    if(sample_hf){
      V_alpha_hf_inv = phi_hf*phi_hf*sigma2_inv*K_hf.t()*K_hf + Sigma_alpha_hf; 
      v_alpha_hf = phi_hf*sigma2_inv*K_hf.t()*(z-X*beta-phi_lf*K_lf*alpha_lf);
      alpha_hf = GCN(V_alpha_hf_inv, v_alpha_hf);
      if(i>=burn){alpha_hf_store.row(i-burn) = phi_hf*alpha_hf.t();}
    }
    
    // update phi_lf
    Kbar = K_lf*alpha_lf;
    V_phi_lf_inv = sigma2_inv*as_scalar(Kbar.t()*Kbar) + 1/(phi_lf_scale*phi_lf_scale);
    v_phi_lf = sigma2_inv*dot(Kbar, z-X*beta-phi_hf*K_hf*alpha_hf);
    phi_lf = as<double>(rnorm(1, v_phi_lf/V_phi_lf_inv, 1/sqrt(V_phi_lf_inv)));
    if(i>=burn) phi_lf_store(i-burn) = sqrt(phi_lf*phi_lf); 
    
    // update phi_hf
    if(sample_hf){
      Kbar = K_hf*alpha_hf;
      V_phi_hf_inv = sigma2_inv*as_scalar(Kbar.t()*Kbar) + 1/(phi_hf_scale*phi_hf_scale);
      v_phi_hf = sigma2_inv*dot(Kbar, z-X*beta-phi_lf*K_lf*alpha_lf);
      phi_hf = as<double>(rnorm(1, v_phi_hf/V_phi_hf_inv, 1/sqrt(V_phi_hf_inv)));
      if(i>=burn) phi_hf_store(i-burn) = sqrt(phi_hf*phi_hf); 
    }
    
    // update sigma2_inv
    mu_z = beta + phi_lf*K_lf*alpha_lf + phi_hf*K_hf*alpha_hf;
    if(sample_sigma){
      sigma2_inv = as<double>(rgamma(1, ntot/2 + a_sigma, as_scalar((z-mu_z).t()*(z-mu_z))/2 + b_sigma));
      if(i>=burn) sigma_store(i-burn) = 1/sqrt(sigma2_inv);
    }
    
    // make prediction
    if(i>=burn){
      pred_store.row(i-burn) = exp(mu_z).t();
    }
    
    
  }
  
  return Rcpp::List::create(
    Rcpp::Named("beta") = beta_store, 
    Rcpp::Named("alpha_lf") = alpha_lf_store,
    Rcpp::Named("alpha_hf") = alpha_hf_store,
    Rcpp::Named("phi_lf")=phi_lf_store,
    Rcpp::Named("phi_hf")=phi_hf_store,
    Rcpp::Named("sigma")=sigma_store,
    Rcpp::Named("pred")=pred_store,
    Rcpp::Named("jump_idx")=jump_idx,
    Rcpp::Named("tune")=tune_store
    );  
}

