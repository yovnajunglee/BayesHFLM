// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppNumerical)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(roptim)]]

// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <RcppNumerical.h>
#include "RcppEigen.h"
#include <roptim.h>
#include <cmath> 
using namespace arma;
using namespace Numer;
using namespace Eigen; 
using namespace roptim;

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//
// another simple example: outer product of a vector, 
// returning a matrix
//
// [[Rcpp::export]]
arma::mat rcpparma_outerproduct(const arma::colvec & x) {
    arma::mat m = x * x.t();
    return m;
}

// and the inner product returns a scalar
//
// [[Rcpp::export]]
double rcpparma_innerproduct(const arma::colvec & x) {
    double v = arma::as_scalar(x.t() * x);
    return v;
}


// and we can use Rcpp::List to return both at the same time
//
// [[Rcpp::export]]
Rcpp::List rcpparma_bothproducts(const arma::colvec & x) {
    arma::mat op = x * x.t();
    double    ip = arma::as_scalar(x.t() * x);
    return Rcpp::List::create(Rcpp::Named("outer")=op,
                              Rcpp::Named("inner")=ip);
}

// and we can use Rcpp::List to return both at the same time
//
// [[Rcpp::export]]
arma::mat myfirstfunc() {
    arma::mat m1 = arma::eye<arma::mat>(3, 3);
    arma::mat m2 = arma::eye<arma::mat>(3, 3);
    
    return m1 + 3 * (m1 + m2);
}

// Sample from a MVN posterior distribution
// of the form [theta|.] ~ N(Qc^-1el, Qc^-1)
// using a Cholesky Decomposition
// [[Rcpp::export]]
arma::mat sampleMVN(arma::vec ell, arma::mat Qc) {
    arma::mat Ql = arma::chol(Qc);
    arma::mat z = arma::randn<arma::mat>(Qc.n_cols);
    arma::mat vals = arma::solve(Ql, arma::solve(Ql.t(), ell) + z);
    return vals;
}

// Make matrix N x T out of a vector NT X 1
// [[Rcpp::export]]
arma::mat data_matrix(arma::mat mymat, int nc, int nr){
  
  int pos = 0;
  arma::mat datmat(nr, nc);
  
  for (int i = 0; i < nr ; ++i){
    //std::cout << mymat.rows(pos, pos + nc - 1) << std::endl;
    datmat.row(i) = mymat.rows(pos, (pos+nc-1)).t();
    pos = pos+nc;
    //std::cout << pos << std::endl;
    
  }
  
  return datmat;
}




// Construct X matrix for 
// HFLM 
// [[Rcpp::export]]
arma::mat constructMatrix(arma::mat Xvec, arma::colvec taus,
                          arma::colvec Ss, arma::mat Fk, arma::mat Fk1, 
                          arma::mat Psi, double delta, int nobs){
    // Evaluate no. of observations
    int ntau = taus.n_elem; //Xmat.n_cols;
    //int nobs = Xmat.n_rows;
    int totals = ntau*nobs;
    int U = Psi.n_cols;
    int K = Fk.n_cols;
    // Initialise matrix
    arma::mat Xmat = data_matrix(Xvec, ntau, nobs);
    arma::mat Xtil(totals, U);
    
    int temp = 0;

    for (int i = 0; i < nobs; ++i){
        for (int j = 0; j < ntau; ++j){
            //double xij = Xmat(i, j) ;
            double tauval = taus(j);
            for (int u = 0; u < U; ++u){
                double intval = 0;
                for (int s = 0; s < ntau; ++s){
                    double psival = Psi(s, u);
                    double sval = Ss(s);
                        if (sval >= (tauval - delta) && sval <= (tauval)){
                            intval = intval + Xmat(i,s)*psival;
                        }else{
                            intval = intval + 0;
                        }
                }
                Xtil(temp,u) = intval/ntau;
            }
           temp++;
        }
    }
    
    arma::mat tau_mat = repmat(taus, nobs ,1);
    //std::cout << "in construct matrix ... " << std::endl;
    
    //for(int i = 0; i < totals; ++i){
      //  double xij = Xvec(i);
        //double tauval = tau_mat(i);
        //for (int u = 0; u < U; ++u){
          //  double intval = 0; 
            //for (int s = 0; s < ntau; ++s){
              //  double psival = Psi(s, u);
              //  double sval = Ss(s);
                //if (sval >= std::max((tauval - delta),0.0) && sval <= (tauval)){
                  //  intval = intval + xij*psival;
              //  }else{
              //      intval = intval + 0;
              //  }
          //  }
          //  Xtil(i,u) = intval/ntau;
      //  }
  //  }
    
   arma::mat FF  = repmat(Fk, nobs, 1);
   arma::mat XX(totals, U*K);
        
   int temp1 = 0;
   for(int k = 0 ; k < K; ++k){
       for(int u = 0; u < U; ++u){
           XX.col(temp1) = Xtil.col(u) % FF.col(k);
           temp1 = temp1 + 1;
       }
   }
   arma::mat FF1  = repmat(Fk1, nobs, 1);
   //arma::mat tau2 = tau_mat%tau_mat;
   //arma::mat tau3 = tau2%tau_mat;
   
   arma::mat Rmat  =  arma::join_horiz(FF1, XX);
    //arma::join_horiz(arma::ones(totals,1), FF1, XX);
   //Rmat = arma::join_horiz(Rmat,FF1, XX);
  return Rmat;
}



// MCMC sampler
// [[Rcpp::export]]
Rcpp::List mcmc_sampler(arma::colvec Yvec, arma::mat Ymat,
                        arma::colvec Xvec, arma::mat Xmat, 
                        arma::colvec taus, arma::colvec Ss,
                        arma::mat Fk, arma::mat Fk1, arma::mat Psi, arma::mat Phi,
                        arma::mat Dmu, arma::mat Db, arma::mat Dalpha,
                        arma::mat Dc, double delta,
                        double i1, double i2, double a, double b, int niter){
  
  // Set parameters
  int nobs = Ymat.n_rows;
  int ntaus = Ymat.n_cols;
  int K = Fk.n_cols;
  int K1 = Fk1.n_cols;
  int U = Psi.n_cols;
  int Kphi  = Phi.n_cols;
  int totals = nobs*ntaus;
  
  // Create Phi matrix
  // arma::mat Phi_rep  = repmat(Phi, nobs, 1);
  int dc = Phi.n_cols;
  // Initialise parameters
  //std::cout << "init Rmat ... " << std::endl;
  arma::mat Rmat = constructMatrix(Xvec, taus, Ss, Fk, Fk1, Psi, delta, nobs); 
  int dr = Rmat.n_cols; 
  //std::cout << "Init alpha ... " << std::endl;
  arma::mat alpha(dr, niter) ; 
  alpha.col(0) = arma::randu(dr);
  //std::cout << "Init c... " << std::endl;
  
  arma::mat csims(dc, niter) ; 
  csims.col(0) = arma::randu(dc);
  //std::cout << "Init sigma ... " << std::endl;
  
  arma::mat sigma_e(1, niter);
  sigma_e.col(0) = 10 ; //arma::randu(1);
  
  arma::mat sigma_b(1, niter);
  sigma_b.col(0) = 10 ; //arma::randu(1);
  
  arma::mat sigma_mu(1, niter);
  sigma_mu.col(0) = 10; // arma::randu(1);
  
  arma::mat sigma_c(1, niter);
  sigma_c.col(0) = 10; // arma::randu(1);
  
  arma::mat sigma_v(1, niter);
  sigma_v.col(0) = 10; // arma::randu(1);
  
  // Define matrices used below
  arma::mat Qc(Kphi, Kphi);
  arma::vec ell(Kphi); 
  arma::mat c_sample(Kphi, 1);
  arma::mat Xsample(totals, 1);
  arma::mat sigs(dr, dr);
  arma::mat Omegaf(dr, dr);
  arma::mat SigmaAlph(dr, dr);
  arma::mat muAlpha(dr, 1);
  arma::mat alpha_sample(dr, 1);
  double shape_e;
  double rate_e;
  double sigma_e_sample;
  double shape_mu;
  double rate_mu;
  double sigma_mu_sample;
  double shape_b;
  double rate_b;
  double sigma_b_sample;
  double shape_c;
  double rate_c;
  double sigma_c_sample;
  double shape_v;
  double rate_v;
  double sigma_v_sample;
  
  for (int iter = 1; iter < niter; ++iter){
    std::cout << "Iteration "<< iter << std::endl;
    // Sample c for smooth X
    std::cout << "Sampling c ... " << std::endl;
    std::cout << "Qc ... " << std::endl;
    
    //this takes long
    Qc = Dc*(1/as_scalar(sigma_c.col(iter - 1))) + (1/as_scalar(sigma_v.col(iter - 1)))*Phi.t()*Phi+ 0.0000001*arma::eye(dc,dc);;
    //  std::cout << "ell... " << std::endl;
    
    ell = 1/(as_scalar(sigma_v.col(iter - 1)))*Phi.t()*Xvec;
    // std::cout << "Sample mvn ... " << std::endl;
    
    c_sample = sampleMVN(ell, Qc);
    Xsample = Phi*c_sample;
    std::cout << "Constructing new R ..." << std::endl;
    Rmat  = constructMatrix(Xsample, taus, Ss, Fk, Fk1, Psi, delta, nobs);
    
    // Create vector sigmas
    sigs = join_vert(ones(K1, 1)*(1/as_scalar(sigma_mu.col(iter - 1))),
                     ones(U*K, 1)*(1/as_scalar(sigma_b.col(iter -1))));
    std::cout << "Sampling alphas ... " << std::endl;
    std::cout << "omegaf ... " << std::endl;
    
    Omegaf = repmat(sigs,1, Dalpha.n_cols)%Dalpha; 
    std::cout << " sigmaalphas ... " << std::endl;
    
    SigmaAlph = Omegaf + ((1/(as_scalar(sigma_e.col(iter - 1))))*Rmat.t()*Rmat) + 0.0000001*arma::eye(Rmat.n_cols, Rmat.n_cols);
    std::cout << "mu  alphas ... " << std::endl;
    std::cout << Omegaf.is_symmetric() << std::endl;
    std::cout << ((1/(as_scalar(sigma_e.col(iter - 1))))*Rmat.t()*Rmat).is_symmetric() << std::endl;
    
    muAlpha = (1/(as_scalar(sigma_e.col(iter - 1))))*Rmat.t()*Yvec;
    alpha_sample = sampleMVN(muAlpha, SigmaAlph);
    std::cout << alpha_sample << std::endl;
    // Generate the sigmas
    //std::cout << "Sampling sigma  ... " << std::endl;
    //std::cout << "shape e ... " << std::endl;
    shape_e = ntaus*nobs/2 + i1;
    
    //std::cout << "rate e ... " << std::endl;
    rate_e = i2 + 0.5*as_scalar((Yvec - (Rmat*alpha_sample)).t()*(Yvec - (Rmat*alpha_sample)));
    sigma_e_sample = 1/arma::randg<double>( distr_param(shape_e,1/rate_e));
    //(Rcpp::rgamma(1, shape_e, 1/rate_e));
    //std::cout << "shape mu ... " << std::endl;
    shape_mu = K/2 + a;
    
    //std::cout << "rate mu ... " << std::endl;
    rate_mu = b + 0.5*as_scalar((alpha_sample.rows(0, K1-1).t()*Dmu*alpha_sample.rows(0, K1-1)));
    sigma_mu_sample =  1/arma::randg<double>(distr_param(shape_mu,1/rate_mu));
    
    std::cout << "shape b ... " << std::endl;
    shape_b = U*K/2 + a;
    std::cout << "rate b ... " << std::endl;
    
    rate_b = b + 0.5*as_scalar((alpha_sample.rows(K1, U*K + K1 - 1).t()*Db*alpha_sample.rows(K1, U*K + K1 - 1)));
    sigma_b_sample =  1/arma::randg<double>( distr_param(shape_b,1/rate_b));
    std::cout << "shape v ... " << std::endl;
    
    shape_v =  ntaus*nobs/2 + i1;
    std::cout << "rate v ... " << std::endl;
    
    rate_v = i2 + 0.5*as_scalar(((Xvec - Xsample).t()*(Xvec - Xsample)));
    sigma_v_sample = 1/arma::randg<double>( distr_param(shape_v,1/rate_v));
    std::cout << "shape c ... " << std::endl;
    
    shape_c = (Kphi - nobs)/2 + a;
    std::cout << "rate c ... " << std::endl;
    
    rate_c = b + 0.5*as_scalar((c_sample.t()*Dc*c_sample));
    sigma_c_sample = 1/arma::randg<double>( distr_param(shape_c,1/rate_c));
    
    // Store values
    alpha.col(iter) = alpha_sample;
    csims.col(iter) = c_sample;
    
    sigma_e.col(iter) = sigma_e_sample;
    sigma_mu.col(iter) = sigma_mu_sample;
    sigma_b.col(iter) = sigma_b_sample;
    sigma_v.col(iter) = sigma_v_sample;
    sigma_c.col(iter) = sigma_c_sample;
    
  }
  
  return Rcpp::List::create(Rcpp::Named("alphas")=alpha,
                            Rcpp::Named("c")=csims,
                            Rcpp::Named("sigmae") = sigma_e,
                            Rcpp::Named("sigmamu") = sigma_mu,
                            Rcpp::Named("sigmab") = sigma_b,
                            Rcpp::Named("sigmac") = sigma_c,
                            Rcpp::Named("sigmav") = sigma_v);
}


// MCMC sampler 
// Parameter reparameterisation to obtain diagonal penalty matrices
// for computation
// [[Rcpp::export]]
Rcpp::List mcmc_sampler2(arma::colvec Yvec, arma::mat Ymat,
                         arma::colvec Xvec, arma::mat Xmat, 
                         arma::colvec taus, arma::colvec Ss,
                         arma::mat Fk, arma::mat Fk1, arma::mat Psi, arma::mat Phi,
                         arma::mat Dmu, arma::mat Db, arma::mat Dalpha,
                         arma::mat Dc, double delta,
                         double i1, double i2, double a, double b, int niter){
  
  // Set parameters
  int nobs = Ymat.n_rows;
  int ntaus = Ymat.n_cols;
  int K = Fk.n_cols;
  int K1 = Fk1.n_cols;
  int U = Psi.n_cols;
  int Kphi  = Phi.n_cols;
  int totals = nobs*ntaus;
  
  // Create Phi matrix
  // arma::mat Phi_rep  = repmat(Phi, nobs, 1);
  int dc = Phi.n_cols;
  // Initialise parameters
  //std::cout << "init Rmat ... " << std::endl;
  arma::mat Rmat = constructMatrix(Xvec, taus, Ss, Fk, Fk1, Psi, delta, nobs); 
  int dr = Rmat.n_cols; 
  //std::cout << "Init alpha ... " << std::endl;
  arma::mat alpha(dr, niter) ; 
  alpha.col(0) = arma::randu(dr);
  //std::cout << "Init c... " << std::endl;
  
  arma::mat csims(dc, niter) ; 
  csims.col(0) = arma::randu(dc);
  //std::cout << "Init sigma ... " << std::endl;
  
  arma::mat sigma_e(1, niter);
  sigma_e.col(0) = 10 ; //arma::randu(1);
  
  arma::mat sigma_b(1, niter);
  sigma_b.col(0) = 10 ; //arma::randu(1);
  
  arma::mat sigma_mu(1, niter);
  sigma_mu.col(0) = 10; // arma::randu(1);
  
  arma::mat sigma_c(1, niter);
  sigma_c.col(0) = 10; // arma::randu(1);
  
  arma::mat sigma_v(1, niter);
  sigma_v.col(0) = 10; // arma::randu(1);
  
  // Define matrices used below
  arma::mat Qc(Kphi, Kphi);
  arma::vec ell(Kphi); 
  arma::mat c_sample(Kphi, 1);
  arma::mat Xsample(totals, 1);
  arma::mat sigs(dr, dr);
  arma::mat Omegaf(dr, dr);
  arma::mat SigmaAlph(dr, dr);
  arma::mat muAlpha(dr, 1);
  arma::mat alpha_sample(dr, 1);
  double shape_e;
  double rate_e;
  double sigma_e_sample;
  double shape_mu;
  double rate_mu;
  double sigma_mu_sample;
  double shape_b;
  double rate_b;
  double sigma_b_sample;
  double shape_c;
  double rate_c;
  double sigma_c_sample;
  double shape_v;
  double rate_v;
  double sigma_v_sample;
  
  for (int iter = 1; iter < niter; ++iter){
    std::cout << "Iteration "<< iter << std::endl;
    // Sample c for smooth X
    std::cout << "Sampling c ... " << std::endl;
    std::cout << "Qc ... " << std::endl;
    
    //this takes long
    Qc = (1/as_scalar(sigma_c.col(iter - 1))) + (1/as_scalar(sigma_v.col(iter - 1)))*Phi.t()*Phi+ 0.0000001*arma::eye(dc,dc);;
    //  std::cout << "ell... " << std::endl;
    
    ell = 1/(as_scalar(sigma_v.col(iter - 1)))*Phi.t()*Xvec;
    // std::cout << "Sample mvn ... " << std::endl;
    
    c_sample = sampleMVN(ell, Qc);
    Xsample = Phi*c_sample;
    std::cout << "Constructing new R ..." << std::endl;
    Rmat  = constructMatrix(Xsample, taus, Ss, Fk, Fk1, Psi, delta, nobs);
    
    // Create vector sigmas
    sigs = join_vert(ones(K1, 1)*(1/as_scalar(sigma_mu.col(iter - 1))),
                     ones(U*K, 1)*(1/as_scalar(sigma_b.col(iter -1))));
    std::cout << "Sampling alphas ... " << std::endl;
    std::cout << "omegaf ... " << std::endl;
    
    Omegaf = repmat(sigs,1, Dalpha.n_cols)%Dalpha; 
    std::cout << " sigmaalphas ... " << std::endl;
    
    SigmaAlph = Omegaf + ((1/(as_scalar(sigma_e.col(iter - 1))))*Rmat.t()*Rmat) + 0.0000001*arma::eye(Rmat.n_cols, Rmat.n_cols);
    std::cout << "mu  alphas ... " << std::endl;
    std::cout << Omegaf.is_symmetric() << std::endl;
    std::cout << ((1/(as_scalar(sigma_e.col(iter - 1))))*Rmat.t()*Rmat).is_symmetric() << std::endl;
    
    muAlpha = (1/(as_scalar(sigma_e.col(iter - 1))))*Rmat.t()*Yvec;
    alpha_sample = sampleMVN(muAlpha, SigmaAlph);
    std::cout << alpha_sample << std::endl;
    // Generate the sigmas
    //std::cout << "Sampling sigma  ... " << std::endl;
    //std::cout << "shape e ... " << std::endl;
    shape_e = ntaus*nobs/2 + i1;
    
    //std::cout << "rate e ... " << std::endl;
    rate_e = i2 + 0.5*as_scalar((Yvec - (Rmat*alpha_sample)).t()*(Yvec - (Rmat*alpha_sample)));
    sigma_e_sample = 1/arma::randg<double>( distr_param(shape_e,1/rate_e));
    //(Rcpp::rgamma(1, shape_e, 1/rate_e));
    //std::cout << "shape mu ... " << std::endl;
    shape_mu = K/2 + a;
    
    //std::cout << "rate mu ... " << std::endl;
    rate_mu = b + 0.5*as_scalar((alpha_sample.rows(0, K1-1).t()*Dmu*alpha_sample.rows(0, K1-1)));
    sigma_mu_sample =  1/arma::randg<double>(distr_param(shape_mu,1/rate_mu));
    
    std::cout << "shape b ... " << std::endl;
    shape_b = U*K/2 + a;
    std::cout << "rate b ... " << std::endl;
    
    rate_b = b + 0.5*as_scalar((alpha_sample.rows(K1, U*K + K1 - 1).t()*Db*alpha_sample.rows(K1, U*K + K1 - 1)));
    sigma_b_sample =  1/arma::randg<double>( distr_param(shape_b,1/rate_b));
    std::cout << "shape v ... " << std::endl;
    
    shape_v =  ntaus*nobs/2 + i1;
    std::cout << "rate v ... " << std::endl;
    
    rate_v = i2 + 0.5*as_scalar(((Xvec - Xsample).t()*(Xvec - Xsample)));
    sigma_v_sample = 1/arma::randg<double>( distr_param(shape_v,1/rate_v));
    std::cout << "shape c ... " << std::endl;
    
    shape_c = (Kphi - nobs)/2 + a;
    std::cout << "rate c ... " << std::endl;
    
    rate_c = b + 0.5*as_scalar((c_sample.t()*Dc*c_sample));
    sigma_c_sample = 1/arma::randg<double>( distr_param(shape_c,1/rate_c));
    
    // Store values
    alpha.col(iter) = alpha_sample;
    csims.col(iter) = c_sample;
    
    sigma_e.col(iter) = sigma_e_sample;
    sigma_mu.col(iter) = sigma_mu_sample;
    sigma_b.col(iter) = sigma_b_sample;
    sigma_v.col(iter) = sigma_v_sample;
    sigma_c.col(iter) = sigma_c_sample;
    
  }
  
  return Rcpp::List::create(Rcpp::Named("alphas")=alpha,
                            Rcpp::Named("c")=csims,
                            Rcpp::Named("sigmae") = sigma_e,
                            Rcpp::Named("sigmamu") = sigma_mu,
                            Rcpp::Named("sigmab") = sigma_b,
                            Rcpp::Named("sigmac") = sigma_c,
                            Rcpp::Named("sigmav") = sigma_v);
}

// Reparam theta(s-t)
//[[Rcpp::export]]
arma::mat constructMatrix2(arma::mat Xvec, arma::colvec taus,
                          arma::colvec Ss, arma::mat Fk, arma::mat Fk1, 
                          arma::mat Psi, double delta, int nobs){
  // Evaluate no. of observations
  int ntau = taus.n_elem; //Xmat.n_cols;
  //int nobs = Xmat.n_rows;
  int totals = ntau*nobs;
  int U = Psi.n_cols;
  int K = Fk.n_cols;
  // Initialise matrix
  arma::mat Xmat = data_matrix(Xvec, ntau, nobs);
  arma::mat Xtil(totals, U);
  
  int temp = 0;
  double sval ;
  double tauval;
  double intval;
  int v; 
  //std::cout << "first for loop ... " << std::endl;
  for (int i = 0; i < nobs; ++i){
    for (int j = 0; j < ntau; ++j){
      //double xij = Xmat(i, j) ;
      tauval= taus(j);
      for (int u = 0; u < U; ++u){
        intval = 0;
        v = 0;
        for (int s = 0; s < ntau; ++s){
          //std::cout << "s "<< s << "v "<< v << std::endl;
          sval = Ss(s);
          if (sval >= (tauval - delta) && sval <= (tauval)){
            intval = intval + Xmat(i,s)*Psi(v,u); // X_i(s)*Psi_u(v=s-t)
            ++v ; 
          }else{
            intval = intval + 0;
          }
        }
        Xtil(temp,u) = intval/ntau;
      }
      temp++;
    }
  }
  
  arma::mat tau_mat = repmat(taus, nobs ,1);
  //std::cout << "in construct matrix ... " << std::endl;
  
  //for(int i = 0; i < totals; ++i){
  //  double xij = Xvec(i);
  //double tauval = tau_mat(i);
  //for (int u = 0; u < U; ++u){
  //  double intval = 0; 
  //for (int s = 0; s < ntau; ++s){
  //  double psival = Psi(s, u);
  //  double sval = Ss(s);
  //if (sval >= std::max((tauval - delta),0.0) && sval <= (tauval)){
  //  intval = intval + xij*psival;
  //  }else{
  //      intval = intval + 0;
  //  }
  //  }
  //  Xtil(i,u) = intval/ntau;
  //  }
  //  }
  
  arma::mat FF  = repmat(Fk, nobs, 1);
  arma::mat XX(totals, U*K);
  
  int temp1 = 0;
  for(int k = 0 ; k < K; ++k){
    for(int u = 0; u < U; ++u){
      XX.col(temp1) = Xtil.col(u) % FF.col(k);
      temp1 = temp1 + 1;
    }
  }
  arma::mat FF1  = repmat(Fk1, nobs, 1);
  //arma::mat tau2 = tau_mat%tau_mat;
  //arma::mat tau3 = tau2%tau_mat;
  
  arma::mat Rmat  =  arma::join_horiz(FF1, XX);
  //arma::join_horiz(arma::ones(totals,1), FF1, XX);
  //Rmat = arma::join_horiz(Rmat,FF1, XX);
  return Rmat;
}


// MCMC sampler 
// estimate b(s,t)= theta(s-t,t)
// with uniform prior
// for computation
// [[Rcpp::export]]
Rcpp::List mcmc_sampler3(arma::colvec Yvec, arma::mat Ymat,
                         arma::colvec Xvec, arma::mat Xmat, 
                         arma::colvec taus, arma::colvec Ss,
                         arma::mat Fk, arma::mat Fk1, arma::mat Psi, arma::mat Phi,
                         arma::mat Dmu, arma::mat Db, arma::mat Dalpha,
                         arma::mat Dc, double delta,
                         double i1, double i2, double a, double b, int niter){
  
  // Set parameters
  int nobs = Ymat.n_rows;
  int ntaus = Ymat.n_cols;
  int K = Fk.n_cols;
  int K1 = Fk1.n_cols;
  int U = Psi.n_cols;
  int Kphi  = Phi.n_cols;
  int totals = nobs*ntaus;
  
  // Create Phi matrix
  // arma::mat Phi_rep  = repmat(Phi, nobs, 1);
  int dc = Phi.n_cols;
  // Initialise parameters
  std::cout << "init Rmat ... " << std::endl;
  arma::mat Rmat = constructMatrix2(Xvec, taus, Ss, Fk, Fk1, Psi, delta, nobs); 
  int dr = Rmat.n_cols; 
  //std::cout << "Init alpha ... " << std::endl;
  arma::mat alpha(dr, niter) ; 
  alpha.col(0) = arma::randu(dr);
  //std::cout << "Init c... " << std::endl;
  
  arma::mat csims(dc, niter) ; 
  csims.col(0) = arma::randu(dc);
  //std::cout << "Init sigma ... " << std::endl;
  
  arma::mat sigma_e(1, niter);
  sigma_e.col(0) = 10 ; //arma::randu(1);
  
  arma::mat sigma_b(1, niter);
  sigma_b.col(0) = 10 ; //arma::randu(1);
  
  arma::mat sigma_mu(1, niter);
  sigma_mu.col(0) = 10; // arma::randu(1);
  
  arma::mat sigma_c(1, niter);
  sigma_c.col(0) = 10; // arma::randu(1);
  
  arma::mat sigma_v(1, niter);
  sigma_v.col(0) = 10; // arma::randu(1);
  
  // Define matrices used below
  arma::mat Qc(Kphi, Kphi);
  arma::vec ell(Kphi); 
  arma::mat c_sample(Kphi, 1);
  arma::mat Xsample(totals, 1);
  arma::mat sigs(dr, dr);
  arma::mat Omegaf(dr, dr);
  arma::mat SigmaAlph(dr, dr);
  arma::mat muAlpha(dr, 1);
  arma::mat alpha_sample(dr, 1);
  double shape_e;
  double rate_e;
  double sigma_e_sample;
  double shape_mu;
  double rate_mu;
  double sigma_mu_sample;
  double shape_b;
  double rate_b;
  double sigma_b_sample;
  double shape_c;
  double rate_c;
  double sigma_c_sample;
  double shape_v;
  double rate_v;
  double sigma_v_sample;
  
  for (int iter = 1; iter < niter; ++iter){
    std::cout << "Iteration "<< iter << std::endl;
    // Sample c for smooth X
    std::cout << "Sampling c ... " << std::endl;
    std::cout << "Qc ... " << std::endl;
    
    //this takes long
    Qc = Dc*(1/as_scalar(sigma_c.col(iter - 1))) + (1/as_scalar(sigma_v.col(iter - 1)))*Phi.t()*Phi+ 0.0000001*arma::eye(dc,dc);;
    //  std::cout << "ell... " << std::endl;
    
    ell = 1/(as_scalar(sigma_v.col(iter - 1)))*Phi.t()*Xvec;
    // std::cout << "Sample mvn ... " << std::endl;
    
    c_sample = sampleMVN(ell, Qc);
    Xsample = Phi*c_sample;
    std::cout << "Constructing new R ..." << std::endl;
    Rmat  = constructMatrix2(Xsample, taus, Ss, Fk, Fk1, Psi, delta, nobs);
    
    // Create vector sigmas
    sigs = join_vert(ones(K1, 1)*(1/as_scalar(sigma_mu.col(iter - 1))),
                     ones(U*K, 1)*(1/as_scalar(sigma_b.col(iter -1))));
    std::cout << "Sampling alphas ... " << std::endl;
    std::cout << "omegaf ... " << std::endl;
    
    Omegaf = repmat(sigs,1, Dalpha.n_cols)%Dalpha; 
    std::cout << " sigmaalphas ... " << std::endl;
    
    SigmaAlph = Omegaf + ((1/(as_scalar(sigma_e.col(iter - 1))))*Rmat.t()*Rmat) + 0.0000001*arma::eye(Rmat.n_cols, Rmat.n_cols);
    std::cout << "mu  alphas ... " << std::endl;
    std::cout << Omegaf.is_symmetric() << std::endl;
    std::cout << ((1/(as_scalar(sigma_e.col(iter - 1))))*Rmat.t()*Rmat).is_symmetric() << std::endl;
    
    muAlpha = (1/(as_scalar(sigma_e.col(iter - 1))))*Rmat.t()*Yvec;
    alpha_sample = sampleMVN(muAlpha, SigmaAlph);
    std::cout << alpha_sample << std::endl;
    // Generate the sigmas
    //std::cout << "Sampling sigma  ... " << std::endl;
    //std::cout << "shape e ... " << std::endl;
    shape_e = ntaus*nobs/2 + i1;
    
    //std::cout << "rate e ... " << std::endl;
    rate_e = i2 + 0.5*as_scalar((Yvec - (Rmat*alpha_sample)).t()*(Yvec - (Rmat*alpha_sample)));
    sigma_e_sample = 1/arma::randg<double>( distr_param(shape_e,1/rate_e));
    //(Rcpp::rgamma(1, shape_e, 1/rate_e));
    //std::cout << "shape mu ... " << std::endl;
    shape_mu = K/2 + a;
    
    //std::cout << "rate mu ... " << std::endl;
    rate_mu = b + 0.5*as_scalar((alpha_sample.rows(0, K1-1).t()*Dmu*alpha_sample.rows(0, K1-1)));
    sigma_mu_sample =  1/arma::randg<double>(distr_param(shape_mu,1/rate_mu));
    
    std::cout << "shape b ... " << std::endl;
    shape_b = U*K/2 + a;
    std::cout << "rate b ... " << std::endl;
    
    rate_b = b + 0.5*as_scalar((alpha_sample.rows(K1, U*K + K1 - 1).t()*Db*alpha_sample.rows(K1, U*K + K1 - 1)));
    sigma_b_sample =  1/arma::randg<double>( distr_param(shape_b,1/rate_b));
    std::cout << "shape v ... " << std::endl;
    
    shape_v =  ntaus*nobs/2 + i1;
    std::cout << "rate v ... " << std::endl;
    
    rate_v = i2 + 0.5*as_scalar(((Xvec - Xsample).t()*(Xvec - Xsample)));
    sigma_v_sample = 1/arma::randg<double>( distr_param(shape_v,1/rate_v));
    std::cout << "shape c ... " << std::endl;
    
    shape_c = (Kphi - nobs)/2 + a;
    std::cout << "rate c ... " << std::endl;
    
    rate_c = b + 0.5*as_scalar((c_sample.t()*Dc*c_sample));
    sigma_c_sample = 1/arma::randg<double>( distr_param(shape_c,1/rate_c));
    
    // Store values
    alpha.col(iter) = alpha_sample;
    csims.col(iter) = c_sample;
    
    sigma_e.col(iter) = sigma_e_sample;
    sigma_mu.col(iter) = sigma_mu_sample;
    sigma_b.col(iter) = sigma_b_sample;
    sigma_v.col(iter) = sigma_v_sample;
    sigma_c.col(iter) = sigma_c_sample;
    
  }
  
  return Rcpp::List::create(Rcpp::Named("alphas")=alpha,
                            Rcpp::Named("c")=csims,
                            Rcpp::Named("sigmae") = sigma_e,
                            Rcpp::Named("sigmamu") = sigma_mu,
                            Rcpp::Named("sigmab") = sigma_b,
                            Rcpp::Named("sigmac") = sigma_c,
                            Rcpp::Named("sigmav") = sigma_v);
}


