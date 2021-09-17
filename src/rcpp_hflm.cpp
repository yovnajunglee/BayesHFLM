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
#include <iostream>
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
    ////std::cout << mymat.rows(pos, pos + nc - 1) << std::endl;
    datmat.row(i) = mymat.rows(pos, (pos+nc-1)).t();
    pos = pos+nc;
    ////std::cout << pos << std::endl;
    
  }
  
  return datmat;
}


// check null
// [[Rcpp::export]]
double is_null(Rcpp::Nullable<Rcpp::NumericVector> r){
  double bound = 0; 
  
  return bound;
}

// check null
// [[Rcpp::export]]
double isLag(Rcpp::Nullable<Rcpp::NumericVector> x_) {
  Rcpp::NumericVector x;
  double bound;
  if (x_.isNotNull()) {
    x = x_;
    bound = x[0];
  }
  else{ bound = 1;}
  
  return bound;
}


// Construct X matrix for 
// HFLM 
// [[Rcpp::export]]
arma::mat constructMatrix(arma::mat Xvec, arma::colvec taus,
                          arma::colvec Ss, arma::mat Fk, arma::mat Fk1, 
                          arma::mat Psi, int nobs, double lower_bound, arma::mat Zmat){
    // Evaluate no. of observations
    int ntau = taus.n_elem; //Xmat.n_cols;
    int ns = Ss.n_elem; //Xmat.n_cols;
    std::cout << lower_bound << std::endl;
    //int nobs = Xmat.n_rows;
    int totals = ntau*nobs;
    int U = Psi.n_cols;
    int K = Fk.n_cols;
    // Initialise matrix
    //std::cout << "reshape" << std::endl;
    arma::mat Xmat =  reshape(Xvec, ns, nobs).t();
    //data_matrix(Xvec, ntau, nobs);
    arma::mat Xtil(totals, U);
    
   // double lower_bound = 1; 
    
    // if (is_null(delta) == TRUE ){
    //   lower_bound = 1; 
    // } else {
    //   lower_bound = delta ; 
    // }
    // 
   // if tensor products
      int temp = 0;
  
      for (int i = 0; i < nobs; ++i){
        //std::cout << i << std::endl;
          for (int j = 0; j < ntau; ++j){
              //double xij = Xmat(i, j) ;
            double tauval = taus(j);
             // std::cout << j << std::endl;
              for (int u = 0; u < U; ++u){
                  double intval = 0;
                if(tauval == Ss(0) & lower_bound == 1){
                  intval = intval + Xmat(i,0)*Psi(0, u)*(taus(1)-taus(0));
                }else{
                  for (int s = 0; s < ns; ++s){
                   // std::cout << s << std::endl;
                      double psival = Psi(s, u);
                      double sval = Ss(s);
                          if (sval < (tauval) & sval >= (tauval-lower_bound)){
                              intval = intval + Xmat(i,s)*psival*(taus(1)-taus(0));
                          }

                          
                          //else if (sval == lower_bound ){
                          //  intval = intval + Xmat(i,s)*psival*(taus(1)-taus(0));
                        //  }
                          else {
                              intval = intval + 0;
                          }
                  } }
                  Xtil(temp,u) = intval;
              }
              temp++;
              
            }
          
      }
     arma::mat FF  = repmat(Fk, nobs, 1);
     arma::mat XX(totals, U*K);
          
     int temp1 = 0;
     for(int k = 0 ; k < K; ++k){
         for(int u = 0; u < U; ++u){
             XX.col(temp1) = Xtil.col(u) % FF.col(k);
             temp1 = temp1 + 1;
         }
     }
     
    // -- End of outer produce
  
    
   //arma::mat FF1  = repmat(Fk1, nobs, 1);
   //arma::mat tau2 = tau_mat%tau_mat;
   //arma::mat tau3 = tau2%tau_mat;

   arma::mat Rmat  = XX*Zmat ; //arma::join_horiz(FF1, XX); // arma::join_horiz(arma::ones(totals,1), FF1, XX); 
   //arma::join_horiz(FF1, XX);
    //
   //Rmat = arma::join_horiz(Rmat,FF1, XX);
  return Rmat;
}


// Function to calculate thinplate regression splines basis
// [[Rcpp::export]]
double tprs(double x1, double x2, arma::colvec knot){
  double dist = sqrt(pow(x1-knot(0), 2) + pow(x2-knot(1), 2));
  double tprs;
  if (dist  == 0){
    tprs = 0;
  }
  else {
    tprs = pow(dist,2)*log(dist);
  }
  
  return tprs;
}

// MCMC sampler
// [[Rcpp::export]]
arma::mat constructTPRS(arma::mat Xvec, arma::colvec taus,
                        arma::colvec Ss, int nobs, arma::mat knots, arma::mat Zmat, double delta){
  // Evaluate no. of observations
  int ntau = taus.n_elem; //Xmat.n_cols;
  int ns = Ss.n_elem;
  //int nobs = Xmat.n_rows;
  int totals = ntau*nobs;
  int U  = knots.n_cols;
  // Initialise matrix
  //std::cout << "reshape" << std::endl;
  arma::mat Xmat =  reshape(Xvec, ns, nobs).t();
  //data_matrix(Xvec, ntau, nobs);
  arma::mat Xtil(totals, U + 3);
  std::cout << totals << std::endl;
  std::cout << U << std::endl;
  
  int temp = 0;
  for(int i = 0; i < nobs; ++i) {
    for(int j = 0; j < ntau ; ++j){
      double tauval = taus[j];
      if(tauval >= delta){
        for(int u = 0; u < (U+3); ++ u ){
            double intval = 0;
            
            for(int s = 0; s < ns; ++s){
              double sval = Ss(s);
              if (sval < tauval & sval >= (tauval - delta)){
                if (u < U){
                  intval = intval + Xmat(i, s)*tprs(sval, tauval, knots.col(u))*(taus(1)-taus(0));
                  
                }
                else if (u == U){
                  intval = intval + Xmat(i, s)*1*(taus(1)-taus(0));
                }
                else if (u == U+1){
                  intval = intval + Xmat(i, s)*sval*(taus(1)-taus(0));
                }
                else if (u == U+2){
                  intval = intval + Xmat(i, s)*tauval*(taus(1)-taus(0));
                }
              }
              else{
                intval = intval + 0 ;
              }
            }
            Xtil(temp, u) = intval; 
            
          }        
        
        }
            else if (tauval < delta){
             for(int u = 0; u < U; ++u){
               double intval = 0;
               intval = intval + Xmat(i, 0)*tprs(tauval, tauval, knots.col(u))*taus(1);
               Xtil(temp, u) =0; //intval; 
               
             }
             Xtil(temp, U) = 0; //Xmat(i, 0)*taus(1) ;
             Xtil(temp,U+1) = 0; //Xmat(i, 0)*tauval*taus(1);
             Xtil(temp,U+2) = 0;//Xmat(i, 0)*tauval*taus(1);
             
          } 
          //std::cout << temp << std::endl;
          
        
      temp++;
    }
  }
  

  //Xtil.cols(0, U-1) =  Xtil.cols(0, U-1)*Zmat;
  return Xtil;
}


// MCMC sampler
// [[Rcpp::export]]
Rcpp::List mcmc_sampler(arma::colvec Yvec, arma::mat Ymat,
                        arma::colvec Xvec, arma::mat Xmat, 
                        arma::colvec taus, arma::colvec Ss,
                        arma::mat Fk, arma::mat Fk1, arma::mat Psi, arma::mat Phi,
                        arma::mat Dmu, arma::mat Db, arma::mat Dalpha,
                        arma::mat Dc,
                        double i1, double i2, double a, double b,  Rcpp::Nullable<Rcpp::NumericVector> lag, int niter, arma::mat Zmat){
  
  
  // =======================================
  // Find dimensions of matrices
  // =======================================
  
  int nobs = Ymat.n_rows; // No. of curves
  int ntaus = Ymat.n_cols; // Frequency of data
  int K = Fk.n_cols; // Size of first basis function (wrt tau) [part of tensor product]
  int K1 = Fk1.n_cols; // Size of basis function for intercept term
  int U = Psi.n_cols; // Size of second basis function (wrt v) [part of tensor product]
  int Kphi  = Phi.n_cols; // Basis functions for X smoothin
  int totals = nobs*ntaus; 
  int dc = Phi.n_cols;  // ??same as Kphi
  double delta = isLag(lag);
    
  // =======================================
  // Initialise parameters and create vectors/
  // matrices to store estimates at each
  // iteration
  // =======================================

  ////std::cout << "initialise Rmat ... " << std::endl;
  arma::mat Rmat =  constructMatrix(Xvec, taus, Ss, Fk, Fk1, Psi,  nobs, delta, Zmat); 
  int dr = Rmat.n_cols; 
  ////std::cout << "Init alpha ... " << std::endl;
  arma::mat alpha(dr, niter) ; 
  alpha.col(0) = arma::randu(dr);
  ////std::cout << "Init c... " << std::endl;
  arma::mat csims(dc, niter) ; 
  csims.col(0) = arma::randu(dc);
  ////std::cout << "Init sigma ... " << std::endl;
  arma::mat sigma_e(1, niter);
  sigma_e.col(0) = 10 ; //arma::randu(1); // temporary
  arma::mat sigma_b(1, niter);
  sigma_b.col(0) = 10 ; //arma::randu(1);
  arma::mat sigma_mu(1, niter);
  sigma_mu.col(0) = 10; // arma::randu(1);
  arma::mat sigma_c(1, niter);
  sigma_c.col(0) = 10; // arma::randu(1);
  arma::mat sigma_v(1, niter);
  sigma_v.col(0) = 10; // arma::randu(1);
  
  // =======================================
  // Define matrices and objects
  // used in the Gibbs sampler
  // =======================================
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
  
  
  // =======================================
  // Start of the Gibbs sampler
  // =======================================
  
  for (int iter = 1; iter < niter; ++iter){
    //std::cout << "Iteration "<< iter << std::endl;
    
    
    
    //  ----- Sample the smooth (x(v)) first
    //std::cout << "Sampling c ... " << std::endl;
    //std::cout << "Qc ... " << std::endl;
    //this takes long
    // Posterior covariance for [c|. ]
    Qc = Dc*(1/as_scalar(sigma_c.col(iter - 1))) + (1/as_scalar(sigma_v.col(iter - 1)))*Phi.t()*Phi ;  //+ 0.0000001*arma::eye(dc,dc);;
    //  //std::cout << "ell... " << std::endl;
    ell = 1/(as_scalar(sigma_v.col(iter - 1)))*Phi.t()*Xvec;
    // //std::cout << "Sample mvn ... " << std::endl;
    // Sample using the Cholesky decomposition
    c_sample = sampleMVN(ell, Qc);
    Xsample = Phi*c_sample;
    // Update R matrix using the new X(v) smooth samples
    //std::cout << "Constructing new R ..." << std::endl;
    Rmat  = constructMatrix(Xsample, taus, Ss, Fk, Fk1, Psi, nobs, delta, Zmat);
    
    
    
    
  
    
    // ----- Sampling the regression coefficients
    // Create matrix of sigma_mu, sigma_b to simultaneously sample 
    // the intercept and regression surface terms
    // alpha = [mu^T B^T]
    
    
    
    sigs = join_vert(ones(K1, 1)*(1/as_scalar(sigma_mu.col(iter - 1))),
                     ones(U*K, 1)*(1/as_scalar(sigma_b.col(iter -1))));
    //std::cout << "Sampling alphas ... " << std::endl;
    //std::cout << "omegaf ... " << std::endl;
    Omegaf = repmat(sigs,1, Dalpha.n_cols)%Dalpha; // Copy the vector Dalpha.n_cols = K1 + UK times 
    // Posterior covariance of [alpha|.  ]
    //std::cout << " sigmaalphas ... " << std::endl;
    SigmaAlph = Omegaf + ((1/(as_scalar(sigma_e.col(iter - 1))))*Rmat.t()*Rmat); //+ 0.0000001*arma::eye(Rmat.n_cols, Rmat.n_cols);
    // Check if symmetric
    ////std::cout << Omegaf.is_symmetric() << std::endl;
    ////std::cout << ((1/(as_scalar(sigma_e.col(iter - 1))))*Rmat.t()*Rmat).is_symmetric() << std::endl;
    
    // NB : This is NOT the posterior mean; the L part in the cholesky stuff
    //std::cout << "mu  alphas ... " << std::endl;
    muAlpha = (1/(as_scalar(sigma_e.col(iter - 1))))*Rmat.t()*Yvec;
    alpha_sample = sampleMVN(muAlpha, SigmaAlph);
    ////std::cout << alpha_sample << std::endl;
    
    
    
    
    
    
    // ----- Generate the variance parameters
    
    
    
    
    // ===== VARIANCE PARAMETERS FOR REGRESSION PART
    
    
    
    
    // ----- 
    ////std::cout << "Sampling sigma  ... " << std::endl;
    ////std::cout << "shape e ... " << std::endl;
    shape_e = ntaus*nobs/2 + i1;
    ////std::cout << "rate e ... " << std::endl;
    rate_e = i2 + 0.5*as_scalar((Yvec - (Rmat*alpha_sample)).t()*(Yvec - (Rmat*alpha_sample)));
    sigma_e_sample = 1/arma::randg<double>( distr_param(shape_e,1/rate_e));
    // ----- 
    
    
    
    // ----- 
    //(Rcpp::rgamma(1, shape_e, 1/rate_e));
    ////std::cout << "shape mu ... " << std::endl;
    shape_mu = K/2 + a;
    ////std::cout << "rate mu ... " << std::endl;
    rate_mu = b + 0.5*as_scalar((alpha_sample.rows(0, K1-1).t()*Dmu*alpha_sample.rows(0, K1-1)));
    sigma_mu_sample =  1/arma::randg<double>(distr_param(shape_mu,1/rate_mu));
    // ----- 
    
    
    // ----- 
    ////std::cout << "shape b ... " << std::endl;
    shape_b = U*K/2 + a;
    ////std::cout << "rate b ... " << std::endl;
    rate_b = b + 0.5*as_scalar((alpha_sample.rows(K1, U*K + K1 - 1).t()*Db*alpha_sample.rows(K1, U*K + K1 - 1)));
    sigma_b_sample =  1/arma::randg<double>( distr_param(shape_b,1/rate_b));
    // ----- 
    
    
    
    
    
    // ===== VARIANCE PARAMETERS FOR X SMOOTHING PART
    
    
    // ----- 
    // [eta]
    ////std::cout << "shape v ... " << std::endl;
    shape_v =  ntaus*nobs/2 + i1;
    //std::cout << "rate v ... " << std::endl;
    rate_v = i2 + 0.5*as_scalar(((Xvec - Xsample).t()*(Xvec - Xsample)));
    sigma_v_sample = 1/arma::randg<double>( distr_param(shape_v,1/rate_v));
    //std::cout << "shape c ... " << std::endl;
    // ----- 
    
    
    // ----- 
    shape_c = (Kphi)/2 + a;
    //std::cout << "rate c ... " << std::endl;
    rate_c = b + 0.5*as_scalar((c_sample.t()*Dc*c_sample));
    sigma_c_sample = 1/arma::randg<double>( distr_param(shape_c,1/rate_c));
    // ----- 
  
    
    //  -------- Store values
    alpha.col(iter) = alpha_sample;
    csims.col(iter) = c_sample;
    sigma_e.col(iter) = sigma_e_sample;
    sigma_mu.col(iter) = sigma_mu_sample;
    sigma_b.col(iter) = sigma_b_sample;
    sigma_v.col(iter) = sigma_v_sample;
    sigma_c.col(iter) = sigma_c_sample;
    
  }
  
  
  // =======================================
  // End of the Gibbs sampler
  // =======================================
  
  return Rcpp::List::create(Rcpp::Named("alphas")=alpha,
                            Rcpp::Named("c")=csims,
                            Rcpp::Named("sigmae") = sigma_e,
                            Rcpp::Named("sigmamu") = sigma_mu,
                            Rcpp::Named("sigmab") = sigma_b,
                            Rcpp::Named("sigmac") = sigma_c,
                            Rcpp::Named("sigmav") = sigma_v);
}





// MCMC sampler 
// reparameterisation to obtain diagonal penalty matrices
// for faster computation
// [[Rcpp::export]]
Rcpp::List mcmc_sampler2(arma::colvec Yvec, arma::mat Ymat,
                         arma::colvec Xvec, arma::mat Xmat, 
                         arma::colvec taus, arma::colvec Ss,
                         arma::mat Fk, arma::mat Fk1, arma::mat Psi, arma::mat Phi,
                         arma::mat Dmu, arma::mat Db, arma::mat Dalpha,
                         arma::mat Dc, 
                         double i1, double i2, double a, double b,  Rcpp::Nullable<Rcpp::NumericVector> lag, int niter, arma::mat Zmat){
  
  // =======================================
  // Find dimensions
  int nobs = Ymat.n_rows; // No. of curves
  int ntaus = Ymat.n_cols; // Frequency of data
  int K = Fk.n_cols; // Size of first basis function (wrt tau) [part of tensor product]
  int K1 = Fk1.n_cols; // Size of basis function for intercept term
  int U = Psi.n_cols; // Size of second basis function (wrt v) [part of tensor product]
  int Kphi  = Phi.n_cols; // Basis functions for X smoothin
  int totals = nobs*ntaus; 
  int dc = Phi.n_cols;  // must remove
  double delta = isLag(lag);
  
  // Initialise parameters
  
  ////std::cout << "init Rmat ... " << std::endl;
  arma::mat Rmat = constructMatrix(Xvec, taus, Ss, Fk, Fk1, Psi,  nobs, delta, Zmat); 
  int dr = Rmat.n_cols; 
  ////std::cout << "Init alpha ... " << std::endl;
  arma::mat alpha(dr, niter) ; 
  alpha.col(0) = arma::randu(dr);
  ////std::cout << "Init c... " << std::endl;
  arma::mat csims(dc, niter) ; 
  csims.col(0) = arma::randu(dc);
  ////std::cout << "Init sigma ... " << std::endl;
  
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
    //std::cout << "Iteration "<< iter << std::endl;
    // Sample c for smooth X
    //std::cout << "Sampling c ... " << std::endl;
    //std::cout << "Qc ... " << std::endl;
    
    //this takes long
    Qc = (1/as_scalar(sigma_c.col(iter - 1))) + (1/as_scalar(sigma_v.col(iter - 1)))*Phi.t()*Phi+ 0.0000001*arma::eye(dc,dc);;
    //  //std::cout << "ell... " << std::endl;
    
    ell = 1/(as_scalar(sigma_v.col(iter - 1)))*Phi.t()*Xvec;
    // //std::cout << "Sample mvn ... " << std::endl;
    
    c_sample = sampleMVN(ell, Qc);
    Xsample = Phi*c_sample;
    //std::cout << "Constructing new R ..." << std::endl;
    Rmat  = constructMatrix(Xsample, taus, Ss, Fk, Fk1, Psi,nobs, delta, Zmat);
    
    // Create vector sigmas
    sigs = join_vert(ones(K1, 1)*(1/as_scalar(sigma_mu.col(iter - 1))),
                     ones(U*K, 1)*(1/as_scalar(sigma_b.col(iter -1))));
    //std::cout << "Sampling alphas ... " << std::endl;
    //std::cout << "omegaf ... " << std::endl;
    
    Omegaf = repmat(sigs,1, Dalpha.n_cols)%Dalpha; 
    //std::cout << " sigmaalphas ... " << std::endl;
    
    SigmaAlph = Omegaf + ((1/(as_scalar(sigma_e.col(iter - 1))))*Rmat.t()*Rmat) + 0.0000001*arma::eye(Rmat.n_cols, Rmat.n_cols);
    //std::cout << "mu  alphas ... " << std::endl;
    //std::cout << Omegaf.is_symmetric() << std::endl;
    //std::cout << ((1/(as_scalar(sigma_e.col(iter - 1))))*Rmat.t()*Rmat).is_symmetric() << std::endl;
    
    muAlpha = (1/(as_scalar(sigma_e.col(iter - 1))))*Rmat.t()*Yvec;
    alpha_sample = sampleMVN(muAlpha, SigmaAlph);
    //std::cout << alpha_sample << std::endl;
    // Generate the sigmas
    ////std::cout << "Sampling sigma  ... " << std::endl;
    ////std::cout << "shape e ... " << std::endl;
    shape_e = ntaus*nobs/2 + i1;
    
    ////std::cout << "rate e ... " << std::endl;
    rate_e = i2 + 0.5*as_scalar((Yvec - (Rmat*alpha_sample)).t()*(Yvec - (Rmat*alpha_sample)));
    sigma_e_sample = 1/arma::randg<double>( distr_param(shape_e,1/rate_e));
    //(Rcpp::rgamma(1, shape_e, 1/rate_e));
    ////std::cout << "shape mu ... " << std::endl;
    shape_mu = K/2 + a;
    
    ////std::cout << "rate mu ... " << std::endl;
    rate_mu = b + 0.5*as_scalar((alpha_sample.rows(0, K1-1).t()*Dmu*alpha_sample.rows(0, K1-1)));
    sigma_mu_sample =  1/arma::randg<double>(distr_param(shape_mu,1/rate_mu));
    
    //std::cout << "shape b ... " << std::endl;
    shape_b = U*K/2 + a;
    //std::cout << "rate b ... " << std::endl;
    
    rate_b = b + 0.5*as_scalar((alpha_sample.rows(K1, U*K + K1 - 1).t()*Db*alpha_sample.rows(K1, U*K + K1 - 1)));
    sigma_b_sample =  1/arma::randg<double>( distr_param(shape_b,1/rate_b));
    //std::cout << "shape v ... " << std::endl;
    
    shape_v =  ntaus*nobs/2 + i1;
    //std::cout << "rate v ... " << std::endl;
    
    rate_v = i2 + 0.5*as_scalar(((Xvec - Xsample).t()*(Xvec - Xsample)));
    sigma_v_sample = 1/arma::randg<double>( distr_param(shape_v,1/rate_v));
    //std::cout << "shape c ... " << std::endl;
    
    shape_c = (Kphi - nobs)/2 + a;
    //std::cout << "rate c ... " << std::endl;
    
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
                          arma::mat Psi,  int nobs,  double delta){
  // Evaluate no. of observations
  int ntau = taus.n_elem; //Xmat.n_cols;
  int ns = Ss.n_elem; //Xmat.n_cols;
  //double delta = is_null(lag);
  
  //int nobs = Xmat.n_rows;
  int totals = ntau*nobs;
  int U = Psi.n_cols;
  int K = Fk.n_cols;
  // Initialise matrix
  arma::mat Xmat =  reshape(Xvec, ns, nobs).t();
  //data_matrix(Xvec, ntau, nobs);
  arma::mat Xtil(totals, U);
  
  int temp = 0;
  double sval ;
  double tauval;
  double intval;
  int v; 
  //std::cout << Xtil.n_cols << std::endl;
  //std::cout << Xtil.n_rows << std::endl;
  
  for (int i = 0; i < nobs; ++i){
    ////std::cout << temp << std::endl;
    for (int j = 0; j < ntau; ++j){
      //double xij = Xmat(i, j) ;
      tauval= taus(j);
      //std::cout << tauval << std::endl;
      
      for (int u = 0; u < U; ++u){
        if(tauval >= delta){
        intval = 0;
        v = 0;
        for (int s = 0; s < ns; ++s){
          //std::cout << s << std::endl;
          
          sval = Ss(s);
          if ( sval <(tauval) & sval >= (tauval - delta)){
            intval = intval + Xmat(i,s)*Psi(v,u); // X_i(s)*Psi_u(v=s-t)
            ++v ; 
            //std::cout << Xmat(i,s) << std::endl;
           // std::cout << Psi(v,u) << std::endl;
            
          }else{
            intval = intval + 0;
          }
        }
        //std::cout <<  temp << std::endl;
        Xtil(temp,u) = intval*(taus(1)-taus(0));
        //std::cout <<  Xtil(temp,u) << std::endl;
        
        }else{Xtil(temp,u) = 0 ;}
      }
      temp++;
      
    }
  }
  
  arma::mat tau_mat = repmat(taus, nobs ,1);
  ////std::cout << "in construct matrix ... " << std::endl;
  
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
  
  arma::mat Rmat  =  XX;
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
                         arma::mat Dc, 
                         double i1, double i2, double a, double b, Rcpp::Nullable<Rcpp::NumericVector> lag, int niter){
  
  
  int nobs = Ymat.n_rows;
  int ntaus = Ymat.n_cols;
  // Find dimensions of basis matrices
  int K = Fk.n_cols;
  int K1 = Fk1.n_cols;
  int U = Psi.n_cols;
  int Kphi  = Phi.n_cols;
  int totals = nobs*ntaus;
  int dc = Phi.n_cols;
  double delta = isLag(lag);
  
  
  // Initialise parameters
  
  //std::cout << "init Rmat ... " << std::endl;
  // Construct the R matrix from initials
  arma::mat Rmat = constructMatrix2(Xvec, taus, Ss, Fk, Fk1, Psi,
                                    nobs, delta); 
  int dr = Rmat.n_cols; 
  ////std::cout << "Init alpha ... " << std::endl;
  arma::mat alpha(dr, niter) ; 
  alpha.col(0) = arma::randu(dr);
  ////std::cout << "Init c... " << std::endl;
  
  arma::mat csims(dc, niter) ; 
  csims.col(0) = arma::randu(dc);
  ////std::cout << "Init sigma ... " << std::endl;
  
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
    //std::cout << "Iteration "<< iter << std::endl;
    // Sample c for smooth X
    //std::cout << "Sampling c ... " << std::endl;
    //std::cout << "Qc ... " << std::endl;
    
    //this takes long
    Qc = Dc*(1/as_scalar(sigma_c.col(iter - 1))) + (1/as_scalar(sigma_v.col(iter - 1)))*Phi.t()*Phi+ 0.0000001*arma::eye(dc,dc);;
    //  //std::cout << "ell... " << std::endl;
    
    ell = 1/(as_scalar(sigma_v.col(iter - 1)))*Phi.t()*Xvec;
    // //std::cout << "Sample mvn ... " << std::endl;
    
    c_sample = sampleMVN(ell, Qc);
    Xsample = Phi*c_sample;
    //std::cout << "Constructing new R ..." << std::endl;
    Rmat  = constructMatrix2(Xsample, taus, Ss, Fk, Fk1, Psi,nobs, delta);
    
    // Create vector sigmas
    sigs = join_vert(ones(K1, 1)*(1/as_scalar(sigma_mu.col(iter - 1))),
                     ones(U*K, 1)*(1/as_scalar(sigma_b.col(iter -1))));
    //std::cout << "Sampling alphas ... " << std::endl;
    //std::cout << "omegaf ... " << std::endl;
    
    Omegaf = repmat(sigs,1, Dalpha.n_cols)%Dalpha; 
    //std::cout << " sigmaalphas ... " << std::endl;
    
    SigmaAlph = Omegaf + ((1/(as_scalar(sigma_e.col(iter - 1))))*Rmat.t()*Rmat) + 0.0000001*arma::eye(Rmat.n_cols, Rmat.n_cols);
    //std::cout << "mu  alphas ... " << std::endl;
    //std::cout << Omegaf.is_symmetric() << std::endl;
    //std::cout << ((1/(as_scalar(sigma_e.col(iter - 1))))*Rmat.t()*Rmat).is_symmetric() << std::endl;
    
    muAlpha = (1/(as_scalar(sigma_e.col(iter - 1))))*Rmat.t()*Yvec;
    alpha_sample = sampleMVN(muAlpha, SigmaAlph);
    //std::cout << alpha_sample << std::endl;
    // Generate the sigmas
    ////std::cout << "Sampling sigma  ... " << std::endl;
    ////std::cout << "shape e ... " << std::endl;
    shape_e = ntaus*nobs/2 + i1;
    
    ////std::cout << "rate e ... " << std::endl;
    rate_e = i2 + 0.5*as_scalar((Yvec - (Rmat*alpha_sample)).t()*(Yvec - (Rmat*alpha_sample)));
    sigma_e_sample = 1/arma::randg<double>( distr_param(shape_e,1/rate_e));
    //(Rcpp::rgamma(1, shape_e, 1/rate_e));
    ////std::cout << "shape mu ... " << std::endl;
    shape_mu = K/2 + a;
    
    ////std::cout << "rate mu ... " << std::endl;
    rate_mu = b + 0.5*as_scalar((alpha_sample.rows(0, K1-1).t()*Dmu*alpha_sample.rows(0, K1-1)));
    sigma_mu_sample =  1/arma::randg<double>(distr_param(shape_mu,1/rate_mu));
    
    //std::cout << "shape b ... " << std::endl;
    shape_b = U*K/2 + a;
    //std::cout << "rate b ... " << std::endl;
    
    rate_b = b + 0.5*as_scalar((alpha_sample.rows(K1, U*K + K1 - 1).t()*Db*alpha_sample.rows(K1, U*K + K1 - 1)));
    sigma_b_sample =  1/arma::randg<double>( distr_param(shape_b,1/rate_b));
    //std::cout << "shape v ... " << std::endl;
    
    shape_v =  ntaus*nobs/2 + i1;
    //std::cout << "rate v ... " << std::endl;
    
    rate_v = i2 + 0.5*as_scalar(((Xvec - Xsample).t()*(Xvec - Xsample)));
    sigma_v_sample = 1/arma::randg<double>( distr_param(shape_v,1/rate_v));
    //std::cout << "shape c ... " << std::endl;
    
    shape_c = (Kphi - nobs)/2 + a;
    //std::cout << "rate c ... " << std::endl;
    
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
// fpca
// for faster computation
// [[Rcpp::export]]
Rcpp::List mcmc_sampler4(arma::colvec Yvec, arma::mat Ymat,
                         arma::colvec Xvec, arma::mat Xmat,  arma::mat Xmat_centered, 
                         arma::colvec taus, arma::colvec Ss,
                         arma::mat Fk, arma::mat Fk1, arma::mat Psi, 
                         arma::mat fpcs, arma::mat mux, arma::mat eigs, arma::mat lambda_start, // chnage here
                         arma::mat Dmu, arma::mat Db, arma::mat Dalpha,
                         arma::mat Dc, 
                         double i1, double i2, double a, double b, Rcpp::Nullable<Rcpp::NumericVector> lag, int niter, arma::mat Zmat){
  
  // =======================================
  // Find dimensions
  int nobs = Ymat.n_rows; // No. of curves
  int ntaus = Ymat.n_cols; // Frequency of data
  int K = Fk.n_cols; // Size of first basis function (wrt tau) [part of tensor product]
  int K1 = Fk1.n_cols; // Size of basis function for intercept term
  int U = Psi.n_cols; // Size of second basis function (wrt v) [part of tensor product]
  int totals = nobs*ntaus; 
  int npc = fpcs.n_cols;
  double delta = isLag(lag);
  //arma::mat Rmat2 = constructTPRS(X2vec, taus, Ss, nobs, knots, Zmat, delta); 
  arma::mat FF1  = repmat(Fk1, nobs, 1);
  
  // Initialise parameters
  
  ////std::cout << "init Rmat ... " << std::endl;
  
  arma::mat Rmat1 = constructMatrix(Xvec, taus, Ss, Fk, Fk1, Psi,  nobs, delta, Zmat); 
  arma::mat Rmat = arma::join_horiz(FF1, Rmat1);
  
  int dr = Rmat.n_cols; 
  ////std::cout << "Init alpha ... " << std::endl;
  arma::mat alpha(dr, niter) ; 
  alpha.col(0) = arma::randu(dr);
  ////std::cout << "Init c... " << std::endl;
  arma::mat xisims(nobs*npc, niter) ; 
  xisims.col(0) = arma::randu(nobs*npc);
  
  ////std::cout << "Init sigma ... " << std::endl;
  
  arma::mat sigma_e(1, niter);
  sigma_e.col(0) = 10 ; //arma::randu(1);
  
  arma::mat sigma_b(1, niter);
  sigma_b.col(0) = 10 ; //arma::randu(1);
  
  arma::mat sigma_mu(1, niter);
  sigma_mu.col(0) = 10; // arma::randu(1);
  
  arma::mat lambda(npc, niter);
  lambda.col(0) = lambda_start;
  
  arma::mat sigma_v(1, niter);
  sigma_v.col(0) = arma::randu(1);
  
  // Define matrices used below
  //arma::vec ell(Kphi); 
  arma::mat xi_sample(npc*nobs, 1);
  arma::mat xi_sample_sq(nobs, npc);
    
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
  double shape_lam;
  double rate_lam;
  double shape_v;
  double rate_v;
  double sigma_v_sample;
  
  arma::mat lambda_sample(npc, 1);
  arma::mat  post_xi_var(npc, 1);
  arma::mat  post_xi_mean(npc, 1);
  
  arma::mat Fmat  = kron(arma::eye(nobs,nobs), fpcs);
  
  for (int iter = 1; iter < niter; ++iter){
    //std::cout << "Iteration "<< iter << std::endl;
    // Sample c for smooth X
    //std::cout << "Sampling c ... " << std::endl;
    
    // ====== Using fpca
    int temp = 0;
    
    //std::cout << "post c varr ... " << std::endl;
    
    post_xi_var = (lambda.col(iter-1)*as_scalar(sigma_v.col(iter - 1)))/(lambda.col(iter-1)+as_scalar(sigma_v.col(iter - 1)));
    ////std::cout << post_xi_var << std::endl;
    
    for(int x = 0; x < nobs; ++x){
     // //std::cout << "post c mean ... " << std::endl;
      post_xi_mean  = post_xi_var%(fpcs.t()*Xmat_centered.row(x).t()*(1/as_scalar(sigma_v.col(iter - 1))));
     // //std::cout << "sample ci " << std::endl;
     ////std::cout << post_xi_mean << std::endl;
       for(int k = 0 ; k < npc; ++k){
         xi_sample.row(temp) = post_xi_mean.row(k) + sqrt(post_xi_var.row(k))*arma::randn<double>();
         xi_sample_sq.col(k).row(x)= pow(xi_sample.row(temp), 2);
         temp = temp + 1;
       }
     
    }
    
    Xsample = mux +  Fmat*xi_sample;
      
  
    // //std::cout << c_sample << std::endl;
    
    // ========
    
   // //std::cout << "Constructing new R ..." << std::endl;
    Rmat  =  arma::join_horiz(FF1, constructMatrix(Xsample, taus, Ss, Fk, Fk1, Psi,nobs, delta, Zmat));
    ////std::cout << Rmat << std::endl;
    // Create vector sigmas
    sigs = join_vert(ones(K1, 1)*(1/as_scalar(sigma_mu.col(iter - 1))),
                     ones(U*K, 1)*(1/as_scalar(sigma_b.col(iter -1))));
    ////std::cout << "Sampling alphas ... " << std::endl;
    ////std::cout << "omegaf ... " << std::endl;

    Omegaf = repmat(sigs,1, Dalpha.n_cols)%Dalpha;
    ////std::cout << " sigmaalphas ... " << std::endl;

    SigmaAlph = Omegaf + ((1/(as_scalar(sigma_e.col(iter - 1))))*Rmat.t()*Rmat) ;//+ 0.0000001*arma::eye(Rmat.n_cols, Rmat.n_cols);
    ////std::cout << "mu  alphas ... " << std::endl;
    //std::cout << Omegaf.is_symmetric() << std::endl;
    //std::cout << ((1/(as_scalar(sigma_e.col(iter - 1))))*Rmat.t()*Rmat).is_symmetric() << std::endl;

    muAlpha = (1/(as_scalar(sigma_e.col(iter - 1))))*Rmat.t()*Yvec;
    alpha_sample = sampleMVN(muAlpha, SigmaAlph);
    ////std::cout << alpha_sample << std::endl;
    // Generate the sigmas
    ////std::cout << "Sampling sigma  ... " << std::endl;
    ////std::cout << "shape e ... " << std::endl;
    shape_e = ntaus*nobs/2 + i1;

    ////std::cout << "rate e ... " << std::endl;
    rate_e = i2 + 0.5*as_scalar((Yvec - (Rmat*alpha_sample)).t()*(Yvec - (Rmat*alpha_sample)));
    sigma_e_sample = 1/arma::randg<double>( distr_param(shape_e,1/rate_e));
    //(Rcpp::rgamma(1, shape_e, 1/rate_e));
    ////std::cout << "shape mu ... " << std::endl;
    shape_mu = K1/2 + a;

    ////std::cout << "rate mu ... " << std::endl;
    rate_mu = b + 0.5*as_scalar((alpha_sample.rows(0, K1-1).t()*Dmu*alpha_sample.rows(0, K1-1)));
    sigma_mu_sample =  1/arma::randg<double>(distr_param(shape_mu,1/rate_mu));

    //std::cout << "shape b ... " << std::endl;
    shape_b = U*K/2 + a;
    //std::cout << "rate b ... " << std::endl;

    rate_b = b + 0.5*as_scalar((alpha_sample.rows(K1, U*K + K1 - 1).t()*Db*alpha_sample.rows(K1, U*K + K1 - 1)));
    sigma_b_sample =  1/arma::randg<double>( distr_param(shape_b,1/rate_b));
    //std::cout << "shape v ... " << std::endl;

    shape_v =  ntaus*nobs/2 + i1;
    //std::cout << "rate v ... " << std::endl;

    rate_v = i2 + 0.5*as_scalar(((Xvec - Xsample).t()*(Xvec - Xsample)));
    sigma_v_sample = 1/arma::randg<double>( distr_param(shape_v,1/rate_v));
    //std::cout << "shape c ... " << std::endl;

    
    shape_lam = nobs/2 + a;
    for(int k = 0 ; k < npc ; ++k){
      rate_lam = b + 0.5*sum(xi_sample_sq.col(k));
      lambda_sample.row(k) = 1/arma::randg<double>( distr_param(shape_lam,1/rate_lam));
    }

    // Store values
    alpha.col(iter) = alpha_sample;
    xisims.col(iter) = xi_sample;

    sigma_e.col(iter) = sigma_e_sample;
    sigma_mu.col(iter) = sigma_mu_sample;
    sigma_b.col(iter) = sigma_b_sample;
    sigma_v.col(iter) = sigma_v_sample;
    lambda.col(iter) = lambda_sample;

  }
  
  return Rcpp::List::create(Rcpp::Named("alphas")=alpha,
                            Rcpp::Named("xi_k")=xisims,
                            Rcpp::Named("Xv")=Xsample,
                            Rcpp::Named("sigmae") = sigma_e,
                            Rcpp::Named("sigmamu") = sigma_mu,
                            Rcpp::Named("sigmab") = sigma_b,
                            Rcpp::Named("lambda") = lambda,
                            Rcpp::Named("sigmav") = sigma_v);
}


// MCMC sampler 
// fpca with horsehoe 
// [[Rcpp::export]]
Rcpp::List mcmc_sampler5(arma::colvec Yvec, arma::mat Ymat,
                         arma::colvec Xvec, arma::mat Xmat,  arma::mat Xmat_centered, 
                         arma::colvec taus, arma::colvec Ss,
                         arma::mat Fk, arma::mat Fk1, arma::mat Psi, 
                         arma::mat fpcs, arma::mat mux, arma::mat eigs, arma::mat lambda_start, // chnage here
                         arma::mat Dmu, arma::mat Db, arma::mat Dalpha,
                         arma::mat Dc, arma::mat Zmat,
                         double i1, double i2, double a, double b, double A,  Rcpp::Nullable<Rcpp::NumericVector> lag, int niter){
  
  // =======================================
  // Find dimensions
  int nobs = Ymat.n_rows; // No. of curves
  int ntaus = Ymat.n_cols; // Frequency of data
  int K = Fk.n_cols; // Size of first basis function (wrt tau) [part of tensor product]
  int K1 = Fk1.n_cols; // Size of basis function for intercept term
  int U = Psi.n_cols; // Size of second basis function (wrt v) [part of tensor product]
  int totals = nobs*ntaus; 
  int npc = fpcs.n_cols;
  double delta = isLag(lag);
  
  // Initialise parameters
  
  ////std::cout << "init Rmat ... " << std::endl;
  arma::mat Rmat = constructMatrix(Xvec, taus, Ss, Fk, Fk1, Psi,  nobs, delta, Zmat); 
  int dr = Rmat.n_cols; 
  ////std::cout << "Init alpha ... " << std::endl;
  arma::mat alpha(dr, niter) ; 
  alpha.col(0) = arma::randu(dr);
  ////std::cout << "Init c... " << std::endl;
  arma::mat xisims(nobs*npc, niter) ; 
  xisims.col(0) = arma::randu(nobs*npc);
  
  ////std::cout << "Init sigma ... " << std::endl;
  
  arma::mat sigma_e(1, niter);
  sigma_e.col(0) = 10 ; //arma::randu(1);
  
  arma::mat sigma_b(1, niter);
  sigma_b.col(0) = 10 ; //arma::randu(1);
  
  arma::mat sigma_mu(1, niter);
  sigma_mu.col(0) = 10; // arma::randu(1);
  
  arma::mat lambda(npc, niter);
  lambda.col(0) = lambda_start;
  
  arma::mat vk(npc, niter);
  vk.col(0) = arma::randu(npc);
  
  arma::mat kappa2(1, niter);
  kappa2.col(0) = arma::randu(1);
  
  arma::mat sigma_v(1, niter);
  sigma_v.col(0) = arma::randu(1);
  
  // Define matrices used below
  //arma::vec ell(Kphi); 
  arma::mat xi_sample(npc*nobs, 1);
  arma::mat xi_sample_sq(nobs, npc);
  
  arma::mat asims(1, niter);
  asims.col(0) = randu(1);
  
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
  double shape_lam;
  double rate_lam;
  double shape_v;
  double rate_v;
  double sigma_v_sample;
  double kappa2_sample; 
  double rate_kappa2; 
  double a_sample;
  double kappa0; 
  
  arma::mat lambda_sample(npc, 1);
  arma::mat vk_sample(npc, 1);
  
  arma::mat  post_xi_var(npc, 1);
  arma::mat  post_xi_mean(npc, 1);
  
  
  arma::mat Fmat  = kron(arma::eye(nobs,nobs), fpcs);
  
  for (int iter = 1; iter < niter; ++iter){
    //std::cout << "Iteration "<< iter << std::endl;
    // Sample c for smooth X
    //std::cout << "Sampling c ... " << std::endl;
    
    // ====== Using fpca
    int temp = 0;
    
    //std::cout << "post c varr ... " << std::endl;
    
    post_xi_var = (as_scalar(kappa2.col(iter-1))*lambda.col(iter-1)*
      as_scalar(sigma_v.col(iter - 1)))/((lambda.col(iter-1)*as_scalar(kappa2.col(iter-1)))+as_scalar(sigma_v.col(iter - 1)));
    ////std::cout << post_xi_var << std::endl;
    
    for(int x = 0; x < nobs; ++x){
      // //std::cout << "post c mean ... " << std::endl;
      post_xi_mean  = post_xi_var%(fpcs.t()*Xmat_centered.row(x).t()*(1/as_scalar(sigma_v.col(iter - 1))));
      // //std::cout << "sample ci " << std::endl;
      ////std::cout << post_xi_mean << std::endl;
      for(int k = 0 ; k < npc; ++k){
        xi_sample.row(temp) = post_xi_mean.row(k) + sqrt(post_xi_var.row(k))*arma::randn<double>();
        xi_sample_sq.col(k).row(x)= pow(xi_sample.row(temp), 2);
        temp = temp + 1;
      }
      
    }
    
    Xsample = mux +  Fmat*xi_sample;
    
    
    // //std::cout << c_sample << std::endl;
    
    // ========
    
    // //std::cout << "Constructing new R ..." << std::endl;
    Rmat  = constructMatrix(Xsample, taus, Ss, Fk, Fk1, Psi,nobs,delta, Zmat);
    ////std::cout << Rmat << std::endl;
    // Create vector sigmas
    sigs = join_vert(ones(K1, 1)*(1/as_scalar(sigma_mu.col(iter - 1))),
                     ones(U*K, 1)*(1/as_scalar(sigma_b.col(iter -1))));
    ////std::cout << "Sampling alphas ... " << std::endl;
    ////std::cout << "omegaf ... " << std::endl;
    
    Omegaf = repmat(sigs,1, Dalpha.n_cols)%Dalpha;
    ////std::cout << " sigmaalphas ... " << std::endl;
    
    SigmaAlph = Omegaf + ((1/(as_scalar(sigma_e.col(iter - 1))))*Rmat.t()*Rmat); // + 0.0000001*arma::eye(Rmat.n_cols, Rmat.n_cols);
    ////std::cout << "mu  alphas ... " << std::endl;
    ////std::cout << Omegaf.is_symmetric() << std::endl;
    ////std::cout << ((1/(as_scalar(sigma_e.col(iter - 1))))*Rmat.t()*Rmat).is_symmetric() << std::endl;
    
    muAlpha = (1/(as_scalar(sigma_e.col(iter - 1))))*Rmat.t()*Yvec;
    alpha_sample = sampleMVN(muAlpha, SigmaAlph);
    ////std::cout << alpha_sample << std::endl;
    // Generate the sigmas
    ////std::cout << "Sampling sigma  ... " << std::endl;
    ////std::cout << "shape e ... " << std::endl;
    shape_e = ntaus*nobs/2 + i1;
    
    ////std::cout << "rate e ... " << std::endl;
    rate_e = i2 + 0.5*as_scalar((Yvec - (Rmat*alpha_sample)).t()*(Yvec - (Rmat*alpha_sample)));
    sigma_e_sample = 1/arma::randg<double>( distr_param(shape_e,1/rate_e));
    //(Rcpp::rgamma(1, shape_e, 1/rate_e));
    ////std::cout << "shape mu ... " << std::endl;
    shape_mu = K1/2 + a;
    
    ////std::cout << "rate mu ... " << std::endl;
    rate_mu = b + 0.5*as_scalar((alpha_sample.rows(0, K1-1).t()*Dmu*alpha_sample.rows(0, K1-1)));
    sigma_mu_sample =  1/arma::randg<double>(distr_param(shape_mu,1/rate_mu));
    
    //std::cout << "shape b ... " << std::endl;
    shape_b = U*K/2 + a;
    //std::cout << "rate b ... " << std::endl;
    
    rate_b = b + 0.5*as_scalar((alpha_sample.rows(K1, U*K + K1 - 1).t()*Db*alpha_sample.rows(K1, U*K + K1 - 1)));
    sigma_b_sample =  1/arma::randg<double>( distr_param(shape_b,1/rate_b));
    //std::cout << "shape v ... " << std::endl;
    
    shape_v =  ntaus*nobs/2 + i1;
    //std::cout << "rate v ... " << std::endl;
    
    rate_v = i2 + 0.5*as_scalar(((Xvec - Xsample).t()*(Xvec - Xsample)));
    sigma_v_sample = 1/arma::randg<double>( distr_param(shape_v,1/rate_v));
    //std::cout << "shape c ... " << std::endl;
    
    //std::cout << "horseshoe stuff" << std::endl;
    
    shape_lam = (nobs+1)/2 ;
    ////std::cout << vk_sample << std::endl;
    for(int k = 0 ; k < npc ; ++k){
      rate_lam = (1/as_scalar(vk.col(iter-1).row(k))) + (1/as_scalar(kappa2.col(iter-1)))*sum(xi_sample_sq.col(k));
      //std::cout << rate_lam << std::endl;
      
      lambda_sample.row(k) = 1/arma::randg<double>( distr_param(shape_lam,1/rate_lam));
      vk_sample.row(k)= 1/arma::randg<double>( distr_param(0.5,1/((1/as_scalar(lambda.col(iter-1).row(k))) + 1/pow(A, 2))));
    }
    
    //std::cout << "Sampling kappa2" << std::endl;
    // 0 is column, 1 is row
    kappa0 = (2.0/(1.0*npc - 2.0))*(sqrt(2/1.0*totals));
    rate_kappa2 = 0.5*sum((1/lambda_sample.t())%(sum(xi_sample_sq, 0))) + as_scalar(asims.col(iter-1)) ;
    kappa2_sample = 1/arma::randg<double>( distr_param(nobs*npc/2 + 0.5, 1/rate_kappa2));
      
    a_sample = 1/arma::randg<double>( distr_param(0.5, 1/((1/kappa2_sample)+ 1/pow(kappa0, 2) )));
      
    // Store values
    alpha.col(iter) = alpha_sample;
    xisims.col(iter) = xi_sample;
    sigma_e.col(iter) = sigma_e_sample;
    sigma_mu.col(iter) = sigma_mu_sample;
    sigma_b.col(iter) = sigma_b_sample;
    sigma_v.col(iter) = sigma_v_sample;
    lambda.col(iter) = lambda_sample;
    vk.col(iter) = vk_sample;
    kappa2.col(iter) = kappa2_sample;
    asims.col(iter) = a_sample;
  }
  
  return Rcpp::List::create(Rcpp::Named("alphas")=alpha,
                            Rcpp::Named("xi_k")=xisims,
                            Rcpp::Named("Xv")=Xsample,
                            Rcpp::Named("sigmae") = sigma_e,
                            Rcpp::Named("sigmamu") = sigma_mu,
                            Rcpp::Named("sigmab") = sigma_b,
                            Rcpp::Named("lambda") = lambda,
                            Rcpp::Named("sigmav") = sigma_v);
}


// vectorise test
// [[Rcpp::export]]
arma::mat vector_test(arma::mat A){
  arma::mat B = vectorise(A.t());
  return B;
}

// make mat
// [[Rcpp::export]]
arma::mat make_mat(arma::mat A, int nrow, int ncol){
  arma::mat B = reshape(A, nrow, ncol).t(); // fills by column
  return B;
}

// make mat
// [[Rcpp::export]]
double sum_mat(arma::mat A, arma::mat B){
  double vec = sum((1/B.t())/sum(A, 0));
  return vec;
}



// MCMC sampler without smoothing x and 2 covariates
// [[Rcpp::export]]
Rcpp::List  mcmc_sampler6(arma::colvec Yvec, arma::mat Ymat,
                         arma::colvec X1vec, arma::mat Xmat, arma::colvec X2vec, arma::mat Xmat_centered, 
                         arma::colvec taus, arma::colvec Ss,
                         arma::mat Fk, arma::mat Fk1, arma::mat Psi, 
                         arma::mat fpcs, arma::mat mux, arma::mat eigs, arma::mat lambda_start, // chnage here
                         arma::mat Dmu, arma::mat Db, arma::mat Dalpha,
                         arma::mat Dc, arma::mat knots, arma::mat Zmat,
                         double i1, double i2, double a, double b,  Rcpp::Nullable<Rcpp::NumericVector> lag, int niter){
  
  // =======================================
  // Find dimensions
  int nobs = Ymat.n_rows; // No. of curves
  int ntaus = Ymat.n_cols; // Frequency of data
  int K = Fk.n_cols; // Size of first basis function (wrt tau) [part of tensor product]
  int K1 = Fk1.n_cols; // Size of basis function for intercept term
  int U = Psi.n_cols; // Size of second basis function (wrt v) [part of tensor product]
  int totals = nobs*ntaus; 
  int npc = fpcs.n_cols;
  double delta = isLag(lag);
  
  // Initialise parameters
  ////std::cout << "init Rmat ... " << std::endl;
  arma::mat Rmat1 = constructMatrix(X1vec, taus, Ss, Fk, Fk1, Psi,  nobs, delta, Zmat); 
  //arma::mat Rmat1 = constructTPRS(X1vec, taus, Ss,  nobs, knots, Zmat, delta); 
  arma::mat Rmat2 =constructMatrix(X2vec, taus, Ss, Fk, Fk1, Psi,  nobs, delta, Zmat); 
  //arma::mat Rmat2 = constructTPRS(X2vec, taus, Ss, nobs, knots, Zmat, delta); 
  arma::mat FF1  = repmat(Fk1, nobs, 1);
  arma::mat Rmat = arma::join_horiz(FF1, Rmat1, Rmat2);
  
  int dr1 = Rmat1.n_cols; 
  int dr2 = Rmat2.n_cols; 
  int dr = Rmat.n_cols; 
  
  ////std::cout << "Init alpha ... " << std::endl;
  arma::mat alpha(dr, niter) ; 
  alpha.col(0) = arma::randu(dr);
  ////std::cout << "Init c... " << std::endl;
  arma::mat xisims(nobs*npc, niter) ; 
  xisims.col(0) = arma::randu(nobs*npc);
  
  ////std::cout << "Init sigma ... " << std::endl;
  
  arma::mat sigma_e(1, niter);
  sigma_e.col(0) = 10 ; //arma::randu(1);
  
  arma::mat sigma_b1(1, niter);
  sigma_b1.col(0) = 10 ; //arma::randu(1);
  arma::mat sigma_b2(1, niter);
  sigma_b2.col(0) = 10 ; //arma::randu(1);
  
  
  arma::mat sigma_mu(1, niter);
  sigma_mu.col(0) = 10; // arma::randu(1);
  
  arma::mat lambda(npc, niter);
  lambda.col(0) = lambda_start;
  
  arma::mat sigma_v(1, niter);
  sigma_v.col(0) = arma::randu(1);
  
  // Define matrices used below
  //arma::vec ell(Kphi); 
  arma::mat xi_sample(npc*nobs, 1);
  arma::mat xi_sample_sq(nobs, npc);
  
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
  double shape_b1;
  double rate_b1;
  double sigma_b1_sample;
  
  double shape_b2;
  double rate_b2;
  double sigma_b2_sample;
  
  double shape_lam;
  double rate_lam;
  double shape_v;
  double rate_v;
  double sigma_v_sample;
  double nknots = dr1;
  
  arma::mat lambda_sample(npc, 1);
  arma::mat  post_xi_var(npc, 1);
  arma::mat  post_xi_mean(npc, 1);
  
  arma::mat Fmat  = kron(arma::eye(nobs,nobs), fpcs);
  
  for (int iter = 1; iter < niter; ++iter){
    //std::cout << "Iteration "<< iter << std::endl;
    // Sample c for smooth X
    //std::cout << "Sampling c ... " << std::endl;
    
    // ====== Using fpca
    //int temp = 0;

    
    
    std::cout << dr1 << std::endl;
    

    // Create vector sigmas
    sigs = join_vert(ones(K1, 1)*(1/as_scalar(sigma_mu.col(iter - 1))),
                     ones(dr1, 1)*(1/as_scalar(sigma_b1.col(iter -1))), 
                     ones(dr1, 1)*(1/as_scalar(sigma_b2.col(iter -1))));
    ////std::cout << "Sampling alphas ... " << std::endl;
    ////std::cout << "omegaf ... " << std::endl;

    Omegaf = repmat(sigs,1, Dalpha.n_cols)%Dalpha;
    std::cout << " sigmaalphas ... " << std::endl;
    
    SigmaAlph = Omegaf + ((1/(as_scalar(sigma_e.col(iter - 1))))*Rmat.t()*Rmat); //+ 0.0000001*arma::eye(Rmat.n_cols, Rmat.n_cols);
    std::cout << "mu  alphas ... " << std::endl;
    //std::cout << Omegaf << std::endl;
    std::cout << Dalpha.is_symmetric() << std::endl;
    std::cout << ((1/(as_scalar(sigma_e.col(iter - 1))))*Rmat.t()*Rmat).is_symmetric() << std::endl;
    //std::cout << SigmaAlph.eigenvalues() <<std::endl;
    muAlpha = (1/(as_scalar(sigma_e.col(iter - 1))))*Rmat.t()*Yvec;
    alpha_sample = sampleMVN(muAlpha, SigmaAlph);
    //std::cout << alpha_sample << std::endl;
    // Generate the sigmas
    ////std::cout << "Sampling sigma  ... " << std::endl;
    ////std::cout << "shape e ... " << std::endl;
    shape_e = ntaus*nobs/2 + i1;
    
    ////std::cout << "rate e ... " << std::endl;
    rate_e = i2 + 0.5*as_scalar((Yvec - (Rmat*alpha_sample)).t()*(Yvec - (Rmat*alpha_sample)));
    sigma_e_sample = 1/arma::randg<double>( distr_param(shape_e,1/rate_e));
    //(Rcpp::rgamma(1, shape_e, 1/rate_e));
    ////std::cout << "shape mu ... " << std::endl;
    shape_mu = K1/2 + a;
    
    ////std::cout << "rate mu ... " << std::endl;
    rate_mu = b + 0.5*as_scalar((alpha_sample.rows(0, K1-1).t()*Dmu*alpha_sample.rows(0, K1-1)));
    sigma_mu_sample =  1/arma::randg<double>(distr_param(shape_mu,1/rate_mu));
    
    std::cout << "shape b1 ... " << std::endl;
    shape_b1 = (dr1)/2 + a;
    //std::cout << "rate b ... " << std::endl;
    
    rate_b1 = b + 0.5*as_scalar((alpha_sample.rows(K1, nknots+ K1 - 1).t()*Db*alpha_sample.rows(K1, nknots + K1 - 1)));
    sigma_b1_sample =  1/arma::randg<double>( distr_param(shape_b1,1/rate_b1));
    
    std::cout << "shape b2 ... " << std::endl;
    shape_b2 =(dr2)/2 + a;
    //std::cout << "rate b ... " << std::endl;
    
    rate_b2 = b + 0.5*as_scalar((alpha_sample.rows(nknots + K1, 2*(nknots) + K1 - 1).t()*Db*alpha_sample.rows(nknots + K1, 2*(nknots) + K1 - 1)));
    sigma_b2_sample =  1/arma::randg<double>( distr_param(shape_b2,1/rate_b2));
    //std::cout << "shape v ... " << std::endl;
   
    
    // Store values
    alpha.col(iter) = alpha_sample;
    xisims.col(iter) = xi_sample;
    
    sigma_e.col(iter) = sigma_e_sample;
    sigma_mu.col(iter) = sigma_mu_sample;
    sigma_b1.col(iter) = sigma_b1_sample;
    sigma_b2.col(iter) = sigma_b2_sample;
    
    
  }
  
  return Rcpp::List::create(Rcpp::Named("alphas")=alpha,
                            Rcpp::Named("sigmae") = sigma_e,
                            Rcpp::Named("sigmamu") = sigma_mu,
                            Rcpp::Named("sigmab1") = sigma_b1,
                            Rcpp::Named("sigmab2") = sigma_b2);
}


// Function to sample regression coefficients 
// and variance parameters using the Grouped
// Horsehoe Prior
// [[Rcpp::export]]
Rcpp::List sampleRegCoeff(arma::mat Rmat, arma::colvec Yvec, int dr1, arma::mat Dmu, int p,
                          double A, double B, double C, double a, double b,  arma::mat curr_delta_guk,
                          arma::mat curr_v_delta,  arma::mat curr_lambda, arma::mat curr_p_g, 
                          double curr_tau , double curr_a_tau,   double curr_sigma_e, double curr_sigma_mu){
  
  int dr =  Dmu.n_cols + p*dr1;
  int K1  = Dmu.n_cols;
  
  arma::mat SigmaAlph(dr, dr);
  arma::mat muAlpha(dr, 1);
  arma::mat alpha_sample(dr, 1);
  
  // arma::colvec curr_delta_guk = curr_parms["delta"];
  // arma::mat curr_v_delta = curr_parms["v_delta"];
  // 
  // 
  // arma::mat curr_lambda = curr_parms["lambda"];
  // arma::mat curr_p_g = curr_parms["p_g"];
  // 
  // 
  // double curr_tau = curr_parms["tau"];
  // double curr_a_tau = curr_parms["a_tau"];
  // 
  // 
  // double curr_sigma_e = curr_parms["sigma_e"];
  // double curr_sigma_mu = curr_parms["sigma_mu"];
  // 
  arma::mat Dmat_inv = diagmat(join_vert(ones(dr1, 1)*(1/as_scalar(curr_lambda(0,0))), 
                           ones(dr1, 1)*(1/as_scalar(curr_lambda(1,0)))))%diagmat(1/curr_delta_guk);
  
  // Initialise parameters
  arma::mat prior_cov(K1 + p*dr1, K1 + p*dr1);
  prior_cov.submat(0, 0, Dmu.n_cols - 1, Dmu.n_cols - 1 ) =  (1/curr_sigma_mu)*Dmu;
  prior_cov.submat(Dmu.n_cols,Dmu.n_cols , dr - 1, dr - 1 ) = (1/curr_tau)*Dmat_inv;
  //must global parameters + predictor param + within pred para,
  
  // ==========================================================
  // Sample parameters for alphas from MVN
  // ==========================================================
  
  SigmaAlph = prior_cov + (1/(curr_sigma_e)*Rmat.t()*Rmat); //+ 0.0000001*arma::eye(Rmat.n_cols, Rmat.n_cols);
  std::cout << "mu  alphas ... " << std::endl;
  //std::cout << Omegaf << std::endl;
 // std::cout << Dalpha.is_symmetric() << std::endl;
  //std::cout << ((1/(as_scalar(sigma_e.col(iter - 1))))*Rmat.t()*Rmat).is_symmetric() << std::endl;
  //std::cout << SigmaAlph.eigenvalues() <<std::endl;
  muAlpha = (1/(curr_sigma_e))*Rmat.t()*Yvec;
  alpha_sample = sampleMVN(muAlpha, SigmaAlph);
  
  
  // ==========================================================
  // Sample Horseshoe prior parameters
  // ==========================================================
  
  // Global scale parameter
  // double myval = 0.009;
  std::cout << "Sampling taus" << std::endl;
  std::cout << 0.5*as_scalar((alpha_sample.rows(K1, dr - 1).t()*Dmat_inv*alpha_sample.rows(K1, dr - 1))) +
    (1/curr_a_tau) << std::endl;
  double tau_rate  = 0.5*as_scalar((alpha_sample.rows(K1, dr - 1).t()*Dmat_inv*alpha_sample.rows(K1, dr - 1))) +
    (1/curr_a_tau);
  double tau_sample =  1/arma::randg<double>( distr_param(0.5*p*dr1 + 0.5,1/tau_rate));
  double a_tau_sample = 1/arma::randg<double>( distr_param(1.676,1/((1/tau_sample)+(1/(pow(A,2))))));
  
  // Local predictor level shrinkage
  arma::mat lambda(p,1);
  arma::mat p_g(p,1);
  std::cout << "Sampling lambda" << std::endl;
  arma::mat curr_delta_g;
  for(int g = 0; g < p ; ++g){
    curr_delta_g = curr_delta_guk.rows(g*dr1, (g+1)*dr1 - 1);
    std::cout << alpha_sample.rows(K1+ g*dr1 ,K1 + (g+1)*dr1 - 1).t()*diagmat(1/(tau_sample*curr_delta_g))*
      alpha_sample.rows(K1+ g*dr1 ,K1 + (g+1)*dr1 - 1) << std::endl;
    // Evaluate the rate parameter for each lambda_g
    double rate_lambdag = 0.5*as_scalar(alpha_sample.rows(K1+ g*dr1 ,K1 + (g+1)*dr1 - 1).t()*diagmat(1/(tau_sample*curr_delta_g))*
      alpha_sample.rows(K1+ g*dr1 ,K1 + (g+1)*dr1 - 1)) + as_scalar(curr_p_g(g, 0));
    
    lambda(g,0) =  1/arma::randg<double>( distr_param(0.5*dr1 + 0.5, 1/rate_lambdag));
    
    p_g(g,0) =  1/arma::randg<double>( distr_param(1, 1/((1/lambda(g,0))+(1/pow(B,2)))));
  }
  
  // Within predictor level shrinkage
  arma::mat delta_guk(p*dr1, 1);
  arma::mat v_delta(p*dr1, 1);
  double g_uk = 0;
  std::cout << "Sampling delta" << std::endl;
  
  for(int gg = 0; gg < p; ++gg){
    for(int uk = 0; uk < dr1 ; ++uk){
      
      // Evaluate rate for delta_g,ij
      double rate_delta_guk = 0.5*(1/(tau_sample*lambda(gg,0)))*pow(as_scalar(alpha_sample(K1+ gg*dr1 + uk)),2) + 
      1/as_scalar(curr_v_delta(g_uk,0));
        
      delta_guk(g_uk,0) = 1; //1/arma::randg<double>( distr_param(1.00, 1/rate_delta_guk));
      
      v_delta(g_uk,0) = 1; // 1/arma::randg<double>(distr_param(1.00, 1/((1/(delta_guk(g_uk,0)))+(1/pow(C,2)))));
      
      g_uk++;
    }
  }
  
  // ==========================================================
  // Sample variance for  mu
  // ==========================================================
  ////std::cout << "shape mu ... " << std::endl;
  ////std::cout << "shape mu ... " << std::endl;
  
  std::cout << "sampling sigma mu" << std::endl;
  
  double rate_mu = b + 0.5*as_scalar((alpha_sample.rows(0, K1-1).t()*Dmu*alpha_sample.rows(0, K1-1)));
  double sigma_mu_sample = 1/arma::randg<double>(distr_param(0.5*K1 + a,1/rate_mu));
  std::cout << "going back to mcmc7" << std::endl;
  
  return Rcpp::List::create(Rcpp::Named("alpha") = alpha_sample,
                            Rcpp::Named("sigmamu") = sigma_mu_sample,
                            Rcpp::Named("tau") = tau_sample,
                            Rcpp::Named("a_tau") = a_tau_sample,
                            Rcpp::Named("lambda") = lambda,
                            Rcpp::Named("p_g") = p_g,
                            Rcpp::Named("delta_guk") = delta_guk, 
                            Rcpp::Named("v_delta") = v_delta);
}


// MCMC sampler without smoothing x and 2 covariates, Grouped horseshoe prior
// [[Rcpp::export]]
Rcpp::List  mcmc_sampler7(arma::colvec Yvec, arma::mat Ymat,
                          arma::colvec X1vec, arma::mat Xmat, arma::colvec X2vec, arma::mat Xmat_centered, 
                          arma::colvec taus, arma::colvec Ss,
                          arma::mat Fk, arma::mat Fk1, arma::mat Psi, 
                          arma::mat fpcs, arma::mat mux, arma::mat eigs, arma::mat lambda_start, // chnage here
                          arma::mat Dmu, arma::mat Db, arma::mat Dalpha,
                          arma::mat Dc, arma::mat knots, arma::mat Zmat,
                          double i1, double i2, double a, double b, double A, double B, double C,
                          Rcpp::Nullable<Rcpp::NumericVector> lag, int niter){
  // =======================================
  // Find dimensions
  int nobs = Ymat.n_rows; // No. of curves
  int ntaus = Ymat.n_cols; // Frequency of data
  int K = Fk.n_cols; // Size of first basis function (wrt tau) [part of tensor product]
  int K1 = Fk1.n_cols; // Size of basis function for intercept term
  int U = Psi.n_cols; // Size of second basis function (wrt v) [part of tensor product]
  int totals = nobs*ntaus; 
  int npc = fpcs.n_cols;
  double delta = isLag(lag);
  
  // Initialise parameters
  ////std::cout << "init Rmat ... " << std::endl;
  arma::mat Rmat1 = constructMatrix(X1vec, taus, Ss, Fk, Fk1, Psi,  nobs, delta, Zmat); 
  //arma::mat Rmat1 = constructTPRS(X1vec, taus, Ss,  nobs, knots, Zmat, delta); 
  arma::mat Rmat2 = constructMatrix(X2vec, taus, Ss, Fk, Fk1, Psi,  nobs, delta, Zmat); 
  //arma::mat Rmat2 = constructTPRS(X2vec, taus, Ss, nobs, knots, Zmat, delta); 
  arma::mat FF1  = repmat(Fk1, nobs, 1);
  arma::mat Rmat = arma::join_horiz(FF1, Rmat1, Rmat2);
  
  int dr1 = Rmat1.n_cols; 
  int dr2 = Rmat2.n_cols; 
  int dr = Rmat.n_cols; 
  
  // Initialise matrices
  arma::mat alpha(dr, niter);
  alpha.col(0) = arma::randu(dr);
  
  arma::mat sigmae(1, niter);
  sigmae.col(0)= arma::randu(1);
  
  arma::mat sigmamu(1, niter);
  sigmamu.col(0) = arma::randu(1);
  
  arma::mat tau(1, niter);
  tau.col(0)  = arma::randu(1);
  
  arma::mat a_tau(1, niter);
  a_tau.col(0) = arma::randu(1);
  
  
  int p = 2;
  
  arma::mat lambda(p, niter);
  lambda.col(0) = arma::randu(p);
  
  arma::mat p_g(p, niter);
  p_g.col(0)  = arma::randu(p);
  
  arma::mat delta_guk(p*dr1, niter);
  delta_guk.col(0) = ones(1, p*dr1).t() ;//arma::randu(p*dr1);
  
  arma::mat v_delta(p*dr1, niter);
  v_delta.col(0) =ones(1, p*dr1).t() ; // arma::randu(p*dr1); 
  
  Rcpp::List reg_sample; 
  double shape_e;
  double rate_e;
  double sigma_e_sample;
  
  for(int iter = 1; iter < niter ; ++iter){
    
    // Sample the regression coefficients
    reg_sample = sampleRegCoeff(Rmat, Yvec, dr1, Dmu, p,
      A, B, C,  a,  b, delta_guk.col(iter - 1),
      v_delta.col(iter-1),  lambda.col(iter-1), p_g.col(iter - 1), 
      as_scalar(tau(0,iter - 1)), as_scalar(a_tau(0,iter - 1)),
      as_scalar(sigmae(0,iter - 1)), as_scalar(sigmamu(0,iter - 1)));
    
    // Store the values
    arma::mat alpha_sample = reg_sample["alpha"];
    alpha.col(iter) = alpha_sample;
    
    sigmamu(0,iter) = reg_sample["sigmamu"];
    tau(0,iter) = reg_sample["tau"];
    a_tau(0,iter) = reg_sample["a_tau"];
    
    arma::mat lambda_sample = reg_sample["lambda"];
    lambda.col(iter) = lambda_sample;
    arma::mat p_g_sample = reg_sample["p_g"];
    p_g.col(iter) = p_g_sample;
    
    arma::mat delta_guk_sample = reg_sample["delta_guk"];
    delta_guk.col(iter) = delta_guk_sample;
    
    arma::mat v_delta_sample = reg_sample["v_delta"];
    v_delta.col(iter) = v_delta_sample;

    // Sample the error terms
    shape_e = ntaus*nobs/2 + i1;
    
    ////std::cout << "rate e ... " << std::endl;
    rate_e = i2 + 0.5*as_scalar((Yvec - (Rmat*alpha.col(iter))).t()*(Yvec - (Rmat*alpha.col(iter))));
    sigmae(0,iter)= 1/arma::randg<double>( distr_param(shape_e,1/rate_e));
    
    
    
  }
  
  return Rcpp::List::create(Rcpp::Named("alpha") = alpha,
    Rcpp::Named("sigmamu") = sigmamu,
    Rcpp::Named("tau") = tau,
    Rcpp::Named("a_tau") = a_tau,
    Rcpp::Named("lambda") = lambda,
    Rcpp::Named("p_g") = p_g,
    Rcpp::Named("delta_guk") = delta_guk, 
    Rcpp::Named("v_delta") = v_delta,
    Rcpp::Named("sigmamu") = sigmae);
  
}