// [[Rcpp::depends(RcppArmadillo)]]

// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
using namespace arma;

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//
// [[Rcpp::export]]
arma::mat rcpparma_hello_world() {
    arma::mat m1 = arma::eye<arma::mat>(3, 3);
    arma::mat m2 = arma::eye<arma::mat>(3, 3);
	                     
    return m1 + 3 * (m1 + m2);
}


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


// Construct X matrix for 
// HFLM 
// [[Rcpp::export]]
Rcpp::List constructMatrix(arma::mat Xmat, arma::colvec taus,
                          arma::colvec Ss, arma::mat Fk, 
                          arma::mat Psi, int delta){
    // Evaluate no. of observations
    int ntau = Xmat.n_cols;
    int nobs = Xmat.n_rows;
    int totals = ntau*nobs;
    int U = Psi.n_cols;
    int K = Fk.n_cols;
    // Initialise matrix
    arma::mat Xtil(totals, U);
    
    int temp = 0;
    
    for (int i = 0; i < nobs; ++i){
        for (int j = 0; j < ntau; ++j){
            double xij = Xmat(i, j) ; 
            double tauval = taus(j);
            for (int u = 0; u < U; ++u){
                double intval = 0; 
                for (int s = 0; s < ntau; ++s){
                    double psival = Psi(s, u);
                    double sval = Ss[s];
                        if (sval >= (tauval - delta) && sval <= (tauval)){
                            intval = intval + xij*psival;
                        }else{
                            intval = intval + 0;
                        }
                }
                Xtil(temp,u) = intval/ntau;
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
   
   return Rcpp::List::create(Rcpp::Named("Xt")=Xtil,
                             Rcpp::Named("FF")=FF,
                             Rcpp::Named("XX")=XX);
}
// works 
