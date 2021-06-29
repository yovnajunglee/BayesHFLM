// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// rcpparma_outerproduct
arma::mat rcpparma_outerproduct(const arma::colvec& x);
RcppExport SEXP _BayesHFLM_rcpparma_outerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_outerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_innerproduct
double rcpparma_innerproduct(const arma::colvec& x);
RcppExport SEXP _BayesHFLM_rcpparma_innerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_innerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_bothproducts
Rcpp::List rcpparma_bothproducts(const arma::colvec& x);
RcppExport SEXP _BayesHFLM_rcpparma_bothproducts(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_bothproducts(x));
    return rcpp_result_gen;
END_RCPP
}
// myfirstfunc
arma::mat myfirstfunc();
RcppExport SEXP _BayesHFLM_myfirstfunc() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(myfirstfunc());
    return rcpp_result_gen;
END_RCPP
}
// sampleMVN
arma::mat sampleMVN(arma::vec ell, arma::mat Qc);
RcppExport SEXP _BayesHFLM_sampleMVN(SEXP ellSEXP, SEXP QcSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type ell(ellSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Qc(QcSEXP);
    rcpp_result_gen = Rcpp::wrap(sampleMVN(ell, Qc));
    return rcpp_result_gen;
END_RCPP
}
// data_matrix
arma::mat data_matrix(arma::mat mymat, int nc, int nr);
RcppExport SEXP _BayesHFLM_data_matrix(SEXP mymatSEXP, SEXP ncSEXP, SEXP nrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type mymat(mymatSEXP);
    Rcpp::traits::input_parameter< int >::type nc(ncSEXP);
    Rcpp::traits::input_parameter< int >::type nr(nrSEXP);
    rcpp_result_gen = Rcpp::wrap(data_matrix(mymat, nc, nr));
    return rcpp_result_gen;
END_RCPP
}
// constructMatrix
arma::mat constructMatrix(arma::mat Xvec, arma::colvec taus, arma::colvec Ss, arma::mat Fk, arma::mat Fk1, arma::mat Psi, double delta, int nobs);
RcppExport SEXP _BayesHFLM_constructMatrix(SEXP XvecSEXP, SEXP tausSEXP, SEXP SsSEXP, SEXP FkSEXP, SEXP Fk1SEXP, SEXP PsiSEXP, SEXP deltaSEXP, SEXP nobsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Xvec(XvecSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type taus(tausSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type Ss(SsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Fk(FkSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Fk1(Fk1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Psi(PsiSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< int >::type nobs(nobsSEXP);
    rcpp_result_gen = Rcpp::wrap(constructMatrix(Xvec, taus, Ss, Fk, Fk1, Psi, delta, nobs));
    return rcpp_result_gen;
END_RCPP
}
// mcmc_sampler
Rcpp::List mcmc_sampler(arma::colvec Yvec, arma::mat Ymat, arma::colvec Xvec, arma::mat Xmat, arma::colvec taus, arma::colvec Ss, arma::mat Fk, arma::mat Fk1, arma::mat Psi, arma::mat Phi, arma::mat Dmu, arma::mat Db, arma::mat Dalpha, arma::mat Dc, double delta, double i1, double i2, double a, double b, int niter);
RcppExport SEXP _BayesHFLM_mcmc_sampler(SEXP YvecSEXP, SEXP YmatSEXP, SEXP XvecSEXP, SEXP XmatSEXP, SEXP tausSEXP, SEXP SsSEXP, SEXP FkSEXP, SEXP Fk1SEXP, SEXP PsiSEXP, SEXP PhiSEXP, SEXP DmuSEXP, SEXP DbSEXP, SEXP DalphaSEXP, SEXP DcSEXP, SEXP deltaSEXP, SEXP i1SEXP, SEXP i2SEXP, SEXP aSEXP, SEXP bSEXP, SEXP niterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec >::type Yvec(YvecSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Ymat(YmatSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type Xvec(XvecSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Xmat(XmatSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type taus(tausSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type Ss(SsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Fk(FkSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Fk1(Fk1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Psi(PsiSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Phi(PhiSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Dmu(DmuSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Db(DbSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Dalpha(DalphaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Dc(DcSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< double >::type i1(i1SEXP);
    Rcpp::traits::input_parameter< double >::type i2(i2SEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    rcpp_result_gen = Rcpp::wrap(mcmc_sampler(Yvec, Ymat, Xvec, Xmat, taus, Ss, Fk, Fk1, Psi, Phi, Dmu, Db, Dalpha, Dc, delta, i1, i2, a, b, niter));
    return rcpp_result_gen;
END_RCPP
}
// mcmc_sampler2
Rcpp::List mcmc_sampler2(arma::colvec Yvec, arma::mat Ymat, arma::colvec Xvec, arma::mat Xmat, arma::colvec taus, arma::colvec Ss, arma::mat Fk, arma::mat Fk1, arma::mat Psi, arma::mat Phi, arma::mat Dmu, arma::mat Db, arma::mat Dalpha, arma::mat Dc, double delta, double i1, double i2, double a, double b, int niter);
RcppExport SEXP _BayesHFLM_mcmc_sampler2(SEXP YvecSEXP, SEXP YmatSEXP, SEXP XvecSEXP, SEXP XmatSEXP, SEXP tausSEXP, SEXP SsSEXP, SEXP FkSEXP, SEXP Fk1SEXP, SEXP PsiSEXP, SEXP PhiSEXP, SEXP DmuSEXP, SEXP DbSEXP, SEXP DalphaSEXP, SEXP DcSEXP, SEXP deltaSEXP, SEXP i1SEXP, SEXP i2SEXP, SEXP aSEXP, SEXP bSEXP, SEXP niterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec >::type Yvec(YvecSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Ymat(YmatSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type Xvec(XvecSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Xmat(XmatSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type taus(tausSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type Ss(SsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Fk(FkSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Fk1(Fk1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Psi(PsiSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Phi(PhiSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Dmu(DmuSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Db(DbSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Dalpha(DalphaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Dc(DcSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< double >::type i1(i1SEXP);
    Rcpp::traits::input_parameter< double >::type i2(i2SEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    rcpp_result_gen = Rcpp::wrap(mcmc_sampler2(Yvec, Ymat, Xvec, Xmat, taus, Ss, Fk, Fk1, Psi, Phi, Dmu, Db, Dalpha, Dc, delta, i1, i2, a, b, niter));
    return rcpp_result_gen;
END_RCPP
}
// constructMatrix2
arma::mat constructMatrix2(arma::mat Xvec, arma::colvec taus, arma::colvec Ss, arma::mat Fk, arma::mat Fk1, arma::mat Psi, double delta, int nobs);
RcppExport SEXP _BayesHFLM_constructMatrix2(SEXP XvecSEXP, SEXP tausSEXP, SEXP SsSEXP, SEXP FkSEXP, SEXP Fk1SEXP, SEXP PsiSEXP, SEXP deltaSEXP, SEXP nobsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Xvec(XvecSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type taus(tausSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type Ss(SsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Fk(FkSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Fk1(Fk1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Psi(PsiSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< int >::type nobs(nobsSEXP);
    rcpp_result_gen = Rcpp::wrap(constructMatrix2(Xvec, taus, Ss, Fk, Fk1, Psi, delta, nobs));
    return rcpp_result_gen;
END_RCPP
}
// mcmc_sampler3
Rcpp::List mcmc_sampler3(arma::colvec Yvec, arma::mat Ymat, arma::colvec Xvec, arma::mat Xmat, arma::colvec taus, arma::colvec Ss, arma::mat Fk, arma::mat Fk1, arma::mat Psi, arma::mat Phi, arma::mat Dmu, arma::mat Db, arma::mat Dalpha, arma::mat Dc, double delta, double i1, double i2, double a, double b, int niter);
RcppExport SEXP _BayesHFLM_mcmc_sampler3(SEXP YvecSEXP, SEXP YmatSEXP, SEXP XvecSEXP, SEXP XmatSEXP, SEXP tausSEXP, SEXP SsSEXP, SEXP FkSEXP, SEXP Fk1SEXP, SEXP PsiSEXP, SEXP PhiSEXP, SEXP DmuSEXP, SEXP DbSEXP, SEXP DalphaSEXP, SEXP DcSEXP, SEXP deltaSEXP, SEXP i1SEXP, SEXP i2SEXP, SEXP aSEXP, SEXP bSEXP, SEXP niterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec >::type Yvec(YvecSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Ymat(YmatSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type Xvec(XvecSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Xmat(XmatSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type taus(tausSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type Ss(SsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Fk(FkSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Fk1(Fk1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Psi(PsiSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Phi(PhiSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Dmu(DmuSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Db(DbSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Dalpha(DalphaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Dc(DcSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< double >::type i1(i1SEXP);
    Rcpp::traits::input_parameter< double >::type i2(i2SEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    rcpp_result_gen = Rcpp::wrap(mcmc_sampler3(Yvec, Ymat, Xvec, Xmat, taus, Ss, Fk, Fk1, Psi, Phi, Dmu, Db, Dalpha, Dc, delta, i1, i2, a, b, niter));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_BayesHFLM_rcpparma_outerproduct", (DL_FUNC) &_BayesHFLM_rcpparma_outerproduct, 1},
    {"_BayesHFLM_rcpparma_innerproduct", (DL_FUNC) &_BayesHFLM_rcpparma_innerproduct, 1},
    {"_BayesHFLM_rcpparma_bothproducts", (DL_FUNC) &_BayesHFLM_rcpparma_bothproducts, 1},
    {"_BayesHFLM_myfirstfunc", (DL_FUNC) &_BayesHFLM_myfirstfunc, 0},
    {"_BayesHFLM_sampleMVN", (DL_FUNC) &_BayesHFLM_sampleMVN, 2},
    {"_BayesHFLM_data_matrix", (DL_FUNC) &_BayesHFLM_data_matrix, 3},
    {"_BayesHFLM_constructMatrix", (DL_FUNC) &_BayesHFLM_constructMatrix, 8},
    {"_BayesHFLM_mcmc_sampler", (DL_FUNC) &_BayesHFLM_mcmc_sampler, 20},
    {"_BayesHFLM_mcmc_sampler2", (DL_FUNC) &_BayesHFLM_mcmc_sampler2, 20},
    {"_BayesHFLM_constructMatrix2", (DL_FUNC) &_BayesHFLM_constructMatrix2, 8},
    {"_BayesHFLM_mcmc_sampler3", (DL_FUNC) &_BayesHFLM_mcmc_sampler3, 20},
    {NULL, NULL, 0}
};

RcppExport void R_init_BayesHFLM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
