// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// DDC_cpp
Rcpp::List DDC_cpp(arma::mat& X, const double& tolProbCell, const double& tolProbRow, const double& tolProbReg, const double& tolProbCorr, const double& corrlim, const int& combinRule, const int& includeSelf, const int& fastDDC, const int& qdim, const int& transFun, unsigned int& k, const unsigned int& numiter, const double& precScale, const int& standType, const int& corrType, const unsigned int& nCorr, const unsigned int& nLocScale, arma::uvec& goodCols, const int& fixedCenter, const arma::vec& center);
RcppExport SEXP _cellWise_DDC_cpp(SEXP XSEXP, SEXP tolProbCellSEXP, SEXP tolProbRowSEXP, SEXP tolProbRegSEXP, SEXP tolProbCorrSEXP, SEXP corrlimSEXP, SEXP combinRuleSEXP, SEXP includeSelfSEXP, SEXP fastDDCSEXP, SEXP qdimSEXP, SEXP transFunSEXP, SEXP kSEXP, SEXP numiterSEXP, SEXP precScaleSEXP, SEXP standTypeSEXP, SEXP corrTypeSEXP, SEXP nCorrSEXP, SEXP nLocScaleSEXP, SEXP goodColsSEXP, SEXP fixedCenterSEXP, SEXP centerSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const double& >::type tolProbCell(tolProbCellSEXP);
    Rcpp::traits::input_parameter< const double& >::type tolProbRow(tolProbRowSEXP);
    Rcpp::traits::input_parameter< const double& >::type tolProbReg(tolProbRegSEXP);
    Rcpp::traits::input_parameter< const double& >::type tolProbCorr(tolProbCorrSEXP);
    Rcpp::traits::input_parameter< const double& >::type corrlim(corrlimSEXP);
    Rcpp::traits::input_parameter< const int& >::type combinRule(combinRuleSEXP);
    Rcpp::traits::input_parameter< const int& >::type includeSelf(includeSelfSEXP);
    Rcpp::traits::input_parameter< const int& >::type fastDDC(fastDDCSEXP);
    Rcpp::traits::input_parameter< const int& >::type qdim(qdimSEXP);
    Rcpp::traits::input_parameter< const int& >::type transFun(transFunSEXP);
    Rcpp::traits::input_parameter< unsigned int& >::type k(kSEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type numiter(numiterSEXP);
    Rcpp::traits::input_parameter< const double& >::type precScale(precScaleSEXP);
    Rcpp::traits::input_parameter< const int& >::type standType(standTypeSEXP);
    Rcpp::traits::input_parameter< const int& >::type corrType(corrTypeSEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type nCorr(nCorrSEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type nLocScale(nLocScaleSEXP);
    Rcpp::traits::input_parameter< arma::uvec& >::type goodCols(goodColsSEXP);
    Rcpp::traits::input_parameter< const int& >::type fixedCenter(fixedCenterSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type center(centerSEXP);
    rcpp_result_gen = Rcpp::wrap(DDC_cpp(X, tolProbCell, tolProbRow, tolProbReg, tolProbCorr, corrlim, combinRule, includeSelf, fastDDC, qdim, transFun, k, numiter, precScale, standType, corrType, nCorr, nLocScale, goodCols, fixedCenter, center));
    return rcpp_result_gen;
END_RCPP
}
// Wrap_cpp
Rcpp::List Wrap_cpp(arma::mat& X, arma::vec& loc, arma::vec& scale, double precScale);
RcppExport SEXP _cellWise_Wrap_cpp(SEXP XSEXP, SEXP locSEXP, SEXP scaleSEXP, SEXP precScaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type loc(locSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< double >::type precScale(precScaleSEXP);
    rcpp_result_gen = Rcpp::wrap(Wrap_cpp(X, loc, scale, precScale));
    return rcpp_result_gen;
END_RCPP
}
// estLocScale_cpp
Rcpp::List estLocScale_cpp(arma::mat& X, unsigned int nLocScale, int type, double precScale, const int center, const double alpha);
RcppExport SEXP _cellWise_estLocScale_cpp(SEXP XSEXP, SEXP nLocScaleSEXP, SEXP typeSEXP, SEXP precScaleSEXP, SEXP centerSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type nLocScale(nLocScaleSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< double >::type precScale(precScaleSEXP);
    Rcpp::traits::input_parameter< const int >::type center(centerSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(estLocScale_cpp(X, nLocScale, type, precScale, center, alpha));
    return rcpp_result_gen;
END_RCPP
}
// unimcd_cpp
Rcpp::List unimcd_cpp(arma::vec& y, const double alpha, const int center);
RcppExport SEXP _cellWise_unimcd_cpp(SEXP ySEXP, SEXP alphaSEXP, SEXP centerSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const int >::type center(centerSEXP);
    rcpp_result_gen = Rcpp::wrap(unimcd_cpp(y, alpha, center));
    return rcpp_result_gen;
END_RCPP
}
// findCellPath_cpp
Rcpp::List findCellPath_cpp(arma::mat& predictors, arma::vec& response, arma::vec& weights, arma::mat& Sigmai, const arma::uvec& naMask);
RcppExport SEXP _cellWise_findCellPath_cpp(SEXP predictorsSEXP, SEXP responseSEXP, SEXP weightsSEXP, SEXP SigmaiSEXP, SEXP naMaskSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type predictors(predictorsSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type response(responseSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Sigmai(SigmaiSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type naMask(naMaskSEXP);
    rcpp_result_gen = Rcpp::wrap(findCellPath_cpp(predictors, response, weights, Sigmai, naMask));
    return rcpp_result_gen;
END_RCPP
}
// allpreds_cpp
Rcpp::List allpreds_cpp(arma::mat& X, arma::mat& S, arma::vec& mu, arma::umat& W);
RcppExport SEXP _cellWise_allpreds_cpp(SEXP XSEXP, SEXP SSEXP, SEXP muSEXP, SEXP WSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type S(SSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::umat& >::type W(WSEXP);
    rcpp_result_gen = Rcpp::wrap(allpreds_cpp(X, S, mu, W));
    return rcpp_result_gen;
END_RCPP
}
// Objective_cpp
double Objective_cpp(arma::mat& X, arma::umat& W, arma::vec& mu, arma::mat& Sigma, arma::mat& Sigmai);
RcppExport SEXP _cellWise_Objective_cpp(SEXP XSEXP, SEXP WSEXP, SEXP muSEXP, SEXP SigmaSEXP, SEXP SigmaiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::umat& >::type W(WSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Sigmai(SigmaiSEXP);
    rcpp_result_gen = Rcpp::wrap(Objective_cpp(X, W, mu, Sigma, Sigmai));
    return rcpp_result_gen;
END_RCPP
}
// updateW_cpp
arma::umat updateW_cpp(const arma::mat& X, arma::umat W, const arma::vec& mu, const arma::mat& Sigma, const arma::mat& Sigmai, const arma::vec& lambda, const arma::uword& h);
RcppExport SEXP _cellWise_updateW_cpp(SEXP XSEXP, SEXP WSEXP, SEXP muSEXP, SEXP SigmaSEXP, SEXP SigmaiSEXP, SEXP lambdaSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::umat >::type W(WSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Sigmai(SigmaiSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const arma::uword& >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(updateW_cpp(X, W, mu, Sigma, Sigmai, lambda, h));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_cellWise_DDC_cpp", (DL_FUNC) &_cellWise_DDC_cpp, 21},
    {"_cellWise_Wrap_cpp", (DL_FUNC) &_cellWise_Wrap_cpp, 4},
    {"_cellWise_estLocScale_cpp", (DL_FUNC) &_cellWise_estLocScale_cpp, 6},
    {"_cellWise_unimcd_cpp", (DL_FUNC) &_cellWise_unimcd_cpp, 3},
    {"_cellWise_findCellPath_cpp", (DL_FUNC) &_cellWise_findCellPath_cpp, 5},
    {"_cellWise_allpreds_cpp", (DL_FUNC) &_cellWise_allpreds_cpp, 4},
    {"_cellWise_Objective_cpp", (DL_FUNC) &_cellWise_Objective_cpp, 5},
    {"_cellWise_updateW_cpp", (DL_FUNC) &_cellWise_updateW_cpp, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_cellWise(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
