// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// DDC_cpp
Rcpp::List DDC_cpp(arma::mat& X, const double& tolProbCell, const double& tolProbRow, const double& tolProbReg, const double& tolProbCorr, const double& corrlim, const int& combinRule, const int& rowdetect, const int& includeSelf, const int& fastDDC, const int& absCorr, const int& qdim, const int& transFun, const int& treetype, const int& searchtype, const double& radius, const double& eps, const int& bruteForce, unsigned int& k, const unsigned int& numiter, const double& precScale);
RcppExport SEXP _cellWise_DDC_cpp(SEXP XSEXP, SEXP tolProbCellSEXP, SEXP tolProbRowSEXP, SEXP tolProbRegSEXP, SEXP tolProbCorrSEXP, SEXP corrlimSEXP, SEXP combinRuleSEXP, SEXP rowdetectSEXP, SEXP includeSelfSEXP, SEXP fastDDCSEXP, SEXP absCorrSEXP, SEXP qdimSEXP, SEXP transFunSEXP, SEXP treetypeSEXP, SEXP searchtypeSEXP, SEXP radiusSEXP, SEXP epsSEXP, SEXP bruteForceSEXP, SEXP kSEXP, SEXP numiterSEXP, SEXP precScaleSEXP) {
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
    Rcpp::traits::input_parameter< const int& >::type rowdetect(rowdetectSEXP);
    Rcpp::traits::input_parameter< const int& >::type includeSelf(includeSelfSEXP);
    Rcpp::traits::input_parameter< const int& >::type fastDDC(fastDDCSEXP);
    Rcpp::traits::input_parameter< const int& >::type absCorr(absCorrSEXP);
    Rcpp::traits::input_parameter< const int& >::type qdim(qdimSEXP);
    Rcpp::traits::input_parameter< const int& >::type transFun(transFunSEXP);
    Rcpp::traits::input_parameter< const int& >::type treetype(treetypeSEXP);
    Rcpp::traits::input_parameter< const int& >::type searchtype(searchtypeSEXP);
    Rcpp::traits::input_parameter< const double& >::type radius(radiusSEXP);
    Rcpp::traits::input_parameter< const double& >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const int& >::type bruteForce(bruteForceSEXP);
    Rcpp::traits::input_parameter< unsigned int& >::type k(kSEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type numiter(numiterSEXP);
    Rcpp::traits::input_parameter< const double& >::type precScale(precScaleSEXP);
    rcpp_result_gen = Rcpp::wrap(DDC_cpp(X, tolProbCell, tolProbRow, tolProbReg, tolProbCorr, corrlim, combinRule, rowdetect, includeSelf, fastDDC, absCorr, qdim, transFun, treetype, searchtype, radius, eps, bruteForce, k, numiter, precScale));
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
Rcpp::List estLocScale_cpp(arma::mat& X, int type, double precScale);
RcppExport SEXP _cellWise_estLocScale_cpp(SEXP XSEXP, SEXP typeSEXP, SEXP precScaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< double >::type precScale(precScaleSEXP);
    rcpp_result_gen = Rcpp::wrap(estLocScale_cpp(X, type, precScale));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_cellWise_DDC_cpp", (DL_FUNC) &_cellWise_DDC_cpp, 21},
    {"_cellWise_Wrap_cpp", (DL_FUNC) &_cellWise_Wrap_cpp, 4},
    {"_cellWise_estLocScale_cpp", (DL_FUNC) &_cellWise_estLocScale_cpp, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_cellWise(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
