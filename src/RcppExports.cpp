// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// sim_swf_cpp
Rcpp::List sim_swf_cpp(Rcpp::List args);
RcppExport SEXP _polySimIBD_sim_swf_cpp(SEXP argsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type args(argsSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_swf_cpp(args));
    return rcpp_result_gen;
END_RCPP
}
// get_arg_cpp
Rcpp::List get_arg_cpp(Rcpp::List args);
RcppExport SEXP _polySimIBD_get_arg_cpp(SEXP argsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type args(argsSEXP);
    rcpp_result_gen = Rcpp::wrap(get_arg_cpp(args));
    return rcpp_result_gen;
END_RCPP
}
// subset_bvtree_cpp
Rcpp::List subset_bvtree_cpp(std::vector<int> c, std::vector<int> t, std::vector<int> m);
RcppExport SEXP _polySimIBD_subset_bvtree_cpp(SEXP cSEXP, SEXP tSEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<int> >::type c(cSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type t(tSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(subset_bvtree_cpp(c, t, m));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_polySimIBD_sim_swf_cpp", (DL_FUNC) &_polySimIBD_sim_swf_cpp, 1},
    {"_polySimIBD_get_arg_cpp", (DL_FUNC) &_polySimIBD_get_arg_cpp, 1},
    {"_polySimIBD_subset_bvtree_cpp", (DL_FUNC) &_polySimIBD_subset_bvtree_cpp, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_polySimIBD(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
