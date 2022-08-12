
#include "misc_v9.h"

#include <vector>

//------------------------------------------------
// [[Rcpp::export]]
Rcpp::List sim_swf_cpp(Rcpp::List args);

//------------------------------------------------
// [[Rcpp::export]]
Rcpp::List get_arg_cpp(Rcpp::List args);

//------------------------------------------------
// [[Rcpp::export]]
Rcpp::List subset_bvtree_cpp(std::vector<int> c,
                             std::vector<int> t,
                             std::vector<int> m);

//------------------------------------------------
// [[Rcpp::export]]
Rcpp::List calc_between_coi_IBD_cpp(Rcpp::List args);
