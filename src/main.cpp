
#include "main.h"
#include "probability_v10.h"

#include <chrono>

using namespace std;

//------------------------------------------------
// draw from forwards-in-time model
//Rcpp::List sim_swf_cpp(Rcpp::List args, Rcpp::List args_functions, Rcpp::List args_progress) {
Rcpp::List sim_swf_cpp(Rcpp::List args) {
  
  // start timer
  chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
  
  // end timer
  chrono_timer(t1);
  
  return Rcpp::List::create(Rcpp::Named("foo") = -9);
}

