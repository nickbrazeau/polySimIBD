
#include "main.h"
#include "probability_v10.h"
#include "Haplotype.h"
#include "Host.h"

#include <chrono>

using namespace std;

//------------------------------------------------
// draw from forwards-in-time model
Rcpp::List sim_swf_cpp(Rcpp::List args) {
  
  // start timer
  chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
  
  // extract arguments
  vector<int> pos = rcpp_to_vector_int(args["pos"]);
  int N = rcpp_to_int(args["N"]);
  double m = rcpp_to_double(args["m"]);
  double rho = rcpp_to_double(args["rho"]);
  double mean_coi = rcpp_to_double(args["mean_coi"]);
  int tlim = rcpp_to_int(args["tlim"]);
  vector<double> odd_prob = rcpp_to_vector_double(args["odd_prob"]);
  int L = int(pos.size());
  double max_pos = pos[pos.size()-1] - pos[0];
  
  // nested vectors, first over time, then individuals
  vector<vector<Host>> pop(tlim, vector<Host>(N));
  
  // draw COI for every individual in every time step
  for (int t = 0; t < tlim; ++t) {
    for (int i = 0; i < N; ++i) {
      pop[t][i].init(mean_coi, L);
    }
  }
  
  // step through time in discrete generations, starting at the second timestep
  for (int t = 1; t < tlim; ++t) {
    for (int i = 0; i < N; ++i) {
      pop[t][i].draw(i, N, m, pop[t-1], odd_prob, L);
    }
  }
  
  // store COI distribution
  vector<int> coi(N);
  for (int i = 0; i < N; ++i) {
    coi[i] = pop[tlim-1][i].coi;
  }
  
  // store ancestry array
  vector<vector<vector<int>>> recomb(tlim);
  vector<vector<int>> parent_host1(tlim);
  vector<vector<int>> parent_host2(tlim);
  vector<vector<int>> parent_haplo1(tlim);
  vector<vector<int>> parent_haplo2(tlim);
  
  for (int t = 0; t < tlim; ++t) {
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < pop[t][i].coi; ++j) {
        recomb[t].push_back(pop[t][i].haplo_vec[j].parent_vec);
        parent_host1[t].push_back(pop[t][i].haplo_vec[j].pat_ind);
        parent_host2[t].push_back(pop[t][i].haplo_vec[j].mat_ind);
        parent_haplo1[t].push_back(pop[t][i].haplo_vec[j].pat_hap);
        parent_haplo2[t].push_back(pop[t][i].haplo_vec[j].mat_hap);
        
      }
    }
  }
  
  // end timer
  chrono_timer(t1);
  
  return Rcpp::List::create(Rcpp::Named("coi") = coi,
                            Rcpp::Named("recomb") = recomb,
                            Rcpp::Named("parent_host1") = parent_host1,
                            Rcpp::Named("parent_host2") = parent_host2,
                            Rcpp::Named("parent_haplo1") = parent_haplo1,
                            Rcpp::Named("parent_haplo2") = parent_haplo2);
}

