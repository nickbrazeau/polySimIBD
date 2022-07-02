
#include "main.h"
#include "probability_v10.h"
#include "Haplotype.h"
#include "Host.h"

#include <chrono>

using namespace std;

//------------------------------------------------
// converts input from Rcpp::List format to vector<vector<vector<vector<int>>>> format.
vector<vector<vector<vector<int>>>> rcpp_to_4d_int(Rcpp::List x) {
  int n1 = int(x.size());
  vector<vector<vector<vector<int>>>> ret(n1);
  for (int i = 0; i < n1; ++i) {
    Rcpp::List x_i = x[i];
    int n2 = int(x_i.size());
    ret[i] = vector<vector<vector<int>>>(n2);
    for (int j = 0; j < n2; ++j) {
      Rcpp::List x_ij = x_i[j];
      int n3 = int(x_ij.size());
      ret[i][j] = vector<vector<int>>(n3);
      for (int k = 0; k < n3; ++k) {
        ret[i][j][k] = Rcpp::as<vector<int>>(x_ij[k]);
      }
    }
  }
  return ret;
}

//------------------------------------------------
// draw from forwards-in-time model
Rcpp::List sim_swf_cpp(Rcpp::List args) {

  // start timer
  chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();

  // extract arguments
  vector<double> m = rcpp_to_vector_double(args["m"]);
  vector<double> mean_coi = rcpp_to_vector_double(args["mean_coi"]);
  int tlim = rcpp_to_int(args["tlim"]);
  vector<double> odd_prob = rcpp_to_vector_double(args["odd_prob"]);
  int L = int(odd_prob.size()) + 1;
  // extract out demes arguements
  vector<int> demesize = rcpp_to_vector_int(args["N"]);
  int demecnt = demesize.size();
  int N = 0;
  for (int i = 0; i < demecnt; ++i) {
    N += demesize[i];
  }

  // migration arguments
  vector<vector<double>> mig_mat_prob = rcpp_to_matrix_double(args["mig_mat_prob"]);

  // nested vectors, first over time, then individuals
  vector<vector<Host>> pop(tlim, vector<Host>(N));

  // draw COI for every individual in every time step as part of init
  for (int t = 0; t < tlim; ++t) {
    // need to use iters to properly count through N
    int citer = 0;
    for (int d = 0; d < demecnt; ++d) {
      for (int i = 0; i < demesize[d]; ++i) {
        pop[t][citer].init(mean_coi[d], L);
        citer++;
      }
    }
  }

  // step through time in discrete generations, starting at the second timestep
  for (int t = 1; t < tlim; ++t) {
    int diter = 0;
    for (int d = 0; d < demecnt; ++d) {
      for (int i = 0; i < demesize[d]; ++i) {
        pop[t][diter].draw(diter, N, m[d], pop[t-1], odd_prob,
                           L, demesize, mig_mat_prob[d]);
        diter++;
      }
    }
  }

  // create objects for storing results
  vector<int> coi(N);
  vector<vector<vector<vector<int>>>> recomb(tlim, vector<vector<vector<int>>>(N));
  vector<vector<vector<int>>> parent_host1(tlim, vector<vector<int>>(N));
  vector<vector<vector<int>>> parent_host2(tlim, vector<vector<int>>(N));
  vector<vector<vector<int>>> parent_haplo1(tlim, vector<vector<int>>(N));
  vector<vector<vector<int>>> parent_haplo2(tlim, vector<vector<int>>(N));

  // store COI distribution at final time point only
  for (int i = 0; i < N; ++i) {
    coi[i] = pop[tlim-1][i].coi;
  }

  // loop through time and individuals
  for (int t = 0; t < tlim; ++t) {
    for (int i = 0; i < N; ++i) {

      // preallocate for this coi
      int this_coi = pop[t][i].coi;
      recomb[t][i] = vector<vector<int>>(this_coi);
      parent_host1[t][i] = vector<int>(this_coi);
      parent_host2[t][i] = vector<int>(this_coi);
      parent_haplo1[t][i] = vector<int>(this_coi);
      parent_haplo2[t][i] = vector<int>(this_coi);

      // loop through all haplotypes in this individual
      for (int j = 0; j < this_coi; ++j) {
        recomb[t][i][j] = pop[t][i].haplo_vec[j].parent_vec;
        parent_host1[t][i][j] = pop[t][i].haplo_vec[j].pat_ind;
        parent_host2[t][i][j] = pop[t][i].haplo_vec[j].mat_ind;
        parent_haplo1[t][i][j] = pop[t][i].haplo_vec[j].pat_hap;
        parent_haplo2[t][i][j] = pop[t][i].haplo_vec[j].mat_hap;
      }
    }
  }

  // end timer
  chrono_timer(t1);

  // return as list
  return Rcpp::List::create(Rcpp::Named("coi") = coi,
                            Rcpp::Named("recomb") = recomb,
                            Rcpp::Named("parent_host1") = parent_host1,
                            Rcpp::Named("parent_host2") = parent_host2,
                            Rcpp::Named("parent_haplo1") = parent_haplo1,
                            Rcpp::Named("parent_haplo2") = parent_haplo2);
}

//------------------------------------------------
// walk back through ancestry, find coalescent events
Rcpp::List get_arg_cpp(Rcpp::List args) {

  // start timer
  chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();

  // extract arguments
  vector<int> coi_final = rcpp_to_vector_int(args["coi"]);
  vector<vector<vector<vector<int>>>> recomb = rcpp_to_4d_int(args["recomb"]);
  vector<vector<vector<int>>> parent_host1 = rcpp_to_array_int(args["parent_host1"]);
  vector<vector<vector<int>>> parent_host2 = rcpp_to_array_int(args["parent_host2"]);
  vector<vector<vector<int>>> parent_haplo1 = rcpp_to_array_int(args["parent_haplo1"]);
  vector<vector<vector<int>>> parent_haplo2 = rcpp_to_array_int(args["parent_haplo2"]);
  vector<int> host_index = rcpp_to_vector_int(args["host_index"]);
  vector<int> haplo_index = rcpp_to_vector_int(args["haplo_index"]);
  int L = recomb[0][0][0].size();
  int tlim = recomb.size();
  int n = host_index.size(); // number of haplotypes within hosts (cumulative)

  // objects for storing results
  vector<vector<int>> coalesce_target_store(L);
  vector<vector<int>> coalesce_time_store(L);

  // loop through loci
  for (int l = 0; l < L; ++l) {

    // initialise ancestry tracking
    vector<int> sample_host = host_index;
    vector<int> sample_haplo = haplo_index;

    // initialise objects for storing coalescent events
    vector<bool> coalesced(n, false);
    vector<int> coalesce_target(n, -1);
    vector<int> coalesce_time(n, -1);

    // loop through time from present going backwards
    for (int t = (tlim-1); t >= 1 ; --t) {

      // loop through samples
      for (int i = 0; i < n; ++i) {

        // skip already coalesced
        if (coalesced[i]) {
          continue;
        }

        // define for convenience
        int this_host = sample_host[i];
        int this_haplo = sample_haplo[i];

        // update ancestry tracking
        if (recomb[t][this_host][this_haplo][l]) {
          sample_host[i] = parent_host1[t][this_host][this_haplo];
          sample_haplo[i] = parent_haplo1[t][this_host][this_haplo];
        } else {
          sample_host[i] = parent_host2[t][this_host][this_haplo];
          sample_haplo[i] = parent_haplo2[t][this_host][this_haplo];
        }

      }

      // check for coalescnce
      for (int i = 1; i < n; ++i) {
        if (coalesced[i]) {
          continue;
        }
        for (int j = 0; j < i; ++j) {
          if (coalesced[j]) {
            continue;
          }
          if (sample_host[j] == sample_host[i]) {
            if (sample_haplo[j] == sample_haplo[i]) {

              // coalescence event
              coalesced[i] = true;
              coalesce_target[i] = j;
              coalesce_time[i] = tlim - t;

            }
          }
        }
      }  // end i loop

    }  // end t loop

    // store the results of this locus
    coalesce_target_store[l] = coalesce_target;
    coalesce_time_store[l] = coalesce_time;

  }  // end l loop

  // end timer
  chrono_timer(t1);

  // return as list
  return Rcpp::List::create(Rcpp::Named("coalesce_target") = coalesce_target_store,
                            Rcpp::Named("coalesce_time") = coalesce_time_store);
}

//------------------------------------------------
// subset bvtree to specified haplotypes only
Rcpp::List subset_bvtree_cpp(vector<int> c,
                             vector<int> t,
                             vector<int> m) {

  // get basic properties
  int n = c.size();

  // fix each haplotype in turn
  for (int i = 0; i < n; ++i) {
    if (m[i] == 0 || c[i] == -1) {
      continue;
    }
    int x = c[i];
    while (m[x] == 0) {
      c[i] = c[x];
      t[i] = t[x];
      if (c[x] == -1) {
        break;
      }
      x = c[x];
    }
  }

  // return list
  return Rcpp::List::create(Rcpp::Named("c") = c,
                            Rcpp::Named("t") = t);
}


//------------------------------------------------
// efficient loop through ARG backwards in time to calculate IBD
Rcpp::List calc_between_IBD_cpp(Rcpp::List args) {
  // extract arguments
  vector<vector<int>> conn = rcpp_to_matrix_int(args["conn"]);
  vector<int> host_haplo_cnt = rcpp_to_vector_int(args["host_haplo_cnt"]);
  int haptot = host_haplo_cnt[0] + host_haplo_cnt[1];
  int L = conn.size(); // loci here are unique from upstream R filtering

  // initialise IBD object for results
    vector<int> ibd_numerator(L, 0);

  // loop through loci
  for (int l = 0; l < L; ++l) {
    // always only two hosts for between IBD comparison
    // looking left to right
    for (int a = 0; a < host_haplo_cnt[0]; ++a) {
      for (int b = host_haplo_cnt[0]; b < haptot; ++b) {
        if (a == conn[l][b]) {
          ibd_numerator[l] = 1;
        }
      }
    } // end left to right loop

    // looking right to left
    for (int b = host_haplo_cnt[0]; b < haptot; ++b) {
      for (int a = 0; a < host_haplo_cnt[0]; ++a) {
        if (b == conn[l][a]) {
          ibd_numerator[l] = 1;
        }
      }
    } // end right to left loop

  } // end loci loop

  // return as list
  return Rcpp::List::create(Rcpp::Named("ibd_numerator") = ibd_numerator);
}
