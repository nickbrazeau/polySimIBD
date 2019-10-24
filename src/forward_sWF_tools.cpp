#include "forward_sWF_tools.h"
#include "misc_v7.h"
#include <algorithm>
using namespace std;

//------------------------------------------------
// find the ARG from our forward wright fisher model
Rcpp::List get_ARG_cpp(Rcpp::List args) {

  // get inputs from Rcpp format to base C++ format
  int L = rcpp_to_int(args("L"));
  int G =  rcpp_to_int(args("G"));
  int coi =  rcpp_to_int(args("coi"));
  vector<int> p1 = rcpp_to_vector_int(args("p1"));
  vector<int> p2 = rcpp_to_vector_int(args("p2"));

  Rcpp::List anctemp = args("anc");
  vector<vector<vector<int>>> anc =rcpp_to_array_int(anctemp);
  vector<vector<int>> coal_times(int(p1.size()), vector<int>(L));


  // nodes for bvtrees (using this for unique https://en.cppreference.com/w/cpp/algorithm/unique)
  vector<int> nodes = p1;
  sort(nodes.begin(), nodes.end());
  auto uninodes = unique(nodes.begin(), nodes.end());
  nodes.erase(uninodes, nodes.end());

  // objects for storing bvtree results
  vector<vector<int>> mincoaltime = vector<vector<int>> (int(nodes.size()), vector<int>(L));
  vector<vector<int>> connectors = vector<vector<int>> (int(nodes.size()), vector<int>(L));
  vector<vector<int>> order = vector<vector<int>> (int(nodes.size()), vector<int>(L));

  //---------------------------------
  // get t and c
  //---------------------------------
  for (int l = 0; l < L; l++) {
    for (int pair = 0; pair < int(p1.size()); pair++){ // for each pair
      // find the coal time by looping back through generations
      for(int g = G; g > 0; g--) {
        int id1 = p1[pair];
        int id2 = p2[pair];
        while (coal_times[pair][l] == 0) {
          if (id1 == id2) {
            coal_times[pair][l] = G - g + 1;
          }
        }
      }
    } // end for loop for find coal_times for a given loci -- have long format for pairwise time to coalescence


    // objects we need quickly loop through nodes to find the min coal times
    // for each node, we know that the p1 rows from our expanded-grid that we are interested are iterated by (COI*COI - COI)/COI
    // we can then iterate through these intervals to get the appropriate p1.
    // Then we can make our connections always look right if we make p2 be greater than p1 which is a further
    // subset of our rows
    vector<int> pair_interval_end = vector<int>(coi - 1); // -1 here for last node being root
    for (int i = 0; i < int(pair_interval_end.size()); i++){
      pair_interval_end[i] = (coi-1)*(i+1) - 1; //-1 at end for 0 based indexing

    }

    vector<int> pair_interval_start = vector<int>(coi);
    for (int i = 0; i < int(pair_interval_start.size()); i++){
      if (i == 0){
        pair_interval_start[i] = 0;
      } else {
        pair_interval_start[i] = pair_interval_end[i] - (coi-2) + i; //-2 for removing self and then for 0 based indexing; +i to skip over the combinations that we have already considered
      }
    }

    // loop through nodes to extract bv_tree information
    for (int n = 0; n < int(nodes.size())-1; n++){
      // subset coal time matrix to find min coal time for a given node
      vector<int> coal_times_loci = vector<int>(int(p1.size())); // copying coal times to a vector
      for (int i = 0; i < int(p1.size()); i++){
        coal_times_loci[i] = coal_times[i][l];
      }
      //subset mintime
      vector<int> coal_times_loci_now = slice(coal_times_loci, pair_interval_start[n], pair_interval_end[n]);
      //subset mintime look ahead
      vector<int> coal_times_loci_ahead = slice(coal_times_loci, pair_interval_start[n+1], *coal_times_loci.end());

      int mintime_now = *min_element(coal_times_loci_now.begin(), coal_times_loci_now.end());
      int mintime_ahead = numeric_limits<int>::max();

      if (n != int(nodes.size())-2){ // if it is the second to last node, don't look ahead at root
        int mintime_ahead = *min_element(coal_times_loci_ahead.begin(), coal_times_loci_ahead.end());
      }

      // now find connector that the min coal time belongs to
      vector<int> lineage = {};
      for (int i = 0; i < int(coal_times_loci_now.size()); i++){
        if (coal_times_loci_now[i] == mintime_now){
          lineage.push_back( p2[pair_interval_start[n]+i] );
        }
      }

      // catch if multiple lineages coalesce at the same time
      int final_lineage;

      if (int(lineage.size()) > 1){
        if(mintime_ahead <= mintime_now){ // look ahead to make sure we don't block future branches
          int final_lineage = lineage[*lineage.end()]; // pick farthest right
        } else {
          int final_lineage = lineage[*lineage.begin()]; //  pick first (most left)
        }
      }
      // store result for coal time (t) and for the lineage connector (c)
      mincoaltime[n][l] = mintime_now;
      connectors[n][l] = final_lineage;

    } // end for loop for nodes

    // store root for each loci -- last node must always be root (if we are always going right)
    mincoaltime[*nodes.end()][l] = -1;
    connectors[*nodes.end()][l] = -1;

    //---------------------------------
    // get z
    //---------------------------------
    vector<int> z = vector<int>(int(nodes.size()), -1);
    // temp copy of t to z
    for (int i = 0; i < int(nodes.size()-1); i++){
      if(mincoaltime[i][l] != -1){
        z[i] = mincoaltime[i][l];
      }
    }
    sort(z.begin(), z.end());

    for (int i = 0; i < int(z.size()); i++){
      for (int j = 0; j < int(nodes.size()); j++){
        if(mincoaltime[i][l] == z[i]){
          order[j][l] = i;
        }
      }
    }

  } // end for loop for loci

  //---------------------------------
  // return as Rcpp list
  //---------------------------------
  // Rcpp::List ret = Rcpp::List::create(Rcpp::Named("t") = mincoaltime,
  //                                     Rcpp::Named("c") = connectors,
  //                                     Rcpp::Named("z") = order,
  //                                     Rcpp::Named("coal_times") = coal_times);
  // return ret;

  Rcpp::List ret = Rcpp::List::create(Rcpp::Named("t") = -9);
  return ret;
}

