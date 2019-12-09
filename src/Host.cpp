
#include "Host.h"
#include "misc_v9.h"
#include "probability_v10.h"

using namespace std;

//------------------------------------------------
// initialise host
void Host::init(double mean_coi, int L) {
  
  // draw new COI
  coi = rztpois1(mean_coi);
  
  // create vector of haplotypes
  haplo_vec = vector<Haplotype>(coi);
  for (int i = 0; i < coi; ++i) {
    haplo_vec[i].init(L);
  }
  
}

//------------------------------------------------
// draws parents of every hapotype 
void Host::draw(int i, int N, double m, const vector<Host> &prev_pop,
                const vector<double> &odd_prob, int L) {
  
  // loop through every haplotype
  for (int j = 0; j < coi; ++j) {
    
    // draw parental host and parental haplos
    int parent_host = rbernoulli1(m) ? sample2(0, N-1) : i;
    int parent_haplo1 = sample2(0, prev_pop[parent_host].coi-1);
    int parent_haplo2 = sample2(0, prev_pop[parent_host].coi-1);
    
    // fill in haplotype details
    haplo_vec[j].populate(parent_host, parent_haplo1,
                          parent_host, parent_haplo2,
                          odd_prob, L);
    
  }
  
  
}

//------------------------------------------------
// print
void Host::print_host() {
  print("COI:", coi);
}
