
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
// draws parents of every hapotype from within deme k
void Host::draw(int i, int N, double m, const vector<Host> &prev_pop,
                const vector<double> &odd_prob, int L) {
  
  // loop through every haplotype
  for (int i = 0; i < coi; ++i) {
    
    // draw parental hosts
    int parent_host1 = rbernoulli1(m) ? sample2(0, N-1) : i;
    int parent_haplo1 = sample2(0, prev_pop[parent_host1].coi-1);
    
    int parent_host2 = rbernoulli1(m) ? sample2(0, N-1) : i;
    int parent_haplo2 = sample2(0, prev_pop[parent_host2].coi-1);
    
    // fill in haplotype details
    haplo_vec[i].populate(parent_host1, parent_haplo1,
                          parent_host2, parent_haplo2,
                          odd_prob, L);
    
  }
  
  
}

//------------------------------------------------
// print
void Host::print_host() {
  print("COI:", coi);
}
