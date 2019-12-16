
#include "Haplotype.h"
#include "misc_v9.h"
#include "probability_v10.h"

using namespace std;


//------------------------------------------------
// constructor
Haplotype::Haplotype(int L) {
  init(L);
}

//------------------------------------------------
// initialise haplotype
void Haplotype::init(int L) {
  
  // populate parental fields
  mat_ind = -1;
  mat_hap = -1;
  
  pat_ind = -1;
  pat_hap = -1;
  
  // initialse parental vector
  parent_vec = vector<int>(L);
  
}

//------------------------------------------------
// populate haplotype fields
void Haplotype::populate(int parent1, int haplo1,
                         int parent2, int haplo2,
                         const vector<double> &odd_prob, int L) {
  
  // populate parental fields
  mat_ind = parent1;
  mat_hap = haplo1;
  
  pat_ind = parent2;
  pat_hap = haplo2;
  
  // draw starting parent
  int switcher = rbernoulli1(0.5);
  parent_vec[0] = switcher;
  
  // draw total number of recombination breakpoints
  for (int i = 1; i < L; ++i) {
    if (rbernoulli1(odd_prob[i-1])) {
      switcher = 1 - switcher;
    }
    parent_vec[i] = switcher;
  }
  
}

//------------------------------------------------
// print haplotype
void Haplotype::print_haplo() {
  
  print_stars();
  print("maternal:", mat_ind, mat_hap);
  print("paternal:", pat_ind, pat_hap);
  print_vector(parent_vec);
  
}
