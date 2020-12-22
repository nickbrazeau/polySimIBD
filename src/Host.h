
#pragma once

#include "Haplotype.h"

#include <vector>

//------------------------------------------------
// class defining a human host
class Host {
  
public:
  
  // PUBLIC OBJECTS
  
  // basic properties
  int coi;
  
  // vector of haplotypes
  std::vector<Haplotype> haplo_vec;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Host() {};
  
  // member functions
  void init(double mean_coi, int L);
  void draw(int i, int Ntot,
            double m, const std::vector<Host> &prev_pop,
            const std::vector<double> &odd_prob, int L,
            const std::vector<int> &demesize, 
            const std::vector<double> &mig_mat_prob);
            void print_host();
            
};
