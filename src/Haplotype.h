
#pragma once

#include <vector>

//------------------------------------------------
// class defining a haplotype
class Haplotype {
  
public:
  
  // PUBLIC OBJECTS
  
  // 'maternal' properties
  int mat_deme;
  int mat_ind;
  int mat_hap;
  
  // 'paternal' properties
  int pat_deme;
  int pat_ind;
  int pat_hap;
  
  // vector along the genome, true for maternal, false for paternal
  std::vector<int> parent_vec;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Haplotype() {};
  Haplotype(int L);
  
  // member functions
  void init(int L);
  void populate(int parent1, int haplo1,
                int parent2, int haplo2,
                const std::vector<double> &odd_prob, int L);
  void print_haplo();
  
};
