
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
// draws parents of every haplotype
void Host::draw(int i, int Ntot,
                double m, const vector<Host> &prev_pop,
                const vector<double> &odd_prob, int L,
                const vector<int> &demesize,
                const vector<double> &mig_mat_prob) {

  // loop through every haplotype
  for (int j = 0; j < coi; ++j) {

    // draw deme to sample parental host and haplotypes
    int this_deme = sample1(mig_mat_prob,1.0);

    // index demes to properly account for which host in N we are considering
    int deme_start = 0;
    if (this_deme != 0) {
      for (int i = 0; i <= (this_deme-1); ++i) {
        deme_start += demesize[i];
      }
    }
    int deme_end = 0;
    for (int i = 0; i <= this_deme; ++i) {
      deme_end += demesize[i];
    }

    // draw parental host and parental haplos
    int parent_host = rbernoulli1(m) ? sample2(deme_start, deme_end-1) : i;
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
