library(simcross)
p1 <- create_parent(L=100, allele=1)
p2 <- create_parent(L=100, allele=2)

# https://github.com/kbroman/simcross/blob/master/src/sim_meiosis.cpp
# as long as we turn off interference this basically just becomes
# a poisson model from sim_meiosis
cross(p1, p2, m=0, p=0,
      xchr = F, obligate_chiasma = F)
