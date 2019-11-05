
set.seed(1)
# define parameters
pos <- seq(1,1e3,5e2)
N <- 10
m <- 0.5
rho <- 1e-3
mean_coi <- 3
tlim = 1e3
max_alleles = 2

swf2vcfR <- function(pos, N, m, rho, mean_coi, tlim, max_alleles){

  swf <- polySimIBD::sim_structured_WF(pos = pos,
                                       N = N,
                                       m = m,
                                       rho = rho,
                                       mean_coi = mean_coi,
                                       tlim = tlim)


}
