#..............................................................................
# Purpose of this script is to generate a series of
# simulations from our structured WF model to compare to msprime
#..............................................................................
library(rslurm)
library(tidyverse)
remotes::install_github("nickbrazeau/polySimIBD"); library(polySimIBD)


# define parameters
mean_coi <- 5e2 # note, we want mean_coi to be large, as for msprime, we simulated 1e3 individuals
K <- 1 # one deme
m <- 1 # panmictic
pos <- seq(0,1e3,1e2)
rho <- 1e-4


paramsdf <- tibble::tibble(
  mean_coi = mean_coi,
  K = K,
  m = m,
  pos = list(pos),
  rho = rep(rho, 5e4) # did 1e3 reps for msprime but doing more here for stochasticity we are entering in for geometric process
)




outdir <- "/pine/scr/n/f/nfb/Projects/polySimIBD/msprimesims/"
dir.create(outdir, recursive = T)
setwd("outdir")
ntry <- 1028 # max number of nodes we can ask for on LL
sjob <- rslurm::slurm_apply(f = sim_structured_WF,
                            params = paramsdf,
                            jobname = 'sim_Forward_sWF',
                            nodes = ntry,
                            cpus_per_node = 1,
                            submit = T,
                            slurm_options = list(mem = 32000,
                                                 array = sprintf("0-%d%%%d",
                                                                 ntry,
                                                                 128),
                                                 'cpus-per-task' = 1,
                                                 error =  "%A_%a.err",
                                                 output = "%A_%a.out",
                                                 time = "1:00:00"))
