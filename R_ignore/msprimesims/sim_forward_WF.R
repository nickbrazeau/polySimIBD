#..............................................................................
# Purpose of this script is to generate a series of
# simulations from our structured WF model to compare to msprime
#..............................................................................
library(rslurm)
library(tidyverse)
remotes::install_github("nickbrazeau/polySimIBD"); library(polySimIBD)


# define parameters
mean_coi <- 1e3 # note, for msprime we used 5e2, but that was diploid. here we are haploid
N <- 1 # one ind/deme
m <- 1 # panmictic
pos <- seq(0,1e3,1e2)
rho <- 1e-4
tlim <- 1e9 # set to something high since we want these to finish on own

paramsdf <- tibble::tibble(
  mean_coi = mean_coi,
  N = N,
  m = m,
  pos = list(pos),
  rho = rep(rho, 3e3), # do 1e3 replicates again
  sample_size = rep(c(2,3,5), 1e3),
  tlim = tlim
)


simwrapper <- function(pos, N, m, rho, mean_coi, sample_size, tlim){
  swf <- sim_structured_WF(pos = pos,
                           N = N,
                           m = m,
                           rho = rho,
                           mean_coi = mean_coi,
                           tlim = tlim)

  # temp
  nodes <- 1:sample_size
  ARG <- get_ARG(swf, nodes = nodes)

  return(ARG)

}



outdir <- "/pine/scr/n/f/nfb/Projects/polySimIBD/"
dir.create(outdir, recursive = T)
setwd(outdir)
ntry <- 1028 # max number of nodes we can ask for on LL
sjob <- rslurm::slurm_apply(f = simwrapper,
                            params = paramsdf,
                            jobname = 'sim_Forward_sWF_runs',
                            nodes = ntry,
                            cpus_per_node = 1,
                            submit = T,
                            slurm_options = list(mem = 32000,
                                                 'cpus-per-task' = 1,
                                                 error =  "%A_%a.err",
                                                 output = "%A_%a.out",
                                                 time = "5-00:00:00"))
