#..............................................................................
# Purpose of this script is to generate a series of
# simulations from our structured WF model to explore
# any bias of the model
#..............................................................................
library(rslurm)
library(tidyverse)
remotes::install_github("nickbrazeau/polySimIBD"); library(polySimIBD)


# define parameters
mean_coi <- c(1, 2, 5, 10, 25, 50, 1e2)
N <- c(seq(from = 0, to = 100, by = 10), seq(from = 100, to = 1e3, by = 1e2))
N <- unique(N)
N[1] <- 1
m <- seq(from = 0, to = 0.5, by = 0.1)
m[1] <- 0.01
pos <- seq(0,1e3,1e2)
rho <- 1e-4
tlim <- 5e2

paramsdf <- expand.grid(mean_coi, N, m)
paramsdf <- lapply(1:1e2, function(x) return(paramsdf)) %>%
  dplyr::bind_rows() %>%
  magrittr::set_colnames(c("mean_coi", "N", "m"))

# add in details
paramsdf$pos <- list(pos)
paramsdf$rho <- rho
paramsdf$tlim <- tlim

# now expand out again for sample size
iters <- nrow(paramsdf)
sample_size = c(2,2,2,3,5,10,15)

paramsdf <- lapply(1:length(sample_size), function(x) return(paramsdf)) %>%
  dplyr::bind_rows()
paramsdf$sample_size <- rep(sample_size, iters)

# remove 1 deme, mean coi of 1 (no pairwise to be made)


# sim wrapper
simwrapper <- function(pos, N, m, rho, mean_coi, sample_size, tlim){
  swf <- sim_structured_WF(pos = pos,
                           N = N,
                           m = m,
                           rho = rho,
                           mean_coi = mean_coi,
                           tlim = tlim)

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
                            jobname = 'edasim_Forward_sWF_runs',
                            nodes = ntry,
                            cpus_per_node = 1,
                            submit = T,
                            slurm_options = list(mem = 8000,
                                                 'cpus-per-task' = 1,
                                                 error =  "%A_%a.err",
                                                 output = "%A_%a.out",
                                                 time = "1-00:00:00"))
