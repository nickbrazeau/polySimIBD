#..............................................................................
# Purpose of this script is to generate a series of
# simulations from our structured WF model to explore
# any bias of the model
#..............................................................................
library(rslurm)
library(tidyverse)
remotes::install_github("nickbrazeau/polySimIBD"); library(polySimIBD)


# define parameters
mean_coi <- c(10, 50, 1e2, 5e2, 1e3)
K <- seq(from = 0, to = 100, by = 10) # one deme
K[1] <- 1
m <- seq(from = 0, to = 1, by = 0.05)
m[1] <- 0.01
pos <- seq(0,1e3,1e2)
rho <- 1e-4

paramsdf <- expand.grid(mean_coi, K, m)
paramsdf <- lapply(1:5e2, function(x) return(paramsdf)) %>%
  dplyr::bind_rows() %>%
  magrittr::set_colnames(c("mean_coi", "K", "m"))

# add in details
paramsdf$pos <- list(pos)
paramsdf$rho <- rho

# now expand out again for sample size
iters <- nrow(paramsdf)
sample_size = c(2,3,5)

paramsdf <- lapply(1:length(sample_size), function(x) return(paramsdf)) %>%
  dplyr::bind_rows()
paramsdf$sample_size <- rep(sample_size, iters)


# sim wrapper

simwrapper <- function(pos, K, m, rho, mean_coi, sample_size){
  swf <- sim_structured_WF(pos = pos,
                           K = K,
                           m = m,
                           rho = rho,
                           mean_coi = mean_coi)

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
                            slurm_options = list(mem = 32000,
                                                 array = sprintf("0-%d%%%d",
                                                                 ntry,
                                                                 128),
                                                 'cpus-per-task' = 1,
                                                 error =  "%A_%a.err",
                                                 output = "%A_%a.out",
                                                 time = "24:00:00"))
