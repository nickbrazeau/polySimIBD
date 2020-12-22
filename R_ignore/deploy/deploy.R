library(tidyverse)
devtools::load_all()


swf <- polySimIBD::sim_swf(pos = seq(0,1e3,1e2),
                           N = 1e3,
                           m = 0.5,
                           rho = 1e-4,
                           mean_coi = 3,
                           tlim = 10)
ARG <- polySimIBD::get_arg(swf, 
                           host_index =  c(1,2))
