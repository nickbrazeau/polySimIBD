#--------------------------------------------------------------------
# Purpose of this script is to RUN the simulations that will allow us to 
# estimate the bias of the
# MLE inbreeding function from the Verity Nat Comm Submission
# using a discrete loci, discrete time structured wright fisher model
#--------------------------------------------------------------------
library(tidyverse)
remotes::install_github("nickbrazeau/polySimIBD")
library(polySimIBD)
library(rslurm)
set.seed(1)
source("R_ignore/NatComms_VerityAB2019_Sims/utils.R")

#---------------------------------
# Magic Numbers for Consideration
#---------------------------------
# Genomic positions based on setting between chromosomes to a very high number.
# Effectively this breaks apart chromosomes for our forward simulation recombination runs
#
# Miles et al. 2016 (PMC5052046) & Taylor et al. 2019 (PMC6707449) gives us a recombination rate by 7.4e-7 M/bp
# Aimee gets this number by taking the inverse of Mile's estiamte of the CO recombination rate of 13.5 kb/cM


pos <- readRDS("R_ignore/NatComms_VerityAB2019_Sims/simparams/sim_POS.rds")
rho <- 7.4e-7
tlim <- 10

# NB, we will consider hosts (N) on a logarithmic-base 10 scale 
N <- round(10^seq(1, 3, l = 11))

# Mean COIs for lambda
coilamdas <- readRDS("R_ignore/NatComms_VerityAB2019_Sims/simparams/optim_lambda.RDS") 

# various levels of M
M <- c(0, 0.25, 0.5, 1)

# expand out combinations
paramsdf <- expand.grid(N, coilamdas, M) %>% 
  tibble::as_tibble() %>% 
  magrittr::set_colnames(c("N", "mean_coi", "m"))

paramsdf <- paramsdf %>% 
  dplyr::mutate(pos = list(pos), 
                rho = rho, 
                tlim = tlim, 
                # take 25% of hosts
                hosts = round(N * 0.25),
                hosts = purrr::map(hosts, function(x){return(1:x)})
  )


# replicates of this framework
#reps <- 25
#paramsdf <- lapply(1:reps, function(x) return(paramsdf)) %>%
#  dplyr::bind_rows()



#..............................................................
# Wrapper Function
#..............................................................

nat_comm_sims_wrapper <- function(pos, N, m, mean_coi, rho, tlim, hosts){
  # run forward
  swfsim <- polySimIBD::sim_swf(pos = pos,
                                N = N,
                                m = m,
                                rho = rho,
                                mean_coi = mean_coi,
                                tlim = tlim)
  
  # extract ARG
  ARG <- polySimIBD::get_arg(swfsim)
  
  
  #..............................................................
  # Downsample ARG
  #..............................................................
  
  this_coi <- swfsim$coi[hosts]
  # find haplotype index within the host
  w <- which(rep(1:length(swfsim$coi), times = swfsim$coi) %in% hosts)
  ARG <- mapply(function(x) polySimIBD::subset_bvtree(x, w), ARG)
  
  
  #..............................................................
  # Extract Haplotype Matrix
  #..............................................................
  hapmat <- polySimIBD::get_haplotype_matrix(ARG)
  
  
  #..............................................................
  # Simulate Reads
  #..............................................................
  WSAF.list <- polySimIBD::sim_biallelic(COIs = this_coi,
                                         haplotypematrix = hapmat,
                                         shape1 = 5.1,
                                         shape2 = 5.1,
                                         coverage = 100,
                                         alpha = 1,
                                         overdispersion = 0.1,
                                         epsilon = 0.05)  
  
  #..............................................................
  # Get True IBD 
  #..............................................................
  trueIBD <- get_truth_from_arg(swfsim = swfsim,
                                arg = ARG,
                                hosts = hosts)
  
  trueIBD.btwn <- trueIBD %>% 
    dplyr::mutate(btwnIBDprop = purrr::map(btwnIBD, "btwn_IBDprop")) %>% 
    dplyr::select(-c("btwnIBD")) %>% 
    tidyr::unnest(cols = btwnIBDprop)
  
  #..............................................................
  # Run Bob's MLE
  #..............................................................
  ret <- wrap_MIPanalyzer_inbreeding_mle_cpp(
    WSAF.list = WSAF.list,
    f = seq(0, 1, l = 100),
    ignore_het = F,
    report_progress = F)
  
  
  
  #..............................................................
  # DATA WRANGLE truth and MLE
  #..............................................................
  ret.long <- broom::tidy(ret$mle) %>%  
    magrittr::set_colnames(c("smpl1", "smpl2", "malecotf")) %>% 
    dplyr::left_join(., y = trueIBD.btwn, by = c("smpl1", "smpl2"))
  
  
  #..............................................................
  # RETURN
  #..............................................................
  ret <- list(sim_out = ret.long,
              ibdmle = ret, 
              WSAF.list = WSAF.list,
              hapmuts = hapmat,
              trueIBD = trueIBD,
              ARG = ARG,
              swfsim = swfsim
  )
  
  return(ret)
  
}


#..............................................................
# Run Sims
#..............................................................
# for slurm on LL
dir.create("R_ignore/NatComms_VerityAB2019_Sims/results/", recursive = T)
setwd("R_ignore/NatComms_VerityAB2019_Sims/results/")
ntry <- 1028 # max number of nodes
sjob <- rslurm::slurm_apply(f = nat_comm_sims_wrapper,
                            params = paramsdf,
                            jobname = 'verity_nat_comm_sims',
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

cat("*************************** \n Submitted Permutation Sims \n *************************** ")







