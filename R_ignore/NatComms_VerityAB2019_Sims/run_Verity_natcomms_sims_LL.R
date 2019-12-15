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
m <- c(0, 0.25, 0.5, 1)

# expand out combinations
paramsdf <- expand.grid(N, coilamdas, m) %>% 
  tibble::as_tibble() %>% 
  magrittr::set_colnames(c("N", "mean_coi", "m"))

paramsdf <- paramsdf %>% 
  dplyr::mutate(pos = list(pos), 
                rho = rho, 
                tlim = tlim, 
                hosts = 1:2)


# replicates of this framework
reps <- 100
paramsdf <- lapply(1:reps, function(x) return(paramsdf)) %>%
  dplyr::bind_rows()



#..............................................................
# Wrapper Function
#..............................................................

nat_comm_sims_wrapper <- function(pos, N, m, mean_coi, rho, tlim, hosts){
  
  #..............................................................
  # Internal Function Start
  #..............................................................
  #------------------------
  # Bob MLE Calculator
  #------------------------
  
  # Cpp wrapper of mipanalyzer inbreeding coeff
  wrap_MIPanalyzer_inbreeding_mle_cpp <- function(WSAF.list,
                                                  f = seq(0, 1, l = 11),
                                                  ignore_het = F,
                                                  report_progress = F){
    
    
    # extract WSAF
    wsaf <- WSAF.list$NRWSAF
    
    # get population allele frequencies based on haplotype biallelic matrix
    # p <- rowMeans(wsaf, na.rm = TRUE)
    
    # get population allele frequencies based on our simulate beta dist
    p <- WSAF.list$rbetaPLAF
    
    
    #------------------------------------
    # liftovers needed for mipanalyzer
    #------------------------------------
    # progress bar
    pb <- txtProgressBar(min = 0, max = nrow(wsaf) - 1, initial = NA,
                         style = 3)
    args_progress <- list(pb = pb)
    
    if (ignore_het) {
      wsaf[wsaf != 0 & wsaf != 1] <- NA
    } else {
      wsaf <- round(wsaf)
    }
    wsaf[is.na(wsaf)] <- -1
    
    # call mipanalyzer internally
    args <- list(x = MIPanalyzer:::mat_to_rcpp(t(wsaf)), f = f, p = p, report_progress = report_progress)
    args_functions <- list(update_progress = MIPanalyzer:::update_progress)
    output_raw <- MIPanalyzer:::inbreeding_mle_cpp(args, args_functions, args_progress)
    ret_ml <- MIPanalyzer:::rcpp_to_mat(output_raw$ret_ml)
    ret_ml[row(ret_ml) >= col(ret_ml)] <- NA
    ret_all <- MIPanalyzer:::rcpp_to_array(output_raw$ret_all)
    ret <- list(mle = ret_ml, loglike = ret_all)
    
    # note Bob stores upper tri not lower,
    ret$mle <- as.dist(t(ret$mle))
    return(ret)
  }
  
  
  
  #------------------------
  # pull out overlap
  #------------------------
  get_truth_from_arg <- function(swfsim, arg, hosts = NULL){
    
    # need these details in order to know
    # which hosts you subsetted to to prune the
    # swfsim since the ARG doesn't store this information
    # choose hosts to subset to
    if(is.null(hosts)){
      hosts <- 1:length(swfsim$coi)
    }
    # find which elements in sim2 bvtrees correspond to haplotypes from these hosts
    this_coi <- swfsim$coi[hosts]
    
    
    # convert trees into matrix of alleles
    allele_mat <- polySimIBD::get_haplotype_matrix(arg)
    
    # split the haplotype matrix into individual (host) matrices 
    hosts.haplotypes <- NULL
    splitter <- rep(x = 1:length(hosts), times = this_coi)
    for (i in 1:length(unique(splitter))) {
      hosthap <- allele_mat[, c( splitter == i ), drop = F]
      hosts.haplotypes <- c(hosts.haplotypes, list(hosthap))
    }
    
    #..............................................................
    # Find true IBD between samples
    #..............................................................
    # expand grid for combinations
    paircompar.long <- expand.grid(list( 1:length(hosts), 1:length(hosts) ) ) %>% 
      magrittr::set_colnames(c("smpl1", "smpl2")) %>% 
      dplyr::filter(smpl1 != smpl2) %>% 
      tibble::as_tibble()
    
    paircompar.long$btwnIBD <- purrr::pmap(paircompar.long[,c("smpl1", "smpl2")], 
                                           function(smpl1, smpl2){
                                             # get mat for pairwise
                                             allele_mat_i <- hosts.haplotypes[[smpl1]]
                                             allele_mat_j <- hosts.haplotypes[[smpl2]]
                                             
                                             # find number of haplotypes that are IBD between hosts
                                             overlap <- mapply(function(x) {
                                               length(intersect(allele_mat_i[x,], allele_mat_j[x,]))
                                             }, 1:nrow(allele_mat_i))
                                             
                                             # just want if any overlap
                                             overlap[overlap >= 1] <- 1
                                             
                                             # return
                                             ret <- list(
                                               locioverlap = overlap,
                                               btwn_IBDprop = sum(overlap)/length(overlap)
                                             )
                                             return(ret)
                                             
                                           })
    
    return(paircompar.long)
  }
  
  #..............................................................
  # Internal Function End
  #..............................................................
  
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
  # ret <- list(sim_out = ret.long,
  #             ibdmle = ret, 
  #             WSAF.list = WSAF.list,
  #             hapmuts = hapmat,
  #             trueIBD = trueIBD,
  #             ARG = ARG,
  #             swfsim = swfsim
  # )
  # 
  # return(ret)
  
  return(ret.long)
  
}

paramsdf$simout <- purrr::pmap(paramsdf, nat_comm_sims_wrapper)

#..............................................................
# Run Sims
#..............................................................
# for slurm on LL
dir.create("R_ignore/NatComms_VerityAB2019_Sims/results/", recursive = T)
setwd("R_ignore/NatComms_VerityAB2019_Sims/results/")
ntry <- 1028 # max number of nodes
#ntry <- nrow(paramsdf)
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







