library(tidyverse)
library(polySimIBD)
set.seed(1)
source("R_ignore/NatComms_VerityAB2019_Sims/utils.R")
#---------------------------------
# Magic Number for Consideration
#---------------------------------
# From Hamilton et al. 2017 (PMC5389722), we can esimate a mutation rate of:
# 2.45 × 10−10 base-pair substitions/generation/basepair --- from Manuscript: The pooled mutation rate from all six P. falciparum isolates was 2.45 × 10−10 (95%CI, 1.70 × 10−10–3.20 × 10−10) BPS/ELC/bp
# therefore we need to scale 2.45e-10 to be sub/gen/recombo block
#
# Miles et al. 2016 (PMC5052046) & Taylor et al. 2019 (PMC6707449) gives us a recombination rate by 7.4e-7 M/bp
# Aimee gets this number by taking the inverse of Mile's estiamte of the CO recombination rate of 13.5 kb/cM
#

pos <- seq(0,1e6,1e4) # assuming 1 million base-pairs and a SNP every 10,000 bp
N <- 50 # small Effective Pop size
m <- 0.5 # intermediate co-transmission, superinfxn
rho <- 1e-6
mean_coi <- 2
tlim <- 10
genome_length <- 23e6
mut_rate <- 0 #2.45e-10 * 23e6/length(pos) * tlim # treating loci as blocks
hosts <- 1:5


#..............................................................
# Run sims
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
# Add noise
#..............................................................
# layer on mutations
hapmat <- polySimIBD::layer_mutations_on_ARG(ARG, 
                                             mutationrate = mut_rate, t_lim = 10)


#..............................................................
# Simulate Reads
#..............................................................
WSAF.list <- polySimIBD::sim_biallelic(COIs = this_coi,
                                       haplotypematrix = hapmat,
                                       shape1 = 0.1,
                                       shape2 = 0.1,
                                       coverage = 100,
                                       alpha = 1,
                                       overdispersion = 0.1,
                                       epsilon = 0.05)

#..............................................................
# Get True IBD 
#..............................................................
trueIBD <- get_truth_from_arg(swfsim = swfsim,
                              arg = ARG,
                              WSAFlist = WSAF.list,
                              t_lim = 10,
                              hosts = hosts)


trueIBD <- trueIBD %>% 
  dplyr::mutate(btwnIBD = purrr::map(btwn_host_comparisons, "btwnIBD"),
                trueIBDprop = purrr::map(btwnIBD, "true_IBDprop"),
                majStrainIBD = purrr:map(btwn_host_comparisons, "majStrainIBD"),
                majIBDprop = purrr::map(majStrainIBD, "majstrain_IBDprop"),
                wthnIBD = purrr::map(wthn_host_comparisons, "wthnIBD"),
                wthnIBDprop = purrr::map(wthn_host_comparisons, "within_IBDprop"),
                ) %>% 
  tidyr::unnest(cols = c("trueIBDprop", "majIBDprop", "within_IBDprop")) %>% 
  dplyr::select(c("smpl1", "smpl2", "trueIBDprop", "majIBDprop", "within_IBDprop"))

# note trueIBD has both relationships
# ignores transitivity which is ok bc of left join later

#plot_coalescence_trees(ARG)

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
  dplyr::left_join(., y = trueIBD, by = c("smpl1", "smpl2"))



