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

pos <- seq(0,1e6,5e5) # assuming 1 million base-pairs and a SNP every 10,000 bp
N <- 100 # small Effective Pop size
m <- 0.1 # intermediate co-transmission, superinfxn
rho <- 1e-3
mean_coi <- 2
tlim <- 10
genome_length <- 23e6
mut_rate <- 2.45e-10 * 23e6/length(pos) * tlim # treating loci as blocks
hosts <- 1:2


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


# get true IBD
trueIBD <- get_truth_from_arg(swfsim = swfsim,
                              arg = ARG,
                              t_lim = 10,
                              hosts = hosts)



# extract relevant elements
trueIBD <- trueIBD %>% 
  dplyr::mutate(IBDprop = purrr::map(IBD, "IBDprop")) %>% 
  tidyr::unnest(cols = IBDprop) %>% 
  dplyr::select(c("smpl1", "smpl2", "IBDprop"))

# note trueIBD has both relationships
# ignores transitivity

#plot_coalescence_trees(ARG)

#..............................................................
# Add noise
#..............................................................
# layer on mutations
hapmat <- polySimIBD::layer_mutations_on_ARG(ARG, 
                                             mutationrate = mut_rate, t_lim = 10)

# simulate biallelic reads
WSAF.list <- polySimIBD::sim_biallelic(COIs = this_coi,
                                       haplotypematrix = hapmat,
                                       shape1 = 0.1,
                                       shape2 = 0.1,
                                       coverage = 100,
                                       alpha = 1,
                                       overdispersion = 0.1,
                                       epsilon = 0.05)


# run Bob's MLE 
ret <- wrap_MIPanalyzer_inbreeding_mle_cpp(
  WSAF.list = WSAF.list,
  f = seq(0, 1, l = 100),
  ignore_het = F,
  report_progress = F)


#..............................................................
# DATA WRANGLE truth and MLE
#..............................................................
colnames(ret$mle) <- rownames(ret$mle) <- 1:length(this_coi)
ret.long <- broom::tidy(as.dist(ret$mle)) %>%  
  magrittr::set_colnames(c("smpl1", "smpl2", "malecotf")) %>% 
  dplyr::left_join(., y = trueIBD, by = c("smpl1", "smpl2"))

