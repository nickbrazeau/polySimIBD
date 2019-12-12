#--------------------------------------------------------------------
# Purpose of this script is to estimate the bias of the
# MLE inbreeding function from the Verity Nat Comm Submission
# using a discrete loci, discrete time structured wright fisher model
#--------------------------------------------------------------------
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
N <- 10 # small Effective Pop size
m <- 0.5 # intermediate co-transmission, superinfxn
rho <- 1e-6
mean_coi <- 2
tlim <- 10
# genome_length <- 23e6
mut_rate <- 0 #2.45e-10 * 23e6/length(pos) 
hosts <- 1:10 

paramsdf <- tibble::tibble(
  N = c(10, 50, 100, 500, 1000),
  m = c(0, 0.25, 0.5, 0.75, 1),
  mean_coi = c(1, 2, 3, 4, 5)
)

paramsdf <- tibble::as_tibble(expand.grid(paramsdf))
paramsdf <- paramsdf %>% 
  dplyr::mutate(pos = list(pos), 
                rho = rho, 
                tlim = tlim, 
                mutationrate = mut_rate,
                hosts = list(hosts))

# replicates of this framework
#reps <- 25
#paramsdf <- lapply(1:reps, function(x) return(paramsdf)) %>% 
#   dplyr::bind_rows()



#..............................................................
# Wrapper Function
#..............................................................

nat_comm_sims_wrapper <- function(pos, N, m, mean_coi, rho, tlim, mutationrate, hosts){
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
  hapmat <- polySimIBD::layer_mutations_on_ARG(ARG, mutationrate = mut_rate, t_lim = tlim)
  
  
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


  trueIBD.btwn <- trueIBD$btwn_host_comparisons %>% 
    dplyr::mutate(btwnIBDprop = purrr::map(btwnIBD, "btwn_IBDprop"),
                  btwnmajIBDprop = purrr::map(majStrainIBD, "majstrain_IBDprop")
    ) %>% 
    tidyr::unnest(cols = c("btwnIBDprop", "btwnmajIBDprop")) %>% 
    dplyr::select(c("smpl1", "smpl2", "btwnIBDprop", "btwnmajIBDprop"))
  # note trueIBD has both relationships
  # ignores transitivity which is ok bc of left join later
 
  trueIBD.wthn <- trueIBD$wthn_host_comparisons %>% 
    dplyr::mutate(
      wthnIBD = purrr::map(wthnIBD, "within_IBDprop") ) %>% 
  tidyr::unnest(cols = wthnIBD) %>% 
  dplyr::select(c("hosts", "wthnIBD"))
  
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
  
  # add in within for sample 1
  trueIBD.wthn <- trueIBD.wthn %>% 
    dplyr::rename(smpl1 = hosts, 
                  wthnIBD.host1 = wthnIBD)
  
  ret.long <- dplyr::left_join(ret.long, trueIBD.wthn, by = "smpl1")
  
  # add in within for sample 2
  trueIBD.wthn <- trueIBD.wthn %>% 
    dplyr::rename(smpl2 = smpl1, 
                  wthnIBD.host2 = wthnIBD.host1)
  
  ret.long <- dplyr::left_join(ret.long, trueIBD.wthn, by = "smpl2")
  
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
paramsdf$nat_com_sims <- purrr::pmap(paramsdf, nat_comm_sims_wrapper)
#saveRDS(paramsdf, "~/Desktop/temp_natcomm/polysim_notmut.RDS")


#..............................................................
# Plot Results
#..............................................................
plot_theme <- theme(plot.title = element_text(family = "Helvetica", face = "bold", hjust = 0.5, size = 14),
                    axis.title = element_text(family = "Helvetica", face = "bold", hjust = 0.5, size = 12),
                    axis.text.y = element_text(family = "Helvetica", hjust = 0.5, size = 11),
                    axis.text.x = element_text(family = "Helvetica", hjust = 0.5, size = 11, angle = 90),
                    legend.position = "right",
                    legend.title = element_text(family = "Helvetica", face = "bold", vjust = 0.85, size = 12),
                    legend.text = element_text(family = "Helvetica", hjust = 0.5, vjust = 0.5, size = 10),
                    axis.line = element_line(color = "#000000", size = 1))

plotdf <- paramsdf %>% 
  dplyr::mutate(simout = purrr::map(paramsdf$nat_com_sims, "sim_out")) %>% 
  tidyr::unnest(cols = simout)

plotObj <- plotdf %>% 
  tidyr::gather(., key = "IBD", value = "IBDest", 12:13) %>% 
  dplyr::group_by(mean_coi, m, N, IBD) %>% 
  dplyr::summarise(
    n = n(),
    meanIBD = mean(IBDest > 0.9)
  ) %>% 
  dplyr::mutate(logN = log10(N)) %>% 
  ggplot() + 
  # geom_pointrange(aes(x = logN, y = meanIBD, ymin = IBDLL, ymax = IBDUL, 
  #                     group = factor(IBD), color = factor(IBD)), alpha = 0.5) +
  geom_point(aes(x = logN, y = meanIBD,
                 group = factor(IBD), color = factor(IBD)), alpha = 0.5) +                     
  scale_color_manual("IBD Measure", values = c("#542788", "#e08214")) +
  facet_grid(mean_coi ~ m) + 
  ylab("Between-Sample IBD > 90") + xlab("Effective Population (log10-transformed)") +
  plot_theme +
  ggtitle("Using Not Mut, WSAF as PLAF")

jpeg("~/Desktop/temp_polysimibd_meanIBD_nomut_rho.jpg", width = 11, height = 8, units = "in", res = 300)
plot(plotObj)
graphics.off()





