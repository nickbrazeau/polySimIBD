#--------------------------------------------------------------------
# Purpose of this script is to estimate the bias of the
# MLE inbreeding function from the Verity Nat Comm Submission
# using a discrete loci, discrete time structured wright fisher model
#--------------------------------------------------------------------
set.seed(1)


# Cpp wrapper of mipanalyzer inbreeding coeff
wrap_MIPanalyzer_inbreeding_mle_cpp <- function(WSAF.list,
                                                f = seq(0, 1, l = 11),
                                                ignore_het = F,
                                                report_progress = F){


  # extract WSAD
  wsaf <- WSAF.list$NRWSAFdf[,2:ncol(WSAF.list$NRWSAFdf)]
  wsaf <- as.matrix(wsaf)

  # get population allele frequencies based on haplotype biallelic matrix
  p <- rowMeans(WSAF.list$haplotypematrix.biall, na.rm = TRUE)

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

  # note because we will always have 2 here, we can just call upper tri (Bob stores upper tri not lower)
  ret <- ret$mle[ upper.tri(ret$mle, diag = F) ] # IBD value

  return(ret)
}




#--------------------------------------------------------------------
# Magic Number for Consideration
#--------------------------------------------------------------------
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
rho <- 1e-3
mean_coi <- 3
tlim <- 1e2

swfsim <- polySimIBD::sim_structured_WF(pos = pos,
                                        N = N,
                                        m = m,
                                        rho = rho,
                                        mean_coi = mean_coi,
                                        tlim = tlim)


smpl1 <- 1:swfsim$coi[1]
smpl2 <- (swfsim$coi[1] + 1):(cumsum(swfsim$coi[1:2])[2])

# extract ARG
ARG <- polySimIBD::get_ARG(swfsim, parasites = c(smpl1, smpl2))

# layer on mutations
hapmat <- polySimIBD::layer_mutations_on_ARG(ARG, mutationrate = 0.1)

# simulate biallelic reads
WSAF.list <- polySimIBD::sim_biallelic(COIs = swfsim$coi[1:2],
                                       haplotypematrix = hapmat,
                                       coverage = 100,
                                       alpha = 1,
                                       overdispersion = 0.1,
                                       epsilon = 0.05)

ret <- wrap_MIPanalyzer_inbreeding_mle_cpp(
  WSAF.list = WSAF.list,
  f = seq(0, 1, l = 11),
  ignore_het = F,
  report_progress = F)



