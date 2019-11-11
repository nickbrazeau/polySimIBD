
set.seed(1)
# define parameters
pos <- seq(1,1e3,1e2)
N <- 2
m <- 0.5
rho <- 1e-3
mean_coi <- 3
tlim = 1e3
mutationrate = 1e-2
coverage = 100
alpha = 1
overdispersion = 0
epsilon = 0
cutoff = 0.5

swf2bivcfR <- function(pos, N, m, rho, mean_coi, tlim, # swf params
                     mutationrate, # mutation rate for recombination block
                     coverage, alpha, overdispersion, epsilon, # params for simulating biallelic reads
                     cutoff # cutoff for vcfR
                     ){
  require(vcfR)
  #.......................
  # SWF & ARG
  #.......................
  swf <- polySimIBD::sim_structured_WF(pos = pos,
                                       N = N,
                                       m = m,
                                       rho = rho,
                                       mean_coi = mean_coi,
                                       tlim = tlim)

  # get COI
  coi <- swf$coi

  # get ARG
  ARG <- polySimIBD::get_ARG(swf) # get all nodes

  # layer on mutations
  hapmat <- polySimIBD::layer_mutations_on_ARG(ARG, mutationrate = mutationrate)

  # split up pairs
  p1 <- hapmat[,1:swf$coi[1]]
  p2 <- hapmat[,swf$coi[1]:ncol(hapmat)]

  #.......................
  # sim biallelic
  #.......................
  # sim data to biallelic
  p1 <- polySimIBD::sim_biallelic(haplotypematrix = p1,
                                  coverage = coverage,
                                  alpha = alpha,
                                  overdispersion = overdispersion,
                                  epsilon = epsilon)

  p2 <- polySimIBD::sim_biallelic(haplotypematrix = p2,
                                  coverage = coverage,
                                  alpha = alpha,
                                  overdispersion = overdispersion,
                                  epsilon = epsilon)

  #.......................
  # vcfR object
  #.......................
  # get pop AF from haplotype mats
  hapmat.bi <- purrr::map(list(p1, p2), "phased")
  hapmat.bi <- do.call("cbind", hapmat.bi)
  PLAF <- apply(hapmat.bi, 1, mean)

  # getwsraf
  gt.ad <- purrr::map(list(p1, p2), "counts")
  gt.ad <- do.call("cbind", gt.ad)
  gt.dp <- purrr::map(list(p1, p2), "coverage")
  gt.dp <- do.call("cbind", gt.dp)

  wsaf <- gt.ad/gt.dp

  # liftover to vcfR
  # round accordin to cutoff
  GT <- ifelse(wsaf > 1-cutoff, "0/0",
               ifelse(wsaf < 0+cutoff, "1/1",
                      ifelse(!is.na(wsaf), "0/1", NA)))

  gt <- matrix(NA, nrow = nrow(GT), ncol = ncol(GT))
  for(i in 1:ncol(wsaf)){
    for(j in 1:nrow(wsaf)){
      gt[j,i] <- paste0(GT[j,i], ":", gt.ad[j,i], ",", gt.dp[j,i] - gt.ad[j,i], ":", gt.dp[j,i])
    }
  }

  # append format column and sample names
  gt <- cbind(FORMAT = "GT:AD:DP", gt)
  colnames(gt)[2:ncol(gt)] <- paste0("smpl", 1:(ncol(gt)-1))

  # getFix
  fix <- matrix(NA, nrow = nrow(gt), ncol = 8 )
  colnames(fix) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
  fix[,1] <- "chrom1"
  fix[,2] <- swf$pos

  # get meta
  meta <- c("##fileformat=VCFv4.3")

  # write out new vcfRobj
  vcfRobj <- new("vcfR", meta = meta, fix = fix, gt = gt)

  #.......................
  # MIPanalyzer
  #.......................
  mipobj <- MIPanalyzer::vcf2mipanalyzer_biallelic(vcfR = vcfRobj)
  Fibd <- MIPanalyzer::inbreeding_mle(mipobj)




}
