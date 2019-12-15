#..............................................................
# Utilities for Nat Comm Simulations
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
