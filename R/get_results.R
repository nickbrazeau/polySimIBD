
#------------------------------------------------
#' internal class
#' @noRd
setClass("bvtree",
         slots=list(c="numeric", t="numeric", z="numeric"))


#' @title Extract Effective COI by Loci from SWF Simulation for a Single Host
#' @inheritParams get_swf
#' @details Only accepts a single host 
#' @description TODO 
#' @return vector of effective COI by loci
#' @export
get_effective_coi <- function(swf, host_index = NULL) {
  
  # checks
  assert_custom_class(swf, "swfsim")
  assert_single_pos_int(host_index)
  
  # get ARG from swf and host_index
  arg <- polySimIBD:::quiet(polySimIBD::get_arg(swf = swf, host_index = host_index))
  
  # get connections
  conn <- purrr::map(arg, "c")
  # get effective COI over loci 
  effCOI <- purrr::map_dbl(conn, function(x){sum(x == -1)})
  return(effCOI)
}



#' @inheritParams get_swf
#' @details Only accepts a single host 
#' @description TODO
#' @return double of within-host IBD
#' @export
get_within_IBD <- function(swf, host_index = NULL) {
  # checks
  assert_custom_class(swf, "swfsim")
  assert_single_pos_int(host_index)
  # get ARG from swf and host_index
  # need to call ARG again to make extraction straightforward (eg don't know what user did to ARG upstream)
  arg <- polySimIBD:::quiet(polySimIBD::get_arg(swf = swf, host_index = host_index))
  # get connections
  conn <- purrr::map(arg, "c")
  # get within host IBD per loci 
  numerator <- purrr::map_dbl(conn, function(x){sum(x != -1)})
  # -1 here for self comparison
  denom <-  (swf$coi[[host_index]]-1)
  # under SNP vs PSMC (Li/Durbin model) don't know begin and end, so treat as missing info
  wi <- diff(swf$pos)/sum(diff(swf$pos))
  # out
  wthnIBD <-sum( (numerator[1:(length(numerator) - 1)] / denom) * wi )
  return(wthnIBD)
}



#------------------------------------------------
#' @title Get Connection Intervals 
#' @description Index where unique connections are in the entire genome for proper weighting
#' @param uniqueconn unique bvtree connections from the ARG
#' @param allconn all bvtree connections from the ARG
#' @noRd
#' @noMd

get_conn_intervals <- function(uniqueconn, allconn){
  names(uniqueconn) <- 1:length(uniqueconn)
  intervals <- lapply(uniqueconn,
                      function(uni){
                        return(sapply(allconn, function(x){paste(uni, collapse = "") == paste(x, collapse = "")}))})
  mint <- rep(NA, length(allconn))
  for(i in 1:length(intervals)) {
    mint[intervals[[i]]] <- names(intervals)[i]
  }
  return(as.numeric(mint))
}



#' @title Between Host (pairwise) IBD 
#' @description Given an object \code{swf}, ***
#' @inheritParams get_swf
#' @description ***
#' @details  Only accepts a pair of hosts (i.e. pairwise). Ignores mutations as interrupting IBD segments. 
#' @return ***
#' @export
get_pairwise_ibd <- function(swf, host_index = NULL) {
  # check inputs and define defaults
  assert_custom_class(swf, "swfsim")
  assert_vector(host_index)
  assert_noduplicates(host_index)
  assert_pos_int(host_index, zero_allowed = FALSE)
  if(length(host_index) != 2) {
    stop("host_index must be of length 2 for pairwise comparison", call. = FALSE)
  }
  # get ARG from swf and host_index
  # need to call ARG again to make extraction straightforward (eg don't know what user did to ARG upstream)
  arg <- polySimIBD:::quiet(polySimIBD::get_arg(swf = swf, host_index = host_index))
  # subset to unique loci for speed 
  conn <- purrr::map(arg, "c")
  uniconn <- unique(conn)
  # store indices
  conn_indices <- get_conn_intervals(uniqueconn = uniconn, allconn = conn)
  # define arguments
  haplo_index <- mapply(function(x) 1:x, swf$coi[host_index], SIMPLIFY = FALSE)
  host_haplo_cnt <- mapply(length, haplo_index)
  argums <- list(conn = uniconn, 
                 host_haplo_cnt = host_haplo_cnt)
  
  # we define btwn_pairwise_ibd as: 
  # (n_{coal-btwn} + n_{coal-win-btwn}) / (n_{strains1} * n{strains2})
  # where w/in coal share a btwn coal 

  # pass to efficient C++ function for quick between tree look up 
  # to determine between host IBD
  output_raw <- calc_between_IBD_cpp(argums)
  
  #tidy raw
  # catch if no btwn, no w/in 
  if(sum(output_raw == 0)) { return(0)}
  
  # within host IBD 
 host1_wibd <- polySimIBD:::get_perloci_within_IBD(swf, host_index = host_index[1])
 host2_wibd <- polySimIBD:::get_perloci_within_IBD(swf, host_index = host_index[2])
 
 # numerator of btwn and w/in IBD
  numerator <- # TODO 
 
 # under SNP vs PSMC (Li/Durbin model) don't know begin and end, so treat as missing info
 wi <- diff(swf$pos)/sum(diff(swf$pos))
  # TODO 
  
  
  
  # out
  out <- # TODO 
  return(out)
}




#' Extract haplotypes from ARG
#' @param arg set of bvtrees
#' @return hapmat numeric matrix; a matrix of mutliallelic haplotypes for each parasite considered. Loci are in
#' rows and parasites (haplotypes) are in columns. 
#' @export
get_haplotype_matrix <- function(arg){
  
  # convert trees into matrix of alleles
  # each column is therefore a haplotype since we consider parasite by parasite
  hap_mat <- t(mapply(function(x) {
    c <- x@c
    ret <- c
    ret[ret == -1] <- 1:sum(ret == -1)
    while (any(c != -1)) {
      w <- which(c == -1)
      c[-w] <- c[c[-w]+1]
      ret[-w] <- ret[ret[-w]+1]
    }
    return(ret)
  }, ARG))
  return(hap_mat)
}
