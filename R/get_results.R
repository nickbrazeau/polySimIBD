
#------------------------------------------------
#' internal class
#' @noRd
setClass("bvtree",
         slots=list(c="numeric", t="numeric", z="numeric"))


#' @title Extract Effective COI by Loci from SWF Simulation for a Single Host
#' @inheritParams get_arg
#' @details Only accepts a single host 
#' @return vector of effective COI by loci
#' @export
get_realized_coi <- function(swf, host_index = NULL) {
  
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


#' @title Get Within-Host IBD from SWF Simulation
#' @inheritParams get_arg
#' @details Only accepts a single host 
#' @return double of within-host IBD
#' @export
get_within_host_IBD <- function(swf, host_index = NULL) {
  
  # checks
  assert_custom_class(swf, "swfsim")
  assert_single_pos_int(host_index)
  
  # get ARG from swf and host_index
  arg <- polySimIBD:::quiet(polySimIBD::get_arg(swf = swf, host_index = host_index))
  
  # get connections
  conn <- purrr::map(arg, "c")
  # get effective IBD over loci 
  numerator <- purrr::map_dbl(conn, function(x){sum(x != -1)})
  # -1 here for the SELF comparison
  denom <-  (swf$coi[[host_index]]-1) * length(conn)
  # out
  wthnIBD <- sum(numerator)/denom
  return(wthnIBD)
}


#' Extract haplotypes from ARG
#' @param ARG set of bvtrees
#' @return hapmat numeric matrix; a matrix of mutliallelic haplotypes for each parasite considered. Loci are in
#' rows and parasites (haplotypes) are in columns. 
#' @export
get_haplotype_matrix <- function(ARG){
  
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


#------------------------------------------------
#' @title Effective IBD by Loci from SWF for a Pair of Hosts
#' @description Given an object \code{SWF}, ***
#' @inheritParams get_arg
#' @description Assumes that the minimum realized COI between the pairs of host determines
#' the denominator for the between realized IBD
#' @details  Only accepts a pair of hosts. Ignores mutations as interrupting IBD segments. 
#' @return ***
#' @export


get_realized_pairwise_ibd <- function(swf, host_index = NULL) {
  
  # check inputs and define defaults
  assert_custom_class(swf, "swfsim")
  if (is.null(host_index)) {
    host_index <- 1:length(swf$coi)
  }
  assert_vector(host_index)
  assert_pos_int(host_index, zero_allowed = FALSE)
  
  # define arguments
  arg <- polySimIBD:::quiet(polySimIBD::get_arg(swf = swf, host_index = host_index))
  conn <- purrr::map(arg, "c")
  haplo_index <- mapply(function(x) 1:x, swf$coi[host_index], SIMPLIFY = FALSE)
  host_haplo_cnt <- mapply(length, haplo_index)
  argums <- list(conn = conn, 
                 host_haplo_cnt = host_haplo_cnt,
                 host_index = rep(host_index, mapply(length, haplo_index)) - 1,
                 haplo_index = unlist(haplo_index) - 1)
  
  # pass to efficient C++ function
  output_raw <- calc_between_IBD_cpp(argums)
  
  # tidy cpp output
  # TODO w/ numerator versus denominator
  
  # out
  return(output_raw)
}






