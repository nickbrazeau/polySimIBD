
#' internal class
#' @noRd
setClass("bvtree",
         slots=list(c="numeric", t="numeric", z="numeric"))


#' @title Extract Effective COI by Loci from SWF Simulation for a Single Host
#' @inheritParams get_arg
#' @description From a single host in a SWF Simulation, extract the effective COI 
#'     for loci within the ARG. Effective COI is defined as the number of non-coalesced genomes 
#'     at the end of \code{tlim}. For example, if a host contains five strains (i.e. five parasites),
#'     then at loci 1 we can determine "mini" coalescent trees have occurred. If the effective COI is
#'     2, then three strains have formed one coalescent tree and two strains have formed a separate 
#'     coalescent tree. Notably, this framework is independent for each recombination event 
#'     (versus loci which are just genetic markers). This means, if loci 10 is located within a different 
#'     recombination block than loci 1 - and therefore is a separate entry in the ARG - we expect the 
#'     effective COI to be (potentially) different.
#' @details Function limited to a single host per "realization" 
#' @return vector of effective COI by loci
#' @export
get_effective_coi <- function(swf, host_index = NULL) {
  
  # checks
  goodegg::assert_class(swf, "swfsim")
  goodegg::assert_single_pos_int(host_index)
  
  # get ARG from swf and host_index
  arg <- polySimIBD::get_arg(swf = swf, host_index = host_index)
  
  # get connections
  conn <- purrr::map(arg, "c")
  # get effective COI over loci 
  effCOI <- purrr::map_dbl(conn, function(x){sum(x == -1)})
  return(effCOI)
}

#' @title Calculate Within-Host IBD 
#' @inheritParams get_arg
#' @description The within-host IBD is calculated as the number of strains that have 
#'     coalesced within the \code{tlim} at each loci divided by the original (i.e. not effective)
#'     COI. As an example, consider that there are three strains (i.e. parasites) within a host and
#'     that the parasite genome has ten equidistant loci with a single recombination breakpoint at loci 5 (i.e.).
#'     Within this framework, we consider at loci 1:5 if 2/3 strains have coalesced, the 
#'     within-host IBD for this section is 2/3. Next, for loci 6:10 if no strains have coalesced
#'     the within-host IBD is 0. Combining these results with-weighting for respective length/portion of the 
#'     genome (weights here equal and therefore negligible) the overall within-host IBD is: 
#'     \deqn{\frac{2 + 0}{Host_{COI} - 1}}, where one is subtracted from the Host-COI for self-comparison,
#'     wwhich gives (3-1) + (3-1) (for each loci). Note, because we consider self comparisons, the 
#'     denominator is always less than the true COI.
#' @details Function limited to a single host per "realization" 
#' @return double of within-host IBD
#' @export
get_within_ibd <- function(swf, host_index = NULL) {
  # checks
  goodegg::assert_class(swf, "swfsim")
  goodegg::assert_single_pos_int(host_index)
  goodegg::assert_gr(swf$coi[[host_index]], 1, message = "Cannot perform within IBD calculation when the host's COI is 1")
  # get ARG from swf and host_index
  # need to call ARG again to make extraction straightforward (eg don't know what user did to ARG upstream)
  arg <- polySimIBD::get_arg(swf = swf, host_index = host_index)
  # get connections
  conn <- purrr::map(arg, "c")
  # get within host IBD per loci 
  numerator <- sapply(conn, function(x){sum(x != -1)})
  denom <-  (swf$coi[[host_index]]-1) # -1 here for self comparison
  wi <- diff(swf$pos)/sum(diff(swf$pos))   # under SNP vs PSMC (Li/Durbin model) don't know begin and end, so treat as missing info
  # out
  wthnIBD <-sum( (numerator[1:(length(numerator) - 1)] / denom) * wi )
  return(wthnIBD)
}
