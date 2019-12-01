
#------------------------------------------------
#' @title The Structured Wright Fisher Model for IBD
#'
#' @description Simulate a population forwards with recombination that approximates the Structured Wright Fisher Process
#' and tracks haplotype identity by descent where individuals represent demes, such that within
#' a deme individual-level COI is considered.
#'
#' @param pos dataframe <factor><numeric>; the genomic coordinates for chromosome and position of the sites
#' @param N numeric; The number of individuals to consider
#' @param m numeric; Probability of migration where m represents the probability of moving from \eqn{host_{origin}} to \eqn{host_{new}} by \eqn{m*(1-1/N)}
#' @param mean_coi numeric; The lambda of a right-shifted Poisson process, \eqn{1 + Pos(lambda)} representing the average COI of an individual deme
#' @param rho numeric; expected recombination rate
#' @param tlim numeric; the maximum number of generations to consider before exitting gracefully if all samples have not coalesced
#' 
#' @return Returns a list of length G (the number of generations it took for all lineages to coalesce). Each
#' list item is a matrix with n x L dimensions where n is the number of parasites considered and L is loci.
#' Parasites are then indexed into demes (individuals) based on the COI distribution.
#'
#' @export

sim_swf <- function(pos, N, m, rho, mean_coi, tlim){

  # assertions
  assert_vector(pos)
  assert_numeric(pos)
  assert_single_pos_int(N, zero_allowed = FALSE)
  assert_bounded(m, left = 0, right = 1)
  assert_bounded(rho, left = 0, right = 1, inclusive_left = FALSE, inclusive_right = FALSE)
  assert_single_pos(mean_coi, zero_allowed = FALSE)
  assert_single_pos_int(tlim, zero_allowed = FALSE)

  # warnings
  if (m == 0 & N > 1) {
    warning("You have set the migration rate to 0 but have more than one deme. As a result, all of your samples can never coalesce and this simulation will be limited by the tlim argument")
  }
  
  # precalculate probability of an odd number of recombination events between
  # any interval
  odd_prob <- exp(-rho*diff(pos))*sinh(rho*diff(pos))
  
  # define argument list
  args <- list(pos = pos,
               N = N,
               m = m,
               rho = rho,
               mean_coi = mean_coi,
               tlim = tlim,
               odd_prob = odd_prob)
  
  # run efficient C++ function
  output_raw <- sim_swf_cpp(args)

  return(output_raw)
}
