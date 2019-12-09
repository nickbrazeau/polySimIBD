#------------------------------------------------
#' @title The Structured Wright Fisher Model for IBD
#'
#' @description Simulate a population forwards with recombination that approximates the Structured Wright Fisher Process
#' and tracks haplotype identity by descent where individuals represent demes, such that within
#' a deme individual-level COI is considered.
#'
#' @param pos vector; the genomic coordinates for chromosome and position of the sites
#' @param N numeric; The number of individuals to consider
#' @param m numeric; Probability of migration where m represents the probability of moving from \italic{host_{origin}} to \italic{host_{new}} by \italic{m*(1-1/N)}
#' @param mean_coi numeric; The lambda of a right-shifted Poisson process, \italic{1 + Pos(lambda)} representing the average COI of an individual deme
#' @param rho numeric; expected recombination rate
#' @param tlim numeric; the maximum number of generations to consider before exitting gracefully if all samples have not coalesced
#' 
#' @return Returns a list of length six that contains the COI of each individual. A recombination list of length of tlim
#' where each element contains the recombination block -- as a boolean -- of the two parental haplotypes.   
#'  (the number of generations it took for all lineages to coalesce). Finally, there are lists for the parental host 
#'  and parental haplotype assignments for the "paternal" and "maternal" haplotypes (1 and 2), respectively. 
#'  
#' @details This function is intended to fed directly into the \code{get_arg} function for interpretability.   
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
               mean_coi = mean_coi,
               tlim = tlim,
               odd_prob = odd_prob)
  
  # run efficient C++ function
  ret <- sim_swf_cpp(args)
  
  # return as custom class
  class(ret) <-"swfsim"
  return(ret)
}

#------------------------------------------------
#' @title Get ancestral recombination graph from forward simulations
#'
#' @description Given an object \code{swf}, which is the result of forward
#'   simulation using the function \code{sim_swf()}, walks backwards through the
#'   ancestry and calculates the coalescent tree at every locus for the
#'   specified haplotypes.
#'
#' @param swf result of forwards simulation using the function \code{sim_swf()}
#' @param host_index a vector of target hosts. Defaults to all hosts
#' @param haplo_index a list of target haplotypes within the hosts specified by
#'   \code{host_index}. Defaults to all haplotypes within the specified hosts
#'
#' @export

get_arg <- function(swf, host_index = NULL, haplo_index = NULL) {
  
  # check inputs and define defaults
  assert_custom_class(swf, "swfsim")
  if (is.null(host_index)) {
    host_index <- 1:length(swf$coi)
  }
  if (is.null(haplo_index)) {
    haplo_index <- mapply(function(x) 1:x, swf$coi[host_index], SIMPLIFY = FALSE)
  }
  assert_vector(host_index)
  assert_pos_int(host_index, zero_allowed = FALSE)
  assert_list(haplo_index)
  assert_same_length(host_index, haplo_index)
  assert_pos_int(unlist(haplo_index), zero_allowed = FALSE)
  
  # define arguments
  args <- c(swf, list(host_index = rep(host_index, mapply(length, haplo_index)) - 1,
                      haplo_index = unlist(haplo_index) - 1))
  
  # pass to efficient C++ function
  output_raw <- get_arg_cpp(args)
  
  # create a list of bvtrees
  L <- length(output_raw$coalesce_target)
  ARG <- lapply(1:L, function(x) return(new("bvtree")))
  for (l in 1:L) {
    ARG[[l]]@c = output_raw$coalesce_target[[l]]
    ARG[[l]]@t = output_raw$coalesce_time[[l]]
    t_temp <- ARG[[l]]@t
    t_temp[which(t_temp == -1)] <- NA
    ARG[[l]]@z = order(t_temp) - 1
  }
  
  return(ARG)
}

#------------------------------------------------
#' @title Subset an object of class bvtree
#'
#' @description Given a bvtree and a vector of indices \code{s}, creates a new
#'   tree which is a subset of the original tree focusing only on the elements
#'   \code{s}.
#'
#' @param bvtree an object of class "bvtree"
#' @param s a vector specifying which elements in the bvtree to focus on
#'
#' @export

subset_bvtree <- function(bvtree, s) {
  
  # check inputs
  assert_custom_class(bvtree, "bvtree")
  assert_vector(s)
  assert_gr(length(s), 1)
  assert_noduplicates(s)
  assert_pos_int(s, zero_allowed = FALSE)
  assert_leq(s, length(bvtree@c))
  
  # create mask vector
  m <- rep(0, length(bvtree@c))
  m[s] <- TRUE
  
  # run efficient C++ function
  output_raw <- subset_bvtree_cpp(bvtree@c, bvtree@t, m)
  
  # update tree
  bvtree@c <- output_raw$c[s]
  bvtree@t <- output_raw$t[s]
  
  # update indices
  bvtree@c <- c(-1, (1:length(s))-1)[match(bvtree@c, c(-1,s-1))]
  
  # update z
  tmp <- bvtree@t
  tmp[tmp == -1] <- NA
  bvtree@z <- order(tmp) - 1
  
  return(bvtree)
}
