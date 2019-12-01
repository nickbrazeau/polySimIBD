
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
#' @description if mean_coi > 10, then we say the zero-trunc approximates the pois
#'
#'@details Data is simulated under the following framework:
#'   \enumerate{
#'        \item ... TODO>..
#'   }
#'
#' @return Returns a list of length G (the number of generations it took for all lineages to coalesce). Each
#' list item is a matrix with n x L dimensions where n is the number of parasites considered and L is loci.
#' Parasites are then indexed into demes (individuals) based on the COI distribution.
#'
#' @examples
#' \dontrun{
#'
#' }

#'
#' @export


sim_structured_WF <- function(pos, N, m, rho, mean_coi, tlim){

  # assertions
  assert_vector(pos)
  assert_bounded(N, left = 1, right = Inf)
  assert_bounded(m, left = 0, right = 1)
  assert_bounded(rho, left = 0, right = 1, inclusive_left = F, inclusive_right = F)

  # warnings
  if(m == 0 & N > 1){
    warning("You have set the migration rate to 0 but have more than one deme. As a result, all of your samples can never coalesce and this simulation will be limited by the tlim argument")
  }
  if(m==1){
    warning("With the migration rate set to 1, you have essentially created a panmictic population.")
  }


  # start
  L <- length(pos)

  # draw initial COIs
  if(mean_coi > 10){
    coi_prev <- rpois(N, mean_coi)
  } else {
    coi_prev <- zero_trunc_poisson(N, mean_coi)
  }

  # create ancestry array
  anc <- list()
  # create haploint matrix
  haploint_prev <- outer(1:sum(coi_prev), rep(1,L))

  # loop through generations
  g <- 1
  uniques <- 2
  iter = 0
  while (uniques != 1 & iter != tlim) {
    g <- g + 1

    # draw COIs
    coi <- zero_trunc_poisson(N, mean_coi)

    # initialize haploints
    haploint <- matrix(NA, sum(coi), L)

    # expand ancestry list
    anc[[g-1]] <- matrix(NA, sum(coi), L)

    # loop through demes and haplotypes within demes
    i <- 0
    for (n in 1:N) {
      for (j in 1:coi[n]) {
        i <- i + 1

        # choose parental demes and parental haplotypes from those demes
        parent_deme1 <- ifelse(rbinom(1, 1, m), sample(1:N, 1), n)
        parent_haplo1 <- sample(coi_prev[parent_deme1], 1)
        parent_index1 <- sum(coi_prev[1:parent_deme1]) - coi_prev[parent_deme1] + parent_haplo1

        parent_deme2 <- ifelse(rbinom(1, 1, m), sample(1:N, 1), n)
        parent_haplo2 <- sample(coi_prev[parent_deme2], 1)
        parent_index2 <- sum(coi_prev[1:parent_deme2]) - coi_prev[parent_deme2] + parent_haplo2

        # apply recombination
        recombo_block <- recombine(parent_index1, parent_index2, rho, pos)
        anc[[g-1]][i,] <- recombo_block
        haploint[i,] <- haploint_prev[cbind(recombo_block, 1:L)]
      }
    }

    # update the haploint generation indices and the coi
    # to draw for the next generation
    haploint_prev <- haploint
    coi_prev <- coi

    # check if we can break loop by get number of unique haploints
    uniques <- length(unique(split(haploint, 1:nrow(haploint))))

    # check if we should exit gracefully
    iter <- iter + 1
  }

  # return out

  ret <- list(coi = coi, anc = anc, pos=pos)
  class(ret) <-"sWFsim"
  return(ret)


}

#------------------------------------------------
#' @title Run forwards simulation from structured Wright-Fisher model
#'
#' @description Run forwards simulation from structured Wright-Fisher model, in
#'   which "demes" correspond to human hosts, and fluctuating population sizes
#'   over generations correspond to changing COIs of individuals.
#'
#' @param pos genomic positions of loci to track
#' @param N number of hosts
#' @param m migration parameter, dictating the balance between co-transmission
#'   (small m) and super-infection (large m)
#' @param rho rate of recombination, in bp per generation
#' @param mean_coi the COI of every host is drawn from a zero-truncated Poisson
#'   distribution with mean \code{mean_coi/(1 -exp(-mean_coi))}
#' @param tlim run for this many generations
#'
#' @export

sim_swf <- function(pos, N, m, rho, mean_coi, tlim) {

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
    ARG[[l]]@z = order(ARG[[l]]@t, decreasing = TRUE) - 1
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
