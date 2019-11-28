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
#' @title TODO
#'
#' @description TODO.
#'
#' @param x TODO
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
