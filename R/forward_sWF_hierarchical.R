#' @title The Structured Wright Fisher Model for IBD
#'
#' @description Simulate a population forwards with recombination that approximates the Structured Wright Fisher Process
#' and tracks haplotype identity by descent where individuals represent demes, such that within
#' a deme individual-level COI is considered
#'
#' @param pos dataframe <factor><numeric>; the genomic coordinates for chromosome and position of the sites
#' @param N numeric; The number of individuals to consider
#' @param m numeric; Probability of migration where m represents the probability of a parasite moving from one individual to another individual
#' @param mean_coi numeric; The lambda of a zero-truncated Poisson process, representing the average COI among inviduals
#' @param rho numeric; expected recombination rate
#' @param tlim numeric; the maximum number of generations to consider before exitting gracefully if all samples have not coalesced
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


sim_hierarchical_structured_WF <- function(pos, N, MM, m, rho, mean_coi, tlim){

  # assertions
  assert_vector(pos)
  assert_bounded(N, left = 1, right = Inf)
  assert_bounded(m, left = 0, right = 1)
  assert_bounded(rho, left = 0, right = 1, inclusive_left = F, inclusive_right = F)

  # warnings
  if(m == 0 & N > 1){
    warning("You have set the migration rate to 0 but have more than one individual As a result, all of your samples can never coalesce and this simulation will be limited by the tlim argument")
  }
  if(m==1){
    warning("With the migration rate set to 1, you have essentially created a panmictic population.")
  }


  # start
  L <- length(pos)

  # draw initial COIs
  coi_prev <- rpois(N, mean_coi) + 1

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
    coi <- rpois(N, mean_coi) + 1

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



