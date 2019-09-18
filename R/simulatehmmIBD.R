#' @title Simulate IBD sections between two samples under the hmmIBD model
#'
#' @description Walks along a vector of genomic locations and
#' swiches between two states that represent IBD and non-IBD using a Markov model.
#' The parameters that dictate the chance of switching state at any point include
#' the average level of relatedness (\code{f}), the physical distance between
#' loci in units of base pairs (from \code{pos}), and the recombination rate (\code{rho}),
#' which is assumed constant over all loci.
#'
#' @param chrompos matrix <factor><numeric>; the genomic coordinates for chromosome and position of the sites
#' @param PLAF numeric vector; the population-level allele frequencies to simulate genotypes from
#' @param f the average relatedness between the two samples
#' @param rho the recombination rate
#' @param k the number of generations separating the two lineages
#' @param pos the genomic positions of the sites of interest
#'
#'#'@details Data is simulated under the following framework:
#'   \enumerate{
#'        \item Draw two genotypes from the \code{PLAF}
#'        \item Encode haplotypes with a \emph{bit}
#'        \item Cross samples under the hmmIBD framework
#'   }
#'
#' @export

sim_hmmIBD <- function(chrompos, PLAF, f, k, rho, pos) {

  #.....................
  # Draw inds
  #.....................
  p1 <- new("simhaplo")
  p2 <- new("simhaplo")
  p1@haploint <- rep(1, nrow(chrompos))
  p2@haploint <- rep(2, nrow(chrompos))

  #.....................
  # Draw Recombination block
  #.....................

  # draw starting state
  n <- nrow(chrompos)
  ret <- rep(NA, n)
  ret[1] <- sample( c(0,1), size = 1, prob = c(1-f,f) )

  # draw subsequent states
  for (i in 2:n) {
    d <- pos[i]-pos[i-1]
    if (ret[i-1] == 0) {    # move from non-IBD state to non-IBD is t11
      t11 <- 1 - f*( 1 - exp( -k*rho*d ) )
      ret[i] <- sample( c(0,1), size = 1, prob = c(t11, 1 - t11) )
    } else {  # move from IBD state to IBD is t22
      t22 <- 1 - (1 - f)*( 1 - exp( -k*rho*d ) )
      ret[i] <- sample( c(0,1), size = 1, prob = c(1 - t22, t22) )
    }
  }

  #.....................
  # Apply Recombination block
  #.....................
  child <- p1
  child@haploint[ret == 1] <- p2@haploint[ret == 1]

  ret <- list(
    parent1 = p1,
    parent2 = p2,
    child = child
  )

  return(ret)
}


