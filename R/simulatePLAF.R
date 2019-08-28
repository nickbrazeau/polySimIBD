
#' @title Simulate data ...todo...
#'
#' @description Simulate a population of samples from an pop AF distribution
#'
#' @param pos numeric vector; genomic coordinates
#' @param nsamps integer; number of samples
#' @param p_shape1 integer; alpha in beta dist
#' @param p_shape2 integer; beta in beta dist
#' @param propMissing integer; prob of missingness per sample
#'
#' @export

simPopData <- function(pos=c(sort(sample(1e5, 1e2))),
                       nsamps = 2,
                       p = NULL,
                       p_shape1 = 0.1,
                       p_shape2 = 0.1,
                       propMissing = 0) {


  n <- length(pos)
  # if p=NULL then simulate frequency of the REF allele at each locus in each contig
  if (is.null(p)) {
    p <- rbeta(n, p_shape1, p_shape2)
  }


  # initialise objects
  CHROMPOS <- data.frame(CHROM = "contig1", POS = pos)
  genmat <- matrix(NA, nrow = length(pos), ncol = nsamps)

  genmat <- apply(genmat, 2,
                  function(x){return(2*rbinom(n = n, size = 1, prob = 1 - p))})

  colnames(genmat) <- paste0("smp", 1:ncol(genmat))
  genmat <- cbind.data.frame(CHROMPOS, genmat)

  # missing data
  if (propMissing > 0) {
    for(i in 2:ncol(genmat)){
      genmat[sample(1:sum(n), round(sum(n)*propMissing)), i] <- -1
    }
  }
  ret <- list(paf = p,
              genmat = genmat)
  return(ret)

}



#' @title Simulate IBD sections between two samples under the hmmIBD model
#'
#' @description Walks along a vector of genomic locations and swiches between two states that represent IBD and non-IBD using a Markov model. The parameters that dictate the chance of switching state at any point include the average level of relatedness (\code{f}), the physical distance between loci in units of base pairs (from \code{pos}), and the recombination rate (\code{rho}), which is assumed constant over all loci.
#'
#' @param smpl1 numeric vector; haplotype vector
#' @param smpl2 numeric vector; haplotype vector
#' @param f the average relatedness between the two samples
#' @param rho the recombination rate
#' @param k the number of generations separating the two lineages
#' @param pos the genomic positions of the sites of interest
#'
#' @export

sim_hmmIBD <- function(smpl1, smpl2, f, k, rho, pos) {

  # draw starting state
  n <- length(pos) # abs number of loci, not pos -- this is consistent with simData below
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

  smpl1[ret == 1] <- smpl2[ret == 1]

  out <- list(
    simIBDvector = ret,
    smpl1 = smpl1,
    smpl2 = smpl2
  )
  return(out)
}


