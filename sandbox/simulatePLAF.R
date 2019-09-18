
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



