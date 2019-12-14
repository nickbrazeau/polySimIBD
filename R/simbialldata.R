#------------------------------------------------
#' @title Simulate biallelic data
#'
#' @description Simulate biallelic data from a simple statistical model. Inputs
#'   include the complexity of infection (COI) and some parameters dicating skew and error
#'   distributions. Outputs include the phased haplotypes and the un-phased read
#'   count and coverage data.
#'
#' @details Simulated data are drawn from a simple statistical model:
#'   \enumerate{
#'     \item Strain proportions are drawn from a symmetric Dirichlet
#'     distribution with shape parameter \code{alpha}.
#'     \item Phased haplotypes are drawn at every locus, one for each
#'     \code{COI}. The allele at each locus is drawn from a Bernoulli
#'     distribution with probability given by the \code{PLAF}.
#'     \item The "true" within-sample allele frequency at every locus is
#'     obtained by multiplying haplotypes by their strain proportions, and
#'     summing over haplotypes. Errors are introduced through the equation
#'     \deqn{wsaf_error = wsaf*(1-e) + (1-wsaf)*e}where \eqn{wsaf} is the WSAF
#'     without error and \eqn{e} is the error parameter \code{epsilon}.
#'     \item Final read counts are drawn from a beta-binomial distribution with
#'     expectation \eqn{w_error}. The raw number of draws is given by the
#'     \code{coverage}, and the skew of the distribution is given by the
#'     \code{overdispersion} parameter. If \code{overdispersion = 0} then the
#'     distribution is binomial, rather than beta-binomial.
#'   }
#'
#' @param COIs integer vector of the number of haplotypes within each sample.
#' @param haplotypematrix matrix of haplotypes that correspond to within individual COI.
#' @param shape1 the alpha value for the beta binomial distribution, which the population-allele frequency
#' for each allele is drawn from.
#' @param shape2 the beta value for the beta binomial distribution, which the population-allele frequency
#' for each allele is drawn from.
#' @param coverage coverage at each locus. If a single value then the same
#'   coverage is applied over all loci.
#' @param alpha shape parameter of the symmetric Dirichlet prior on strain
#'   proportions.
#' @param overdispersion the extent to which counts are over-dispersed relative
#'   to the binomial distribution. Counts are Beta-binomially distributed, with
#'   the beta distribution having shape parameters \code{p/overdispersion} and
#'   \code{(1-p)/overdispersion}.
#' @param epsilon the probability of a single read being mis-called as the other
#'   allele. Applies in both directions.
#'
#' @return List of non-referent within-sample allele frequency dataframe (unphased) and phased vectors
#'         of non-referent within-sample allele counts, overall coverage, strain proportions, and the biallelic haplotype matrix.
#' @export

sim_biallelic <- function(COIs = c(1,1),
                          haplotypematrix = matrix(1, 2, 2),
                          shape1 = 1, 
                          shape2  = 1, 
                          coverage = 100,
                          alpha = 1,
                          overdispersion = 0,
                          epsilon = 0) {

  # check inputs
  if (length(coverage) == 1) {
    coverage <- rep(coverage, nrow(haplotypematrix))
  }
  assert_matrix(haplotypematrix)
  assert_vector(coverage)
  assert_pos_int(coverage)
  assert_single_pos(alpha, zero_allowed = FALSE)
  assert_single_pos(overdispersion, zero_allowed = TRUE)
  assert_single_pos(epsilon, zero_allowed = TRUE)
  assert_bounded(epsilon)
  assert_eq(x = sum(COIs), y = ncol(haplotypematrix),
            message = "The COIsum must be equal to the number of columns in your haplotype matrix")


  
  # make copy of hapmat because we are going to
  # overwrite it with the alleles that we draw to be biallelic
  m <- haplotypematrix
  
  # generate biallelic table from PLAF of haplotype
  PLAF <- rbeta(n = nrow(m), shape1 = shape1, shape2 = shape2)
  
  for(i in 1:nrow(m)){
    uniqueAllele <- unique( m[i, ] )
    liftoverAlleles <- sample(x = c(0,1), 
                              size = length(uniqueAllele),
                              prob = c(PLAF[i], (1-PLAF[i])), 
                              replace = T)
    names(liftoverAlleles) <- uniqueAllele
    
    for (j in 1:length(liftoverAlleles)) {
      # now convert multiallelic to biallelic
      m[i,][ m[i,] == names(liftoverAlleles[j]) ] <- liftoverAlleles[j]
    }
  }
  
  # split the haplotype matrix into individual (host) matrices 
  splitter <- rep(1:length(COIs), times = COIs)
  hosts.haplotypes <- NULL
  for (i in 1:length(unique(splitter))) {
    hosthap <- m[, c( splitter == i ), drop = F]
    hosts.haplotypes <- c(hosts.haplotypes, list(hosthap))
  }
  
  # generate strain proportions by drawing from dirichlet
  w <- polySimIBD::rdirichlet(rep(alpha, sum(COIs)))
  # make this a list that corresponds to within host COIs from above
  w.list <- split(w, factor(splitter))
  w.list <- lapply(w.list, function(x){x/sum(unlist(x))})


  # true WSAF levels by summing binomial draws over strain proportions
  get_wsaf <- function(x, y){
    if(length(y) == 1){
      ret <- x * y
    } else {
      ret <- rowSums(sweep(x, 2, y, "*"))
    }
  }
  host.wsaf <- mapply(get_wsaf, hosts.haplotypes, w.list)
  
  # add in gentoyping error
  # Note, genotyping error is fixed
  host.wsaf.genotypeerror <- host.wsaf * (1-epsilon) + (1-host.wsaf)*epsilon
  

  # draw read counts, taking into account overdispersion
  get_read_counts <- function(L, coverage, p, overdispersion){
    if (overdispersion == 0) {
      counts <- rbinom(L, size = coverage, prob = p)
    } else {
      counts <- rbetabinom(L, k = coverage, alpha = p/overdispersion, 
                           beta = (1-p)/overdispersion)
    }
    return(counts)
  }

  # get counts
  counts <- apply(host.wsaf.genotypeerror, 2, get_read_counts, 
                  L = nrow(haplotypematrix), 
                  coverage = coverage, 
                  overdispersion = overdispersion)

  # split into samples
  smpls.list <- lapply(splitter, function(x){
    return(counts[,x])
  })


  # convert counts and coverage to NRWSAF
  coverage <- matrix( rep(coverage, times = length(COIs)), ncol = length(COIs) )
  NRWSAF <- counts/coverage

  # return list
  ret <- list(rbetaPLAF = PLAF,
              strain_proportions = w.list,
              hosts.haplotypes = hosts.haplotypes,
              NRWSAcounts = counts,
              WS.coverage = coverage,
              NRWSAF = NRWSAF)

  return(ret)
}
