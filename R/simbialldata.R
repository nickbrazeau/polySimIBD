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
                          coverage = 100,
                          alpha = 1,
                          overdispersion = 0,
                          epsilon = 0) {

  # get loci
  L <- nrow(haplotypematrix)

  # check inputs
  if (length(coverage) == 1) {
    coverage <- rep(coverage, L)
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



  # generate biallelic table from PLAF of haplotype
  m <- haplotypematrix
  for(i in 1:nrow(m)){
    uniqueAllele <- unique( m[i, ] )

    if (!all(uniqueAllele %in% c(0,1))) {
      # multiallelic (multiallelic to biallelic 50-50 for now)
      liftoverAlleles <- sample(x = c(0,1), size = length(uniqueAllele[! uniqueAllele %in% c(0,1)]), replace = T)
      names(liftoverAlleles) <- uniqueAllele[! uniqueAllele %in% c(0,1)]
    }
    for (j in 1:length(liftoverAlleles)) {
      # now convert multiallelic to biallelic
      m[i,][ m[i,] == names(liftoverAlleles[j]) ] <- liftoverAlleles[j]
    }
  }

  # get COI matrices
  start <- c(1, (COIs+1)[1:(length(COIs)-1)])
  end <- cumsum(COIs)
  COIsplitter <- cbind(start, end)
  COIsplitter <- apply(COIsplitter, 1, function(x){return(x[1]:x[2])})
  smpls.list <- lapply(COIsplitter, function(x){
    return(m[,x])
  })

  # generate strain proportions
  w <- polySimIBD::rdirichlet(rep(alpha, sum(COIs)))
  w.list <- lapply(COIsplitter, function(x){
    return(w[x])
  })


  # true WSAF levels by summing binomial draws over strain proportions
  get_p_levels <- function(m, w, L){

    if(length(m) == L) { # COI of 1
      p_levels <- m * w
    } else {
      p_levels <- rowSums(sweep(m, 2, w, "*"))

    }

    return(p_levels)
  }

  p_levels.list <- mapply(m = smpls.list, w = w.list, L = L, get_p_levels)

  # add in genotyping error
  addgenotype_error <- function(p_levels, epsilon){
    p_error <- p_levels*(1-epsilon) + (1-p_levels)*epsilon
    return(p_error)
  }

  p_error.list <- mapply(p_levels = p_levels.list, epsilon = epsilon, addgenotype_error)

  # draw read counts, taking into account overdispersion
  get_read_counts <- function(L, coverage, p_error, overdispersion){
    if (overdispersion == 0) {
      counts <- rbinom(L, size = coverage, prob = p_error)
    } else {
      counts <- rbetabinom(L, k = coverage, alpha = p_error/overdispersion, beta = (1-p_error)/overdispersion)
    }
    return(counts)

  }

  # get counts
  counts <- mapply(p_error = p_error.list, L = L, coverage = coverage, overdispersion = overdispersion,
                   get_read_counts)

  # split into samples
  smpls.list <- lapply(COIsplitter, function(x){
    return(counts[,x])
  })


  # convert counts and coverage to NRWSAF
  NRWSAF <- matrix(NA, nrow = L, ncol = length(COIs))
  # anon functions to summarize counts and divide by coverage to get NSWRAF
  NRWSAF <- sapply(smpls.list,
                   function(x, coverage){
                     if (is.null(dim(x))) {
                       ret <- x/coverage
                     } else {
                       ret <- rowSums(x)/(coverage*ncol(x))
                     }
                     return(as.matrix(ret))
                   },
                   coverage = coverage)

  # save out positions
  NRWSAFdf <- cbind.data.frame(POS = pos, NRWSAF)
  colnames(NRWSAFdf)[2:ncol(NRWSAFdf)] <- paste0("smpl", 1:(ncol(NRWSAFdf)-1))


  # return list
  ret <- list(strain_proportions = w,
              haplotypematrix.biall = m,
              NRWSAcounts = counts,
              WS.coverage = coverage,
              NRWSAFdf = NRWSAFdf)

  return(ret)
}
