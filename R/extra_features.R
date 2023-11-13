
#------------------------------------------------
#' @title Extract haplotypes from ARG
#' @param arg set of bvtrees
#' @return hapmat numeric matrix; a matrix of multiallelic haplotypes for each parasite considered. Loci are in
#' rows and parasites (haplotypes) are in columns. 
#' @export
extract_haplotype_matrix <- function(arg){
  
  # convert trees into matrix of alleles
  # each column is therefore a haplotype since we consider parasite by parasite
  hap_mat <- t(mapply(function(x) {
    c <- x@c
    ret <- c
    ret[ret == -1] <- 1:sum(ret == -1)
    while (any(c != -1)) {
      w <- which(c == -1)
      c[-w] <- c[c[-w]+1]
      ret[-w] <- ret[ret[-w]+1]
    }
    return(ret)
  }, ARG))
  return(hap_mat)
}

#-------------------------------------------------------------------------------------------
#' @title Layer Mutations onto the ARG for Each Loci
#' @inheritParams extract_haplotype_matrix 
#' @param mutationrate numeric; the genome-wide per-generation mutation rate
#' @description The mutation model approximates an infinite allele model in time which is then collapsed into a single loci. 
#' Mutations are drawn with respect to the mutation rate and overall tree length from a poisson model. Mutations are
#' then "droppped" onto the tree following a uniform distribution. Mutations that happen upstream (i.e. are ancestral)
#' are carried along the tree branch to produce the final haplotypes for each parasite node.
#' @return hapmat numeric matrix; a matrix of mutliallelic haplotypes for each parasite considered. Loci are in
#' rows and parasites (haplotypes) are in columns. 
#' @importFrom stats rpois runif
#' @export


layer_mutations_on_ARG <- function(arg, mutationrate){

  # assertions
  goodegg::assert_numeric(mutationrate)
  goodegg::assert_class(ARG[[1]], "bvtree", message  =  "Elements within the %s must inherit from class '%s'")
  
  # get haplotype matrix
  hapmat <- polySimIBD::extract_haplotype_matrix(arg = ARG)
  
  
  # NB, to keep alleles unique, we will start counting after the point of the genealogy alleles
  mut_allelestart <- max(hapmat)
  
  # mutate trees for each loci
  for (t in 1:length(ARG)) { 
    # extract objs
    coalbvtree <- ARG[[t]]
    
    # init hapmat and gt vector for loci
    gt_l <- rep(NA, length(coalbvtree@t))
    
    # branch lengths, make root equal to max
    brnchlngth <- coalbvtree@t
    
    # catch if no coalescence (all roots)
    if (max(brnchlngth) == -1 ) {
      alltreesbrnchlngth <- unlist( purrr::map(ARG, "t") )
      brnchlngth[brnchlngth == -1] <- max(alltreesbrnchlngth) # just set to max time in simulation
    } else {
      brnchlngth[brnchlngth == -1] <- max(brnchlngth) + 1e-5 # plus small bit 
    }
    
    Tbrnchlngth <- sum(brnchlngth) # total branch length
    mutn <- rpois(n = 1, lambda = mutationrate*Tbrnchlngth) # number of mutations
    
    if (mutn > 0) {
      
      mutg <-  sort( runif(n=mutn, min = 0, max = Tbrnchlngth ), decreasing = T )  # need runif for infinite allele model
      # run back in time to present to allow for the overwriting of branches with multiple mutations
      
      # keep track of orig position in vec but sort for easier looping through mutations
      names(brnchlngth) <- 1:length(brnchlngth)
      sorted.coaltimes <- sort(brnchlngth)
      
      for (i in 1:length(mutg)) {
        # what lineage does this occur on
        mutlin.sorted <-   sort( which( cumsum(sorted.coaltimes) > mutg[i] ) ) # get min, but use sort to keep names 
        mutlin <- which(names(brnchlngth) == names(mutlin.sorted)[1])
        
        # does this affect any other lineages
        if(any(coalbvtree@t[mutlin] > coalbvtree@t[coalbvtree@t != -1])){ # does this coalesce later than other samples
          if(any(coalbvtree@c == mutlin - 1)){ # does this sample have connections; -1 for cpp to R 
            mut_brnchlngth <- coalbvtree@t
            mut_brnchlngth[mut_brnchlngth == -1] <- .Machine$double.xmax
            extmutlin <- which(coalbvtree@t[mutlin] > mut_brnchlngth &
                                 coalbvtree@c == mutlin - 1)
            mutlin <- c(mutlin, extmutlin)
          }
        } 
        
        gt_l[ mutlin ] <- i
      } # end for loop for mutations
      
      
      # fill in hap matrix
      hapmat[t, !is.na(gt_l)] <- gt_l[!is.na(gt_l)] + mut_allelestart
      
      
      # liftover haplotype matrix so that we only carry unique ints and always start at 0 (and then count up 1, 2 ...)
      hapmat <- t( apply(hapmat, 1, function(x){
        unq <- unique(x)
        for(i in 1:length(unq)){
          x[x == unq[i]] <- i - 1 # zero-based
        }
        return(x)
      }) )
      
    } 
    
  } # end for loop for trees (i.e. a tree per loci)
  
  # return
  return(hapmat)
  
}



#-------------------------------------------------------------------------------------------
#' @title Simulate biallelic data
#'
#' @description Simulate biallelic data from a simple statistical model. Inputs
#'   include the complexity of infection and some parameters dictating skew and error
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
#'     \deqn{wsaf_{error} = wsaf*(1-e) + (1-wsaf)*e} where \eqn{wsaf} is the WSAF
#'     without error and \eqn{e} is the error parameter \code{epsilon}.
#'     \item Final read counts are drawn from a beta-binomial distribution with
#'     expectation \eqn{w_{error}}. The raw number of draws is given by the
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
#' @importFrom stats rbeta
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
  goodegg::assert_matrix(haplotypematrix)
  goodegg::assert_vector(coverage)
  goodegg::assert_pos_int(coverage)
  goodegg::assert_single_pos(alpha, zero_allowed = FALSE)
  goodegg::assert_single_pos(overdispersion, zero_allowed = TRUE)
  goodegg::assert_single_pos(epsilon, zero_allowed = TRUE)
  goodegg::assert_bounded(epsilon)
  goodegg::assert_eq(x = sum(COIs), y = ncol(haplotypematrix),
            message = "The COIsum must be equal to the number of columns in your haplotype matrix")
  
  
  
  # make copy of hapmat because we are going to
  # overwrite it with the alleles that we draw to be biallelic
  m <- haplotypematrix
  
  # generate biallelic table from PLAF of haplotype
  PLAF <- rbeta(n = nrow(m), shape1 = shape1, shape2 = shape2)
  
  for(i in 1:nrow(m)){
    uniqueAllele <- unique( m[i, ] )
    liftoverAlleles <- sample(x = c(1,0), 
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
