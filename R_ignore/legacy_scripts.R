

#' @title Make child seuqence from two parent sequence, e.g. a chimeric sequence
#' @param p1 numceric;
#' @param p2 numeric;
#' @param rho numeric;
#' @param pos int vector; Genomic positions along the chromosome
#' @noRd
#' @noMd
#' @details Internal Function

recombine <- function(p1, p2, rho, pos) {
  
  # special case - no recombination
  L <- length(pos)
  if (identical(p1, p2)) {
    return(rep(p1, L))
  }
  
  # draw number of recombination events that occur over the observed genome
  n_breaks <- rpois(1, rho*pos[L])
  
  # special case if no recombination - all originate from one parent
  if (n_breaks == 0) {
    current_i <- sample(c(p1, p2), 1)
    ret <- rep(current_i, L)
    return(ret)
  }
  
  # draw parents over recombination blocks
  breakpoints <- sort(runif(n_breaks, 0, pos[L]))
  recombo_block <- (findInterval(pos, breakpoints) %% 2) + 1
  ret <- c(p1, p2)[recombo_block]
  return(ret)
}



#' Zero Truncated Poisson Model
#' @param n numeric; number of draws to consider
#' @param lambda numeric; the lambda parameter in a poisson distribtuion
#'
#' @return integer vector of sucesses
#' @export

zero_trunc_poisson <- function(n, lambda){
  # giocc.com/zero_truncated_poisson_sampling_algorithm.html
  ret <- rep(NA, n)
  for(i in 1:n){
    # do sampling
    k <- 1
    s <- t <- exp(-lambda)/(1-exp(-lambda)) * lambda
    u <- runif(1)
    
    while(s < u){
      k <- k + 1
      t <- t * lambda/k
      s <- s + t
    }
    ret[i] <- k
  }
  return(ret)
}



#' @title The Structured Wright Fisher Model for IBD
#'
#' @description Simulate a population forwards with recombination that approximates the Structured Wright Fisher Process
#' and tracks haplotype identity by descent where individuals represent demes, such that within
#' a deme individual-level COI is considered
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
#' Find Coalescence
#' @param swf S4 object; A discrete-loci, discrete-time structured Wright Fisher Simulation
#'            called from `sim_structured_WF`.
#' @param parasites int vector; Terminal nodes, or parasites from the Structured Wright Fisher
#'        Simulation for consideration.
#' @export

get_ARG <- function(swf, parasites = NULL){
  
  # graceful exits
  if(sum(swf$coi) == 1){
    warning("Your simulation returned one parasite among all hosts (did you run N = 1 and a low mean_moi?). As such, no pairwise comparisons can be made.")
    return("Only one parasite, no inference can be made")
  }
  
  # assertions
  assert_custom_class(x = swf, c = "sWFsim")
  assert_gr(length(parasites), 1, message = "Parasites must be a vector of terminal nodes for consideration that is longer than 1.")
  
  # tidy up params
  L <- length(swf$pos)
  anc <- swf$anc
  coi <- sum(swf$coi)
  
  if(is.null(parasites)){
    parasites <- 1:coi
  }
  
  
  ARG <- lapply(1:L, function(x) return(new("bvtree")))
  coaltime.squaremat <- lapply(1:L, function(x) return(matrix(NA, nrow = length(parasites), ncol = length(parasites))))
  
  #...................
  # get t and c
  #...................
  for(l in 1:L){ # for each loci find time that pair coalesces
    pairs <- as.data.frame( expand.grid(parasites, parasites) )
    pairs <- pairs[pairs[,1] != pairs[,2], ]
    pairs$coaltime <- rep(NA, nrow(pairs))
    for(i in 1:nrow(pairs)){ # find when each pair coalesces
      
      p1 <- pairs[i,1] # first "pointer"
      p2 <- pairs[i,2] # second "pointer"
      
      for (g in length(anc):1) {
        p1 <- anc[[g]][cbind(p1, l)]
        p2 <- anc[[g]][cbind(p2, l)]
        if(is.na(pairs$coaltime[i])){
          if(p1 == p2) {
            pairs$coaltime[i] <- length(anc)-g + 1 # 1-based
          }
        }
      } # end for loop for g
    } # end for loop for pairs
    
    # if pairs do not have a common ancestor, then set to tlim
    if(any(is.na(pairs$coaltime))){
      pairs$coaltime[is.na(pairs$coaltime)] <- -1
      warning("Your simulation did not fully coalesce. Pairs that did not reach a common ancestor have been set to the t-limit of -1.")
    }
    
    
    # extract lightweight infromation for bvtree
    # make connections for parasites always go right
    for(i in 1:(length(parasites)-1)){
      mintime <- min( pairs[ pairs$Var1 == parasites[i] & pairs$Var2 > i, ]$coaltime )
      
      # only do look ahead when we aren't one right of root
      if(i == length(parasites)-1){
        mintime.forward <- Inf
      } else {
        # loop through all pairs ahead
        mintime.ahead <- c()
        for(j in (i+1):(length(parasites)-1)){
          mintime.ahead <- c(mintime.ahead, min( pairs[ pairs$Var1 == parasites[j] & pairs$Var2 > j, ]$coaltime ))
        }
        mintime.ahead <- min(mintime.ahead)
      }
      
      lineage <- pairs$Var2[ pairs$Var1 == parasites[i] & pairs$coaltime == mintime & pairs$Var2 > i ] # always make trees look "right" via pairs$Var2 > i
      
      if(length(lineage) > 1){ # multiple coal
        if(mintime >= mintime.ahead){ # look ahead to make sure we don't block future branches
          lineage <- max(lineage) # pick farthest right
        } else {
          lineage <- lineage[1] # pick first (most left)
        }
      }
      ARG[[l]]@c[i] <- lineage
      ARG[[l]]@t[i] <- mintime
    }
    
    # note last node must always be root (if we are always going right)
    ARG[[l]]@c[length(parasites)] <- -1
    ARG[[l]]@t[length(parasites)] <- -1
    
    
    #...................
    # get z
    #...................
    ord <- sort(ARG[[l]]@t)
    ord.un <- unique(ord)
    ord.un <- ord.un[ord.un != -1] # remove root
    ord.num <- 1:length(ord.un)
    
    ARG[[l]]@z <-  ARG[[l]]@t # temporary overwrite
    
    for(i in 1:length(ord.un)){
      ARG[[l]]@z[ ARG[[l]]@z == ord.un[i] ] <- ord.num[i]
    }
    # end section for Z fill in
    
    #...................
    # get coal time square matrix
    #...................
    coaltime.squaremat[[l]] <-  pairs %>%
      tidyr::spread(., key = "Var2", value = "coaltime") %>%
      dplyr::select(-c("Var1")) %>%
      stats::as.dist()
    
  } #end for loop for loci
  
  
  ARGlist <- list(ARG = ARG,
                  coal_times = coaltime.squaremat,
                  coi = swf$coi,
                  parasites = parasites)
  
  class(ARGlist) <- "ARGsim"
  
  return(ARGlist)
}




