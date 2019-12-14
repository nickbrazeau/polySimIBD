#' Layer Mutations onto the ARG for Each Loci
#' @param ARG set of bvtree; a set of bvtrees 
#' @param mutationrate numeric; the genome-wide per-generation mutation rate
#' 
#' @details The mutation model approximates an infinite allele model in time which is then collapsed into a single loci. 
#' Mutations are drawn with respect to the mutation rate and overall tree length from a poisson model. Mutations are
#' then "droppped" onto the tree following a uniform distribution. Mutations that happen upstream (i.e. are ancestral)
#' are carried along the tree branch to produce the final haplotypes for each parasite node.
#'
#'
#' @return hapmat numeric matrix; a matrix of mutliallelic haplotypes for each parasite considered. Loci are in
#' rows and parasites (haplotypes) are in columns. 
#' @export


layer_mutations_on_ARG <- function(mutationrate, ARG){
  
  # assertions
  assert_numeric(mutationrate)
  assert_custom_class(ARG[[1]], "bvtree", message  =  "Elements within the %s must inherit from class '%s'")
  
  # get haplotype matrix
  hapmat <- polySimIBD::get_haplotype_matrix(ARG = ARG)
  
  
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



