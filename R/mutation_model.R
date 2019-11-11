#' Layer Mutations onto the ARG for Each Loci
#' @param ARGsim S4 object;
#' @param mutationrate numeric; the genome-wide per-generation mutation rate
#'
#'
#' @return GTmatrix numeric matrix;
#' @export


layer_mutations_on_ARG <- function(mutationrate, ARGsim){
  # assetions
  assert_custom_class(x = ARGsim, c = "ARGsim")
  assert_numeric(mutationrate)

  # init
  hapmat <- matrix(NA, nrow = length(ARGsim$ARG), ncol = length(ARGsim$ARG[[1]]@c)) # loci by samples

  # find breakpoints
  # going to assume that mutational rate among recombination blocks is constant
  brkpts <- c()
  for (l in 2:length(ARGsim$ARG)) {
    if (! isTRUE(all.equal(ARGsim$ARG[[l-1]], ARGsim$ARG[[l]]))) { # if breakpoint between loci
      brkpts <- c(brkpts, l)
    }
  }

  # mutate trees for given breakpoint
  for (t in c(1, brkpts)) {
    # extract objs
    coalbvtree <- ARGsim$ARG[[t]]
    pairw <- ARGsim$coal_times[[t]] %>%
      broom::tidy(.)

    # init gt vector for loci
    gt_l <- rep(0, length(coalbvtree@t))

    # branch lengths, make root equal to max
    brnchlngth <- coalbvtree@t
    brnchlngth[brnchlngth == -1] <- max(brnchlngth)

    Tbrnchlngth <- sum(brnchlngth) # total branch length
    mutn <- rpois(n = 1, lambda = mutationrate*Tbrnchlngth) # number of mutations
    if (mutn > 0){

      mutg <-  sort( runif(n=mutn, min = 0, max = Tbrnchlngth ) )  # need runif for infinite allele model
      # run back in time to present to allow for the overwriting of branches with multiple mutations

      for (i in 1:length(mutg)) {
        # what lineage does this occur on
        mutlin <- min( which(  mutg[i] < cumsum(brnchlngth) ))
        # are any other lineages affected
        for (j in 1:nrow(pairw)) {
          # find nodes that are descendants of the branch that is currently being mutated
          # Note, item2 due to the fact that the swf object has cols 1:(n-1) and rows, 2:n
          mutlindesc <- pairw$item1[ (pairw$distance < mutg[i] & pairw$item2 == mutlin) ]
        }
        mutlin <- c(mutlin, mutlindesc)
        gt_l[ mutlin ] <- i
      } # end if for mutations
    } # end for loop for mutations

    # fill in hap matrix
    if (is.null(brkpts)) { # same tree throughout
      hapmat[1:nrow(hapmat), ] <- gt_l
    } else if (max(brkpts) == t) {
      hapmat[t:nrow(hapmat), ] <- gt_l
    } else {
      hapmat[t: brkpts[ which(t == c(1, brkpts)) ], ] <- gt_l # note, we temporarily include beginning of next block, is immediately overwritten in next gen
    }
  }

  # liftover haplotype matrix so that we only carry unique ints that start at 0, 1, 2 ...
  hapmat <- t( apply(hapmat, 1, function(x){
    unq <- unique(x)
    for(i in 1:length(unq)){
      x[x == unq[i]] <- i - 1 # zero-based
    }
    return(x)
  }) )

  # return
  return(hapmat)

}



