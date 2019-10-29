#' Layer Mutations onto the ARG for Each Loci
#' @param ARGsim S4 object;
#' @param mutationrate numeric; the genome-wide per-generation mutation rate
#'
#' @importFrom magrittr %>%
#'
#' @return GTmatrix numeric matrix;
#' @export


layer_mutations_on_ARG <- function(mutationrate, ARGsim){
  # assetions
  assert_custom_class(x = ARGsim, c = "ARGsim")
  assert_numeric(mutationrate)

  GTmat <- matrix(NA, nrow = length(ARGsim$ARG), ncol = length(ARGsim$ARG[[1]]@c)) # loci by samples

  for (l in 1:length(ARGsim$ARG)){

    # extract objs
    coalbvtree <- ARGsim$ARG[[l]]
    pairw <- ARGsim$coal_times[[l]] %>%
      broom::tidy(.)

    # init gt vector for loci
    gt_l <- rep(0, length(coalbvtree@t))

    Tb <- sum(coalbvtree@t)+1 # total branch length, +1 to account for root
    mutn <- rpois(n = 1, lambda = mutationrate*Tb) # number of mutations
    if (mutn > 0){

      mutg <- rev( sort( runif(n=mutn, min = 0, max = max(coalbvtree@t)) ) ) # need runif for infinite allele model
      # liftover root
      coalbvtree@t[length(coalbvtree@t)] <- .Machine$integer.max


      for (i in 1:mutn){
        # catch node mutation
        if (min(coalbvtree@t) > mutg[i]){
          # mutate a node
          b <- sample(1:length(coalbvtree@t), size = 1)
          gt_l[b] <- length(unique(gt_l)) + 1

        } else {
          b <- which(coalbvtree@t < mutg[i]) # branches to mutate
          b <- sample(b, size = 1)

          # downstream branch effects, loop through tree always going left
          nb <- pairw$item1[ (pairw$distance < mutg[i] & pairw$item2 == b) ] # find nodes that are descendants of the branch that is currently being mutated. Note, item2 due to the fact that the swf object has cols 1:(n-1) and rows, 2:n
          b <- c(b, nb)

          # mutate descendant nodes
          gt_l[b] <- length(unique(gt_l)) + 1

        }
      }

      # because gt_l can be written over many times as we move down the tree, let's restore these ints to base ints
      gt_l <- as.numeric(factor(gt_l))

    }

    # write to larger GTmat
    GTmat[l, ] <- gt_l

  }

  return(GTmat)
}



