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

  hapmat <- matrix(NA, nrow = length(ARGsim$ARG), ncol = length(ARGsim$ARG[[1]]@c)) # loci by samples

  for (l in 1:length(ARGsim$ARG)){

    # extract objs
    coalbvtree <- ARGsim$ARG[[l]]
    pairw <- ARGsim$coal_times[[l]] %>%
      broom::tidy(.)

    # init gt vector for loci
    gt_l <- rep(0, length(coalbvtree@t))

    Tb <- sum(coalbvtree@t) + max(coalbvtree@t) + 1 # total branch length, max(t) to account for root length and then +1 to account for root -1 index
    mutn <- rpois(n = 1, lambda = mutationrate*Tb) # number of mutations
    if (mutn > 0){

      mutg <- rev( sort( runif(n=mutn, min = 0, max = Tb ) ) ) # need runif for infinite allele model
      # run back in time to present to allow for the overwriting of branches with multiple mutations

      for(i in 1:length(mutg)){
        # which branch to mutate
        blen <- coalbvtree@t
        blen[bp == -1] <- max(coalbvtree@t) # overwite wrote branch length

      }



    }

    # write to larger GTmat
    GTmat[l, ] <- gt_l

  }

  return(GTmat)
}



