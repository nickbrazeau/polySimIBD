#'
#' @noRd
#' Internal use. "Hazard" for recombination. K fixed at 1

prob.recombbo.k1 <- function(rho, d){return( 1 - exp(-rho * (d)) )}

#'
#' @noRd
#' Internal Use. Run out the prob of recombination between two
#' strains -- resulting IBD or not....

findrecombination <- function(chrompos, rho){

  recombo.block <- rep(NA, nrow(chrompos))

  # first state
  # randomly pick a strain to start with
  # then see if those two strains have a switch/recbomination
  # event as we move through genome
  recombo.block[1] <- sample( c("I","U"), size = 1)


  # draw subsequent states of recombination block
  for (i in 2:length(recombo.block)) {
    d <- chrompos$POS[i] - chrompos$POS[i-1]

    if(recombo.block[i-1] == "I"){
      recombo.block[i] <- sample( c("I","U"), size = 1, prob = c(1-prob.recombbo.k1(rho,d), prob.recombbo.k1(rho,d)) )
    } else {
      recombo.block[i] <- sample(c("U","I"), size = 1, prob = c(1-prob.recombbo.k1(rho,d), prob.recombbo.k1(rho,d)) )
    }
  }
  recombo.block
  return(recombo.block)
}


#' @noRd
#' Internal use
#' This function takes in two simhaplos and the recombination block
#' and makes haplobits and haplogt with the simulated crossover from
#' the recombination block
makecrossover <- function(smpl1, smpl2, recombo.block){
  out1 <- out2 <- new("simhaplo")

  if( length(unique(recombo.block)) != 1){
    # crossover sister chromatid 1
    out1@haplogt <- smpl1@haplogt
    out1@haplogt[recombo.block == "I"] <- smpl2@haplogt[recombo.block == "I"]
    out1@haplobit <- smpl1@haplobit
    out1@haplobit[recombo.block == "I"] <- smpl2@haplobit[recombo.block == "I"]

    # crossover sister chromatid 2
    out2@haplogt <- smpl2@haplogt
    out2@haplogt[recombo.block == "I"] <- smpl1@haplogt[recombo.block == "I"]
    out2@haplobit <- smpl2@haplobit
    out2@haplobit[recombo.block == "I"] <- smpl1@haplobit[recombo.block == "I"]
  } else { # no recombination occured
    out1@haplogt <- smpl1@haplogt
    out2@haplogt <- smpl2@haplogt

  }

  ret <- list(out1, out2)
  return(ret)

}

