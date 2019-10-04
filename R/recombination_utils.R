
#' @title Make child seuqence from two parent sequence, e.g. a chimeric sequence
#' @param p1 index
#' @param p2 index
#' @param chrompos
#' @param rho
#' @param interference numeric; probability that double crossover occurs
#' @noRd
#' @noMd
#' @details Internal Function

recombine <- function(p1, p2, rho, pos) {

  # special case - no recombination
  L <- length(pos)
  if (identical(p1, p2)) {
    return(rep(p1, L))
  }

  # draw number of recombination events that occur over the genome
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









