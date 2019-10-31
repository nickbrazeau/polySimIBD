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

