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

rpois(n=100, lambda = 1e-4)
zero_trunc_poisson(n=100, 1e-4)
