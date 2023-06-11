# COI_equilibrium.R
#
# Author: Bob Verity
# Date: 2022-08-16
#
# Purpose:
# Derive the equilibrium COI distribution under the PolySimIBD model for
# infinite population size.
#
# ------------------------------------------------------------------

library(copula) # needed for Stirling numbers

# define model parameters
lambda <- 3
m <- 0.1
K_max <- 20

# calculate transition matrix
K_mat <- matrix(NA, K_max, K_max)
for (K_prev in 1:K_max) {
  
  p2 <- rep(0, K_max)
  for (i in 1:K_max) {
    N_clonal <- i
    k <- 1:min(N_clonal, K_prev)
    p1 <- exp(lfactorial(K_prev) - lfactorial(K_prev - k) - N_clonal*log(K_prev)) * Stirling2.all(N_clonal)[k]
    p2 <- p2 + dpois(i, lambda * (1 - m) / K_prev) * c(p1, rep(0, K_max - length(p1)))
  }
  p2 <- c(dpois(0, lambda * (1 - m) / K_prev), p2)
  
  p3 <- rep(NA, K_max + 1)
  lambda2 <- lambda * m + lambda * (1 - m) * (1 - 1 / K_prev)
  for (i in 0:K_max) {
    p3[i + 1] <- sum(p2[(0:i) + 1] * dpois(i:0, lambda2))
  }
  
  K_mat[K_prev,] <- p3[-1]
}
K_mat <- K_mat / (1 - dpois(0, lambda))

# calculate equilibrium solution by applying transition a large number of times
# (I'm having trouble getting the Eigenvalue method to work)
K_init <- rep(1/K_max, K_max)
for (i in 1:100) {
  K_init <- K_init %*% K_mat
}

# barplot of equilibrium distribution. Red dots show the raw Poisson(lambda)
# distribtion, i.e. what we would expect under pure super-infection from an
# infinitely large population
barplot(K_init, space = 0, names.arg = 1:K_max, col = grey(0.6))
points(1:K_max - 0.5, dpois(1:K_max, lambda) / (1 - exp(-lambda)), pch = 20, col = 2)
