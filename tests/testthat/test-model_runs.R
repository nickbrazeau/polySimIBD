test_that("model will run", {
  N <- c(3, 2, 7)
  m <- c(0.5, 0.5, 0.5)
  mean_coi <- c(1, 1, 1)
  #migr_dist_mat <- matrix(runif(9), 3, 3)
  migr_dist_mat <- matrix(0, 3, 3)
  diag(migr_dist_mat) <- 1
  swf <- polySimIBD::sim_swf(pos = sample(1:1e3, size = 50), 
                             N = N, 
                             m = m, 
                             mean_coi = mean_coi, 
                             migr_dist_mat = migr_dist_mat, 
                             dist_scalar = 1,
                             rho = 1e-2, 
                             tlim = 5,
                             verbose = FALSE)
  
  ARG <- polySimIBD::get_arg(swf = swf, host_index = c(1,2))
  
  
})

ret$coi
ret$parent_host1[[5]]

unlist(ret$parent_host1[[5]])
unlist(ret$parent_host2[[5]])
unlist(ret$parent_haplo1[[5]])
unlist(ret$parent_haplo2[[5]])



swf <- polySimIBD::sim_swf(pos = sample(1:1e3, size = 50), 
                           N = 20, 
                           m = 0.5, 
                           mean_coi = 1, 
                           migr_dist_mat = matrix(1,1,1), 
                           dist_scalar = 1,
                           rho = 1e-2, 
                           tlim = 5,
                           verbose = FALSE)
