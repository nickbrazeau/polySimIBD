testthat::test_that("model runs", {
  N <- c(3, 2, 7)
  m <- c(0.5, 0.5, 0.5)
  mean_coi <- c(1, 1, 1)
  migr_dist_mat <- matrix(runif(9), 3, 3)
  swf <- polySimIBD::sim_swf(pos = sample(1:1e3, size = 50), 
                             N = N, 
                             m = m, 
                             mean_coi = mean_coi, 
                             migr_dist_mat = migr_dist_mat, 
                             rho = 1e-2, 
                             tlim = 2)
  
  ARG <- polySimIBD::get_arg(swf = swf, host_index = c(1,2))
  testthat::expect_gt(length(ARG), 0)
})
