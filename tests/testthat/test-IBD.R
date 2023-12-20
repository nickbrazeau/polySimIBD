testthat::test_that("IBD is zero", {
  N <- c(10, 10)
  m <- c(0.5, 0.5)
  mean_coi <- c(3, 3)
  migr_dist_mat <- matrix(c(1,0,0,1), ncol = 2)
  swf <- polySimIBD::sim_swf(pos = sort(sample(1:1e3, size = 50)), 
                             N = N, 
                             m = m, 
                             mean_coi = mean_coi, 
                             migr_mat = migr_dist_mat, 
                             rho = 1e-2, 
                             tlim = 2)
  
  IBD <- polySimIBD::get_bvibd(swf = swf, host_index = c(1,2))
  testthat::expect_equal(IBD, 0)
})
