test_that("effective COI returns vector", {
  N <- c(3, 2, 7)
  m <- c(0.5, 0.5, 0.5)
  mean_coi <- c(1, 1, 1)
  migr_dist_mat <- matrix(runif(9), 3, 3)
  swf <- sim_swf(pos = sort(sample(1:1e3, size = 50)), 
                             N = N, 
                             m = m, 
                             mean_coi = mean_coi, 
                             migr_mat = migr_dist_mat, 
                             rho = 1e-2, 
                             tlim = 2)
  
  coi <- get_effective_coi(swf = swf, host_index = 1)
  testthat::expect_vector(coi)
})


test_that("withinIBD returns double", {
  N <- c(3, 2, 7)
  m <- c(0.5, 0.5, 0.5)
  mean_coi <- c(20, 20, 20)
  migr_dist_mat <- matrix(runif(9), 3, 3)
  swf <- sim_swf(pos = sort(sample(1:1e3, size = 50)), 
                 N = N, 
                 m = m, 
                 mean_coi = mean_coi, 
                 migr_mat = migr_dist_mat, 
                 rho = 1e-2, 
                 tlim = 2)
  
  wibd <- get_within_ibd(swf = swf, host_index = 1)
  testthat::expect_is(wibd, "numeric")
})
