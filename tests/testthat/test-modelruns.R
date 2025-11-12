testthat::test_that("model runs", {
  N <- c(3, 2, 7)
  m <- c(0.5, 0.5, 0.5)
  mean_coi <- c(20, 20, 20)
  migr_dist_mat <- matrix(runif(9), 3, 3)
  swf <- polySimIBD::sim_swf(pos = sort(sample(1:1e3, size = 50)), 
                             N = N, 
                             m = m, 
                             mean_coi = mean_coi, 
                             migr_mat = migr_dist_mat, 
                             rho = 1e-2, 
                             tlim = 2)
  # swf test
  testthat::expect_is(swf, "swfsim")
  # arg test
  ARG <- polySimIBD::get_arg(swf = swf, host_index = c(1,2))
  testthat::expect_is(ARG, "argraph")
  testthat::expect_gt(length(ARG), 0)
  # bv tree sub
  bvsubtree <- polySimIBD::subset_bvtree(bvtree = ARG[[1]], s = c(1,2))
  testthat::expect_is(bvsubtree, "bvtree")
})

