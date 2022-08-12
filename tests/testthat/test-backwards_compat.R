testthat::test_that("Confirming deme spatial model is backwards compatible with verison 0.5", {
  swf <- polySimIBD::sim_swf(pos = sort(sample(1:1e3, size = 50)), 
                             N = 20, 
                             m = 0.5, 
                             mean_coi = 1, 
                             migr_dist_mat = 1, 
                             rho = 1e-2, 
                             tlim = 2)
  ARG <- polySimIBD::get_arg(swf = swf, host_index = c(1,2))
  testthat::expect_gt(length(ARG), 0)
})
