test_that("ARG plot", {
  swfsim <- polySimIBD::sim_swf(pos = sort(sample(1:1e3, size = 100)),
                                N = 20,
                                m = 0.5,
                                rho = 5e-2,
                                mean_coi =  2,
                                tlim = 10)
  argplotlist <- polySimIBD::plotARG(swf = swfsim,  
                                     loci = 1:5, 
                                     host_index = c(1:4), 
                                     tlim = 10)
  
  testthat::expect_gt(length(argplotlist), 0)
})
