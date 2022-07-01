test_that("accurate IBD for clone***", {
  
  # simulate data 
  swfsim <- polySimIBD::sim_swf(pos = sample(1:1e3, size = 100),
                                N = 20,
                                m = 0.5,
                                rho = 5e-2,
                                mean_coi =  2,
                                tlim = 10)
  ARG <- polySimIBD::get_arg(swfsim, host_index = c(1:2))
  
  # get IBD on individual level is fast
  ind <- polySimIBD::get_realized_pairwise_ibd(swfsim, host_index = c(1,2))
  
  # get IBD for all indiv in a pop is slow
  
