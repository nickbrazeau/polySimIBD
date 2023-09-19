testthat::test_that("Plotting functions work", {
  swf <- polySimIBD::sim_swf(pos = sort(sample(1:1e3, size = 50)), 
                             N = 20, 
                             m = 0.5, 
                             mean_coi = 1, 
                             migr_mat = 1, 
                             rho = 1e-2, 
                             tlim = 2)
  ARG <- polySimIBD::get_arg(swf = swf, host_index = c(1,2))
  # plot coal trees
  plotObj <- plot_coalescence_trees(ARG = ARG, loci = 1)
  testthat::expect_true( any(class(plotObj) %in% "ggplot" ))
    
  })

