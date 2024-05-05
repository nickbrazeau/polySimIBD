test_that("test Cpp IBD calculation by hand", {
  # NB, cpp IBD calculation already checked in migration matrix form multiple times
  # here handwriting ARG to confirm IBD
  swf <- polySimIBD::sim_swf(pos = c(0,1e3,2e3), 
                             N = 2, 
                             m = 0.5, 
                             mean_coi = 1, 
                             migr_mat = 1, 
                             rho = 1e-2, 
                             tlim = 2)
  # call bvibd
  bvIBD <- polySimIBD::get_bvibd(swf = swf, host_index = c(1,2))
  # calculate by hand from ARG
  ARG <- polySimIBD::get_arg(swf = swf, host_index = c(1,2))
  coi <- swf$coi
  bvtree1 <- ARG[[1]]@c
  bvtree2 <- ARG[[1]]@c 
  bvtree3 <- ARG[[1]]@c 
  
  # calculate bvtree1 - always look left 
  h1 <- sum( bvtree1[coi[2]] %in% (coi[1]:coi[2] - 1) ) 
  h2 <- sum( bvtree2[coi[2]] %in% (coi[1]:coi[2] - 1) ) 
  h3 <- sum( bvtree2[coi[2]] %in% (coi[1]:coi[2] - 1) ) 
  handIBD <- sum(h1,h2,h3)/sum(coi)
  
  # confirm
  testthat::expect_equal(bvIBD, handIBD)
  
})
