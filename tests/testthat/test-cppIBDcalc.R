test_that("test Cpp IBD calculation by hand", {
  # NB, cpp IBD calculation already checked in migration matrix form multiple times
  # here handwriting ARG to confirm IBD
  swf <- polySimIBD::sim_swf(pos = c(0,1001,2001,3001), 
                             N = 2, 
                             m = 0.5, 
                             mean_coi = 1, 
                             migr_mat = 1, 
                             rho = 1e-2, 
                             tlim = 25)
  # call bvibd
  wi <- c(min(swf$pos), diff(swf$pos))
  wi <- wi/sum(wi)
  bvIBD <- polySimIBD::get_bvibd(swf = swf, host_index = c(1,2), weight_loci = wi)
  
  
  # calculate by hand from ARG
  ARG <- polySimIBD::get_arg(swf = swf, host_index = c(1,2))
  (coi <- swf$coi)
  # always burn first loci
  bvtree1 <- ARG[[2]]@c
  bvtree2 <- ARG[[3]]@c 
  bvtree3 <- ARG[[4]]@c 
  
  # calculate bvtree1 - always look left 
  h1 <- any( bvtree1[(coi[1]+1):sum(coi)] %in% (1:coi[1] - 1) ) 
  h2 <- any( bvtree2[(coi[1]+1):sum(coi)] %in% (1:coi[1] - 1) ) 
  h3 <- any( bvtree3[(coi[1]+1):sum(coi)] %in% (1:coi[1] - 1) ) 
  # under SNP vs PSMC (Li/Durbin model)  
  wi <- c(min(swf$pos), diff(swf$pos))
  wi <- wi/sum(wi)
  # weighted average (each loci, denom is 1)
  handIBD <- sum(c(h1,h2,h3) * wi[2:4]) # 3 loci, equidistant and weighted thus
  
  
  # confirm
  testthat::expect_equal(bvIBD, handIBD)
  
})
