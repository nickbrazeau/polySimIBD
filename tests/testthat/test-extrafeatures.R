testthat::test_that("hapmat extraction returns matrix and it is unique haplotypes in the island model", {
  # ISLAND MODEL
  demesizes <- c(10, 15)
  coimeans <- c(1, 1)
  m <- c(0.25, 0.25)
  migr_dist_mat <- matrix(0, ncol = 2, nrow = 2)
  diag(migr_dist_mat) <- 100
  swf <- polySimIBD::sim_swf(pos = sort(sample(1:1e3, size = 50)), 
                             N = demesizes, 
                             m = m, 
                             mean_coi = coimeans, 
                             migr_mat = migr_dist_mat, 
                             rho = 1e-2, 
                             tlim = 10)
  
  ARG <- polySimIBD::get_arg(swf = swf, host_index = c(1,11))
  # because island model, hapmat should have independent colums
  hapmat <- polySimIBD::extract_haplotype_matrix(ARG)
  testthat::expect_is(hapmat, "matrix")
  # first column will be 1 (because 1 haplotype)
  # first column of next host will be 2 (because 2 haplotype)
  testthat::expect_equal(
    1,
    unique(hapmat[,sum(swf$coi[c(1,11)])]) -  unique(hapmat[,1]) 
  )
  
})



testthat::test_that("hapmat extraction returns matrix and it is unique haplotypes in the island model", {
  # ISLAND MODEL
  demesizes <- c(10, 15)
  coimeans <- c(1, 1)
  m <- c(0.25, 0.25)
  migr_dist_mat <- matrix(0, ncol = 2, nrow = 2)
  diag(migr_dist_mat) <- 100
  swf <- polySimIBD::sim_swf(pos = sort(sample(1:1e3, size = 50)), 
                             N = demesizes, 
                             m = m, 
                             mean_coi = coimeans, 
                             migr_mat = migr_dist_mat, 
                             rho = 1e-2, 
                             tlim = 10)
  
  ARG <- polySimIBD::get_arg(swf = swf, host_index = c(1,11))
  hapmat <- polySimIBD::extract_haplotype_matrix(ARG)
  bisimdata <- sim_biallelic(COIs = swf$coi[c(1,11)],
                             haplotypematrix = hapmat) # default parameters otherwise
  # expect a list of the data
  testthat::expect_is(bisimdata, "list")
  
})
