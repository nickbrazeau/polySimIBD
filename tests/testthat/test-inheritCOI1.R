test_that("inheritnance of COI monoclonal", {
  # per the original model, we first consider the "deme" (host) and then draw parents
  # prior to allowing new progeny to migrate 
      # Host.cpp line 50 
  # As a result, in instances where we force the COI to be 1, we should always
  # draw from the same "parent" strain within the host. Recombination will 
  # thus occur but it will be undetected from a coalescence POV
  
  ### Magic Numbers 
  rep <- 1:2
  tlim <- 10  
  rho <- 7.4e-7
  pflen <- 1.664e6
  pos <- list(sort(sample(pflen, 1e3)))
  mean_coi <- 1.666424e-07 # force COI=1
  N <- 5
  m <- 0.5
  #............................................................
  # run for basic model
  #...........................................................
  migr_mat <- 1
  snglmt <- tidyr::expand_grid(rep, pos, N, m, mean_coi, rho, tlim, migr_mat)
  snglmt$swfsim <- purrr::pmap(snglmt[2:ncol(snglmt)], polySimIBD::sim_swf)
  
  #............................................................
  # now with migration  matrix
  #...........................................................
  migr_mat <- matrix(data = 0, nrow = 3, ncol = 3)
  diag(migr_mat) <- 1
  # lift over prev sims
  dblmt <- tidyr::expand_grid(rep, pos, N, m, mean_coi, rho, tlim)
  dblmt <- dblmt %>%
    dplyr::mutate(N = purrr::map(N, function(x) rep(x,3)),
                  m = purrr::map(m, function(x) rep(x,3)),
                  mean_coi = purrr::map(mean_coi, function(x) rep(x,3)))
  dblmt <- dblmt %>%
    dplyr::mutate(migr_mat = purrr::map(rep, function(x) return(migr_mat)))
  
  dblmt$swfsim <- purrr::pmap(dblmt[,2:ncol(dblmt)], polySimIBD::sim_swf)
  
  #............................................................
  # now compare/test that parents are always the same in original model 
  # and more complex spatial structure model
  #...........................................................
  testthat::expect_equal(dblmt$swfsim[[1]]$parent_host1[[4]][[1]], 
                         dblmt$swfsim[[1]]$parent_host2[[4]][[1]]) # ind 1 gen 4
  
  testthat::expect_equal(dblmt$swfsim[[2]]$parent_host1[[8]][[2]], 
                         dblmt$swfsim[[2]]$parent_host2[[8]][[2]]) # ind 2 gen 8
  
  testthat::expect_equal(snglmt$swfsim[[1]]$parent_host1[[7]][[3]], 
                         snglmt$swfsim[[1]]$parent_host2[[7]][[3]]) # ind 3 gen 7
  
})
