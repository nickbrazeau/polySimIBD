test_that("migration matrix behaves as from-to format", {
  #............................................................
  # make framework
  #...........................................................
  # Lattice
  nCell <- 121
  coords <- round(seq(1, nCell, by = 11))
  latticemodel <- expand.grid(coords, coords)
  colnames(latticemodel) <- c("longnum", "latnum")
  latticemodel$deme <- 1:nrow(latticemodel))
  latticemodel <- latticemodel[c(1,11,58,61,64,111,121),]
  ogdistmat <- matrix(0, nrow = 7, ncol = 7)
  rownames(ogdistmat) <- colnames(ogdistmat) <- latticemodel$deme
  
  #......................
  # spring
  #......................
  spring <- ogdistmat
  spring[4,] <- 2 # middle deme has connections OUT only
  diag(spring) <- 1 # if not out, stat in home deme
  # liftover to rates
  spring <- spring/sum(spring)
  
  #......................
  # blackhole
  #......................
  blackhole <- ogdistmat
  blackhole[,4] <- 2 # middle deme has connections IN only
  diag(blackhole) <- 1 # if not out, stat in home deme
  # liftover to rates
  blackhole <- blackhole/sum(blackhole)
  
  #......................
  # simulate
  #......................
  spring_swfsim <- polySimIBD::sim_swf(pos =  sort(sample(1.664e6, 1e3)),
                                migr_mat = spring,
                                N =         rep(5, 7),
                                m =         rep(0.25, 7),
                                rho =       7.4e-7,
                                mean_coi =  rep(2, 7),
                                tlim =      10)
  
  blackhole_swfsim <- polySimIBD::sim_swf(pos =  sort(sample(1.664e6, 1e3)),
                                          migr_mat = blackhole,
                                          N =         rep(5, 7),
                                          m =         rep(0.25, 7),
                                          rho =       7.4e-7,
                                          mean_coi =  rep(2, 7),
                                          tlim =      10)
  
  
  #............................................................
  # get ibd
  #...........................................................
  comb_hosts_df <- t(combn(1:35, 2))
  comb_hosts_list <- split(comb_hosts_df, 1:nrow(comb_hosts_df))
  # spring
  spring_ibd <- purrr::map_dbl(comb_hosts_list, function(hosts, swf) {
    return(polySimIBD::get_pairwise_bv_ibd(swf = swf, host_index = hosts))
  }, swf = spring_swfsim)
  
  # blackhole
  blackhole_ibd <- purrr::map_dbl(comb_hosts_list, function(hosts, swf) {
    return(polySimIBD::get_pairwise_bv_ibd(swf = swf, host_index = hosts))
  }, swf = blackhole_swfsim)
  
  #............................................................
  # tidyout
  #...........................................................
  j1 <- tibble::tibble(smpl1 = 1:35, 
                       deme1 = sort(rep(1:7, 5)))
  j2 <- tibble::tibble(smpl2 = 1:35, 
                       deme2 = sort(rep(1:7, 5)))
  comb_hosts_df <- tibble::as_tibble(comb_hosts_df)
  colnames(comb_hosts_df) <- c("smpl1", "smpl2")
  comb_hosts_df <- comb_hosts_df %>% 
    dplyr::left_join(., j1) %>% 
    dplyr::left_join(., j2)
  
  # spring 
  spring_comb_hosts_df <- comb_hosts_df %>% 
    dplyr::mutate(ibd = spring_ibd) %>% 
    dplyr::filter(deme1 != deme2)
  
  # only deme 4 can have positive genetic connections
  spring_comb_hosts_df <- spring_comb_hosts_df %>% 
    dplyr::filter(deme1 != 4) %>% 
    dplyr::filter(deme2 != 4) %>% 
    dplyr::filter(ibd > 0)
  
  testthat::expect_equal(nrow(spring_comb_hosts_df), 0)
  
  
  # blackhole 
  blackhole_comb_hosts_df <- comb_hosts_df %>% 
    dplyr::mutate(ibd = blackhole_ibd) %>% 
    dplyr::filter(deme1 != deme2)
  # only deme 4 can have positive genetic connections
  blackhole_comb_hosts_df <- blackhole_comb_hosts_df %>% 
    dplyr::filter(deme1 != 4) %>% 
    dplyr::filter(deme2 != 4) %>% 
    dplyr::filter(ibd > 0)
  
  testthat::expect_equal(nrow(blackhole_comb_hosts_df), 0)
  
  
  
})
