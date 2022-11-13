test_that("migration matrix as a diagonal is only within IBD", { 
  diagnl <- matrix(0, 5, 5)
  diag(diagnl) <- 1
  
  diagnl_sim <- polySimIBD::sim_swf(pos =  sort(sample(1.664e6, 1e3)),
                                    migr_mat = diagnl,
                                    N =         rep(5, 5),
                                    m =         rep(0.25, 5),
                                    rho =       7.4e-7,
                                    mean_coi =  rep(2, 5),
                                    tlim =      10)
  
  comb_hosts_df <- t(combn(1:25, 2))
  comb_hosts_list <- split(comb_hosts_df, 1:nrow(comb_hosts_df))
  diagnl_ibd <- purrr::map_dbl(comb_hosts_list, function(hosts, swf) {
    return(polySimIBD::get_pairwise_bv_ibd(swf = swf, host_index = hosts))
  }, swf = diagnl_sim)
  
  
  
  #............................................................
  # tidyout
  #...........................................................
  j1 <- tibble::tibble(smpl1 = 1:25, 
                         deme1 = sort(rep(1:5, 5)))
  j2 <- tibble::tibble(smpl2 = 1:25, 
                         deme2 = sort(rep(1:5, 5)))
  comb_hosts_df <- tibble::as_tibble(comb_hosts_df)
  colnames(comb_hosts_df) <- c("smpl1", "smpl2")
  comb_hosts_df <- comb_hosts_df %>% 
    dplyr::left_join(., j1) %>% 
    dplyr::left_join(., j2)
  
  # diagonal  
  diagnl_comb_hosts_df <- comb_hosts_df %>% 
    dplyr::mutate(ibd = diagnl_ibd) %>% 
    dplyr::filter(deme1 != deme2) %>% 
    dplyr::filter(ibd > 0)
  
  testthat::expect_equal(nrow(diagnl_comb_hosts_df), 0)
  
  
  })




test_that("migration matrix with only one", { 
  obione <- matrix(0, 7, 7)
  obione[3,4] <- 2
  diag(obione) <- 1
  
  obione_sim <- polySimIBD::sim_swf(pos =  sort(sample(1.664e6, 1e3)),
                                    migr_mat = obione,
                                    N =         rep(5, 7),
                                    m =         rep(0.25, 7),
                                    rho =       7.4e-7,
                                    mean_coi =  rep(2, 7),
                                    tlim =      10)
  
  comb_hosts_df34 <- t(combn(1:35, 2))
  comb_hosts_list34 <- split(comb_hosts_df34, 1:nrow(comb_hosts_df34))
  obione_ibd <- purrr::map_dbl(comb_hosts_list34, function(hosts, swf) {
    return(polySimIBD::get_pairwise_bv_ibd(swf = swf, host_index = hosts))
  }, swf = obione_sim)
  
  
  
  #............................................................
  # tidyout
  #...........................................................
  j134 <- tibble::tibble(smpl1 = 1:35, 
                         deme1 = sort(rep(1:7, 5)))
  j234 <- tibble::tibble(smpl2 = 1:35, 
                         deme2 = sort(rep(1:7, 5)))
  comb_hosts_df34 <- tibble::as_tibble(comb_hosts_df34)
  colnames(comb_hosts_df34) <- c("smpl1", "smpl2")
  comb_hosts_df34 <- comb_hosts_df34 %>% 
    dplyr::left_join(., j134) %>% 
    dplyr::left_join(., j234)
  
  # single deme move 
  obione_comb_hosts_df <- comb_hosts_df34 %>% 
    dplyr::mutate(ibd = obione_ibd) %>% 
    dplyr::filter(deme1 != deme2) %>% 
    dplyr::filter(deme1 != 3 & deme2 != 4) %>% 
    dplyr::filter(ibd > 0)
  
  testthat::expect_equal(nrow(obione_comb_hosts_df), 0)
  
  
})



test_that("migration matrix behaves as from-to format", {
  #............................................................
  # make framework
  #...........................................................
  # Lattice
  nCell <- 121
  coords <- round(seq(1, nCell, by = 11))
  latticemodel <- expand.grid(coords, coords)
  colnames(latticemodel) <- c("longnum", "latnum")
  latticemodel$deme <- 1:nrow(latticemodel)
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
  # simulate
  #......................
  spring_swfsim <- polySimIBD::sim_swf(pos =  sort(sample(1.664e6, 1e3)),
                                migr_mat = spring,
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
  
})



test_that("migration matrix in sink only allows for transitivity", { 
  obitwo <- matrix(0, 7, 7)
  obitwo[2,4] <- 2
  obitwo[7,4] <- 2
  diag(obitwo) <- 1
  
  
  obitwo_sim <- polySimIBD::sim_swf(pos =  sort(sample(1.664e6, 1e3)),
                                    migr_mat = obitwo,
                                    N =         rep(5, 7),
                                    m =         rep(0.25, 7),
                                    rho =       7.4e-7,
                                    mean_coi =  rep(2, 7),
                                    tlim =      10)
  
  comb_hosts_df34 <- t(combn(1:35, 2))
  comb_hosts_list34 <- split(comb_hosts_df34, 1:nrow(comb_hosts_df34))
  obitwo_ibd <- purrr::map_dbl(comb_hosts_list34, function(hosts, swf) {
    return(polySimIBD::get_pairwise_bv_ibd(swf = swf, host_index = hosts))
  }, swf = obitwo_sim)
  
  
  
  #............................................................
  # tidyout
  #...........................................................
  j134 <- tibble::tibble(smpl1 = 1:35, 
                         deme1 = sort(rep(1:7, 5)))
  j234 <- tibble::tibble(smpl2 = 1:35, 
                         deme2 = sort(rep(1:7, 5)))
  comb_hosts_df34 <- tibble::as_tibble(comb_hosts_df34)
  colnames(comb_hosts_df34) <- c("smpl1", "smpl2")
  comb_hosts_df34 <- comb_hosts_df34 %>% 
    dplyr::left_join(., j134) %>% 
    dplyr::left_join(., j234)
  
  # single deme move 
  ret <- comb_hosts_df34 %>% 
    dplyr::mutate(ibd = obitwo_ibd) %>% 
    dplyr::filter(deme1 != deme2) %>% 
    dplyr::filter(deme1 != 4) %>% 
    dplyr::filter(deme2 != 4) %>% 
    dplyr::filter(deme1 != 2) %>% # allow for transitivity  
    dplyr::filter(ibd > 0)
  
  
  # sink causes transitivity 
  testthat::expect_equal(nrow(ret), 0)
  
  
  })
