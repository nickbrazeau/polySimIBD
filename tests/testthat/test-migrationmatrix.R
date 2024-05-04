test_that("migration matrix with one on the diagonal is only within-deme IBD", { 
  diagnl <- matrix(0, 5, 5) # 5 demes
  diag(diagnl) <- 1
  
  diagnl_sim <- polySimIBD::sim_swf(pos =  sort(sample(1.664e6, 1e3)),
                                    migr_mat = diagnl,
                                    N =         c(5, 5, 5, 5, 5),
                                    m =         rep(0.25, 5),
                                    rho =       7.4e-7,
                                    mean_coi =  rep(2, 5),
                                    tlim =      10)
  
  comb_hosts_df <- t(combn(1:25, 2))
  comb_hosts_list <- split(comb_hosts_df, 1:nrow(comb_hosts_df))
  diagnl_ibd <- purrr::map_dbl(comb_hosts_list, function(hosts, swf) {
    return(polySimIBD::get_bvibd(swf = swf, host_index = hosts))
  }, swf = diagnl_sim)
  
  #............................................................
  # tidyout
  #...........................................................
  j1 <- tibble::tibble(smpl1 = 1:25, 
                       deme1 = sort(rep(1:5, 5)))
  j2 <- tibble::tibble(smpl2 = 1:25, 
                       deme2 = sort(rep(1:5, 5)))
  comb_hosts_df <- tibble::as_tibble(comb_hosts_df, .name_repair = "minimal")
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



test_that("migration matrix with zero on diagonal is only between-deme IBD", {
  # NB in large population, shouldn't draw same parent 
  no_diagnl <- matrix(1, 5, 5)
  diag(no_diagnl) <- 0
  
  no_diagnl_sim <- polySimIBD::sim_swf(pos =  sort(sample(1.664e6, 1e3)),
                                       migr_mat = no_diagnl,
                                       N =         rep(1e3, 5),
                                       m =         rep(0.25, 5),
                                       rho =       7.4e-7,
                                       mean_coi =  rep(2, 5),
                                       tlim =      10)
  
  comb_hosts_df <- t(combn(1:1e3, 2))
  # down sample for memory
  rws <- sample(1:nrow(comb_hosts_df), size = 50)
  comb_hosts_df <- comb_hosts_df[rws,]
  comb_hosts_list <- split(comb_hosts_df, 1:nrow(comb_hosts_df))
  no_diagnl_ibd <- purrr::map_dbl(comb_hosts_list, function(hosts, swf) {
    return(polySimIBD::get_bvibd(swf = swf, host_index = hosts))
  }, swf = no_diagnl_sim)
  
  #............................................................
  # tidyout
  #...........................................................
  j1 <- tibble::tibble(smpl1 = 1:1e3, 
                       deme1 = sort(rep(1:5, 1e3/5)))
  j2 <- tibble::tibble(smpl2 = 1:1e3, 
                       deme2 = sort(rep(1:5, 1e3/5)))
  comb_hosts_df <- tibble::as_tibble(comb_hosts_df, .name_repair = "minimal")
  colnames(comb_hosts_df) <- c("smpl1", "smpl2")
  comb_hosts_df <- comb_hosts_df %>% 
    dplyr::left_join(., j1) %>% 
    dplyr::left_join(., j2)
  
  # no diagonal  
  no_diagnl_comb_hosts_df <- comb_hosts_df %>% 
    dplyr::mutate(ibd = no_diagnl_ibd) %>% 
    dplyr::filter(deme1 == deme2) %>% # looking for within-deme when only between is allowed
    dplyr::filter(ibd > 0)
  
  testthat::expect_lt(nrow(no_diagnl_comb_hosts_df), 2) # will add a tolerance of 2 for very unlucky case of drawing same parent in large pop
  
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
    return(polySimIBD::get_bvibd(swf = swf, host_index = hosts))
  }, swf = obione_sim)
  
  #............................................................
  # tidyout
  #...........................................................
  j134 <- tibble::tibble(smpl1 = 1:35, 
                         deme1 = sort(rep(1:7, 5)))
  j234 <- tibble::tibble(smpl2 = 1:35, 
                         deme2 = sort(rep(1:7, 5)))
  comb_hosts_df34 <- tibble::as_tibble(comb_hosts_df34, .name_repair = "minimal")
  colnames(comb_hosts_df34) <- c("smpl1", "smpl2")
  comb_hosts_df34 <- comb_hosts_df34 %>% 
    dplyr::left_join(., j134) %>% 
    dplyr::left_join(., j234)
  
  # single deme move 
  obione_comb_hosts_df_dontcatch <- comb_hosts_df34 %>% 
    dplyr::mutate(ibd = obione_ibd) %>% 
    dplyr::filter(deme1 != deme2) %>% 
    dplyr::filter(deme1 != 3 & deme2 != 4) %>% 
    dplyr::filter(ibd > 0) # look for where we exclude one instance of between, so should be no sharing 
  
  testthat::expect_equal(nrow(obione_comb_hosts_df_dontcatch), 0)
  
  # now long for where we INCLUDE only one instance of between 
  obione_comb_hosts_df_dontcatch <- comb_hosts_df34 %>% 
    dplyr::mutate(ibd = obione_ibd) %>% 
    dplyr::filter(deme1 != deme2) %>% 
    dplyr::filter(deme1 == 3 & deme2 == 4) %>% 
    dplyr::filter(ibd > 0)  
  
  testthat::expect_gt(nrow(obione_comb_hosts_df_dontcatch), 0)
  
})



test_that("migration matrix behaves as to-from format", {
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
  diag(spring) <- 1 # if not out, stay in home deme
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
    return(polySimIBD::get_bvibd(swf = swf, host_index = hosts))
  }, swf = spring_swfsim)
  #............................................................
  # tidyout
  #...........................................................
  j1 <- tibble::tibble(smpl1 = 1:35, 
                       deme1 = sort(rep(1:7, 5)))
  j2 <- tibble::tibble(smpl2 = 1:35, 
                       deme2 = sort(rep(1:7, 5)))
  comb_hosts_df <- tibble::as_tibble(comb_hosts_df, .name_repair = "minimal")
  colnames(comb_hosts_df) <- c("smpl1", "smpl2")
  comb_hosts_df <- comb_hosts_df %>% 
    dplyr::left_join(., j1) %>% 
    dplyr::left_join(., j2)
  
  # spring 
  spring_comb_hosts_df <- comb_hosts_df %>% 
    dplyr::mutate(ibd = spring_ibd) %>% 
    dplyr::filter(deme1 != deme2)
  
  # only deme 4 can have positive genetic connections
  # include 4 
  spring_comb_hosts_df_include4 <- spring_comb_hosts_df %>% 
    dplyr::filter(deme1 == 4 | deme2 == 4) %>% 
    dplyr::filter(ibd > 0)
  testthat::expect_gt(nrow(spring_comb_hosts_df_include4), 0)
  # exclude 4 
  spring_comb_hosts_df_exclude4 <- spring_comb_hosts_df %>% 
    dplyr::filter(deme1 != 4 & deme2 != 4) %>% 
    dplyr::filter(ibd > 0)
  testthat::expect_equal(nrow(spring_comb_hosts_df_exclude4), 0)
  
})


