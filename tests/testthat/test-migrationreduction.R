test_that("Migration with only within host migration reduces back to original model", {
  ### Magic Numbers: immutable throughout simulations
  rep <- 1:100 # number of simulation realizations to perform
  tlim <- 10 # assume IBD to 10 generations for recent coalescent
  # Miles et al. 2016 (PMC5052046) & Taylor et al. 2019 (PMC6707449) gives us a recombination rate by 7.4e-7 M/bp
  # Aimee gets this number by taking the inverse of Mile's estimate of the CO recombination rate of 13.5 kb/cM
  rho <- 7.4e-7
  # approximate average of Pf3d7 Chromosome Lengths
  pflen <- 1.664e6
  # assuming single chromosome for ease
  # assuming 1e3 loci
  pos <- list(sort(sample(pflen, 1e3)))
  mean_coi <- 1.666424e-07 # from optim lamnbda
  # small N for convenience
  N <- 10
  m <- seq(0.05, 0.95, by = 0.3)
  #............................................................
  # run for basic model
  #...........................................................
  migr_mat <- 1
  snglmt <- tidyr::expand_grid(rep, pos, N, m, mean_coi, rho, tlim, migr_mat)
  snglmt$swfsim <- purrr::pmap(snglmt[2:ncol(snglmt)], polySimIBD::sim_swf)
  
  get_ibd <- function(N, swf) {
    comb_hosts_df <- t(combn(N, 2))
    comb_hosts_list <- split(comb_hosts_df, 1:nrow(comb_hosts_df))
    # get pairwise IBD
    ibd <- purrr::map_dbl(comb_hosts_list, function(hosts, swf) {
      return(polySimIBD::get_bvibd(swf = swf, host_index = hosts))
    }, swf = swf)
    return(ibd)
  }
  # get ibd
  snglmt$ibdcalc <- purrr::map(snglmt$swfsim, get_ibd, N = N)
  # summarize
  snglmt_ret <- snglmt %>%
    dplyr::mutate(meanIBD = purrr::map_dbl(ibdcalc, mean)) %>%
    dplyr::select(c("rep", "meanIBD"))
  
  
  
  #............................................................
  # now with migration  matrix
  #...........................................................
  migr_mat <- matrix(data = 0, nrow = 2, ncol = 2)
  diag(migr_mat) <- 1
  # lift over prev sims
  dblmt <- tidyr::expand_grid(rep, pos, N, m, mean_coi, rho, tlim)
  dblmt <- dblmt %>%
    dplyr::mutate(N = purrr::map(N, function(x) rep(x,2)),
                  m = purrr::map(m, function(x) rep(x,2)),
                  mean_coi = purrr::map(mean_coi, function(x) rep(x,2)))
  dblmt <- dblmt %>%
    dplyr::mutate(migr_mat = purrr::map(rep, function(x) return(migr_mat)))
  
  dblmt$swfsim <- purrr::pmap(dblmt[2:ncol(dblmt)], polySimIBD::sim_swf)
  
  
  wrapibdcalc <- function(swfsim, migr_mat){
    # expand hosts
    comb_hosts_df <- t(combn(N*nrow(migr_mat), 2))
    comb_hosts_list <- split(comb_hosts_df, 1:nrow(comb_hosts_df))
    # get pairwise IBD
    ibdmig <- purrr::map_dbl(comb_hosts_list, function(hosts, swf) {
      return(polySimIBD::get_bvibd(swf = swf, host_index = hosts))
    }, swf = swfsim)
    # tidyout
    j1 <- tibble::tibble(smpl1 = 1:(N*nrow(migr_mat)),
                         deme1 = sort(rep(1:nrow(migr_mat), N)))
    j2 <- tibble::tibble(smpl2 = 1:(N*nrow(migr_mat)),
                         deme2 = sort(rep(1:nrow(migr_mat), N)))
    comb_hosts_df <- tibble::as_tibble(comb_hosts_df, .name_repair = "minimal")
    colnames(comb_hosts_df) <- c("smpl1", "smpl2")
    comb_hosts_df <- comb_hosts_df %>%
      dplyr::left_join(., j1) %>%
      dplyr::left_join(., j2) %>%
      dplyr::mutate(ibd = ibdmig)
    return(comb_hosts_df)
  }
  # get ibd
  dblmt$ibdcalc <- purrr::pmap(dblmt[,c("swfsim", "migr_mat")], wrapibdcalc)
  
  getmeanmig <- function(ibdcalc){
    ibdcalc %>%
      dplyr::filter(deme1 == deme2) %>% 
      dplyr::summarise(meanibd = mean(ibd)) %>% 
      dplyr::pull(meanibd)
  }
  # summarize
  dblmt_ret <- dblmt %>%
    dplyr::mutate(meanIBD = purrr::map_dbl(ibdcalc, getmeanmig)) %>%
    dplyr::select(c("rep", "meanIBD"))
  
  #............................................................
  # now compare
  #...........................................................
  snglmt_ret <- snglmt_ret %>%
    dplyr::mutate(state = "orig")
  
  dblmt_ret <- dblmt_ret %>%
    dplyr::mutate(state = "migmat")
  
  ret_comb <- dplyr::bind_rows(snglmt_ret, dblmt_ret)
  
  # ggplot() + 
  #   geom_histogram(data = ret_comb, aes(x = meanIBD, fill = state), alpha = 0.5) + 
  #   scale_fill_viridis_d()
  
  # KW test to ensure that it is from the same distribution 
  checksame <- kruskal.test(ret_comb$meanIBD, g = factor(ret_comb$state))
  testthat::expect_gt(checksame$p.value, 0.05)
})
