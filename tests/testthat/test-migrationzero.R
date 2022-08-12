test_that("migration zero model with 1 deme", {
  swf <- polySimIBD::sim_swf(pos = sort(sample(1:1e3, size = 50)), 
                             N = 25, 
                             m = 0, 
                             mean_coi = 2, 
                             migr_dist_mat = 1, 
                             rho = 1e-2, 
                             tlim = 10)
  
  adj_graph <- t(combn(x = c(1:25), m = 2)) %>% 
    tibble::as_tibble(., .name_repair = "minimal") %>% 
    magrittr::set_colnames(c("smpl1", "smpl2")) %>% 
    dplyr::mutate(pairwiseIBD = purrr::map2_dbl(.x = smpl1, .y = smpl2, .f = function(x,  y){
      polySimIBD:::quiet(
        polySimIBD::get_pairwise_coi_ibd(swf = swf, host_index = c(x,y))
      )
    }))
  # no migration, no co-mingling of hosts, no ibd sharing
  testthat::expect_equal(unique(adj_graph$pairwiseIBD), 0)
})


test_that("migration zero model with multiple demes deme", {
  migmat <- matrix(0, 2, 2)
  diag(migmat) <- 1
  swf <- polySimIBD::sim_swf(pos = sort(sample(1:1e3, size = 50)), 
                             N = c(10, 10), 
                             m = c(0, 0), 
                             mean_coi = c(2, 2),  
                             migr_dist_mat = migmat,
                             rho = 1e-2, 
                             tlim = 10)
  
  adj_graph <- t(combn(x = c(1:20), m = 2)) %>% 
    tibble::as_tibble(., .name_repair = "minimal") %>% 
    magrittr::set_colnames(c("smpl1", "smpl2")) %>% 
    dplyr::mutate(pairwiseIBD = purrr::map2_dbl(.x = smpl1, .y = smpl2, .f = function(x,  y){
      polySimIBD:::quiet(
        polySimIBD::get_pairwise_coi_ibd(swf = swf, host_index = c(x,y))
      )
    }))
  # no migration, no co-mingling of hosts, no ibd sharing
  testthat::expect_equal(unique(adj_graph$pairwiseIBD), 0)
})
