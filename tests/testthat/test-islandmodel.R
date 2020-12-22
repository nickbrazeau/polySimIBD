test_that("island model returns two islands", {
  # vectors must be ordered for population A and B
  demesizes <- c(10, 15)
  coimeans <- c(1.25, 1.25)
  m <- c(0.25, 0.25)
  migr_dist_mat <- matrix(0, ncol = 2, nrow = 2)
  diag(migr_dist_mat) <- 100
  swf <- polySimIBD::sim_swf(pos = sample(1:1e3, size = 50), 
                             N = demesizes, 
                             m = m, 
                             mean_coi = coimeans, 
                             migr_dist_mat = migr_dist_mat, 
                             rho = 1e-2, 
                             tlim = 10)
  
  # get adj graph
  adj_graph <- t(combn(x = c(1:sum(demesizes)), m = 2)) %>% 
    tibble::as_tibble(., .name_repair = "minimal") %>% 
    magrittr::set_colnames(c("smpl1", "smpl2")) %>% 
    dplyr::mutate(pairwiseIBD = purrr::map2_dbl(.x = smpl1, .y = smpl2, .f = function(x,  y){
      polySimIBD::quiet(
        polySimIBD::get_realized_pairwise_ibd(swf = swf, host_index = c(x,y))
      )
    }))
  
  # convert to network and do community detection
  nodesdf <- tidygraph::as_tbl_graph(adj_graph, directed = FALSE) %>% 
    dplyr::mutate(community = as.factor(tidygraph::group_louvain())) %>% 
    tidygraph::activate("nodes") %>% 
    tibble::as_tibble()
  
  # should be length two for two isolated populations
  testthat::expect_length(unique(nodesdf$community), 2)
  
})
