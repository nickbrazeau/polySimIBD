test_that("accurate IBD for clone***", {
  
  # simulate data 
  swfsim <- polySimIBD::sim_swf(pos = sort(sample(1:1e3, size = 100)),
                                N = 20,
                                m = 0.5,
                                rho = 5e-2,
                                mean_coi =  1,
                                migr_dist_mat = 1,
                                tlim = 20)
  
  n <- Sys.time()
  polySimIBD::get_realized_pairwise_ibd(swf = swfsim, host_index = c(1,2))
  tdiff <- Sys.time() - n
  (tdiff * choose(2500,2)) / 3600
  
  
# arg for exploration
ARG <- polySimIBD::get_arg(swfsim, host_index = c(1:2))

#............................................................
# playground ind loci
#...........................................................
uniqueloci <- unique(ARG)

unlist(lapply(ARG, function(x){identical(x, uniqueloci[[1]])}))
plot(unlist(lapply(ARG, function(x){identical(x, uniqueloci[[1]])})), type = "l")
sum(unlist(lapply(ARG, function(x){identical(x, uniqueloci[[1]])})))

#............................................................
# IBD
#...........................................................

# get IBD on individual level is fast
ind <- polySimIBD::get_realized_pairwise_ibd(swfsim, host_index = c(1,2))

# get IBD for all indiv in a pop is slow
n <- Sys.time()
polySimIBD::get_realized_pairwise_ibd(swfsim, host_index = c(1,2))
tdiff <- Sys.time() - n
(tdiff * choose(2500,2)) / 3600




} 
)
