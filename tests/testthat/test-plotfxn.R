test_that("more color for plot", {
  # simple vector of colors from Rcolorbrewer
  vec <- polySimIBD:::more_colors()
  testthat::expect_vector(vec)
})


test_that("check plot of barcodes from hapmat", {
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
  barcodeplot <- plot_barcodes(hapmat = hapmat,
                               coi = swf$coi[c(1,11)])
  
  # expect that it is a plot
  testthat::expect_is(barcodeplot, "gg")
  
})
