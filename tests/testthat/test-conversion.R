testthat::test_that("conversion returns newick string", {
  # make custom bvtree
  bvtree_ex <- new("bvtree", "c" = c(-1,-1,1), "t" = c(-1,-1,1), "z" = c(2,0,1))
 out <- bvtreeToNewick(bvtree = bvtree_ex, tlim = 10) 
 # make sure liftover is correct
 testthat::expect_equal(out,
                        "(Node0:10,(Node1:1,Node2:1):9);")
})

