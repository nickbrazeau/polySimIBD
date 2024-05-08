testthat::test_that("conversion returns newick string for a simple tree", {
  # make custom bvtree
  bvtree_ex <- new("bvtree", "c" = c(-1,-1,1), "t" = c(-1,-1,1), "z" = c(2,1,0))
 out <- bvtreeToNewick(bvtree = bvtree_ex, tlim = 10) 
 # make sure liftover is correct
 testthat::expect_equal(out,
                        "(Node0:10,(Node1:1,Node2:1):9);")
})


testthat::test_that("conversion returns newick string for multi-overlapping tree", {

  # make custom bvtree
  bvtree_ex <- new("bvtree", "c" = c(-1,-1,0,-1,1), "t" = c(-1,-1,4,-1,6), "z" = c(2,3,0,4,1))
  out <- bvtreeToNewick(bvtree = bvtree_ex, tlim = 10) 
  # make sure liftover is correct
  testthat::expect_equal(out,
                         "((Node0:4,Node2:4):6,(Node1:6,Node4:6):4,Node3:10);")
  
  
})
