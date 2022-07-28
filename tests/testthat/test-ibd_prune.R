test_that("accurate IBD prune when within coalescent does not contribute to between coalescence", { 
  #......................
  # test left side 
  #......................
  # host dist
  coi <- c(2, 3)
  # conn and time
  t <- c(-1,5, 8,3-1,5)
  c <- c(-1,0, 0,-1,3)
  # update z
  tmp <- t
  tmp[tmp == -1] <- NA
  z <- order(tmp) - 1
  # make tree
  win_not_part_of_btwn1 <- new("bvtree", c = c, t  = t, z = z)
  out <- polySimIBD:::get_withinIBD_bvtree_subset(c = win_not_part_of_btwn1@c,
                                                  coi1 = 2,
                                                  coi2 = 3)
  testthat::expect_equal(out, c(1,2,3))
  #......................
  # test  right side 
  #......................
  # host dist
  coi <- c(4, 1)
  # conn and time
  t <- c(-1,5,-1,3, 7)
  c <- c(-1,0,-1,2, 2)
  # update z
  tmp <- t
  tmp[tmp == -1] <- NA
  z <- order(tmp) - 1
  # make tree
  win_not_part_of_btwn2 <- new("bvtree", c = c, t = t, z = z)
  out <- polySimIBD:::get_withinIBD_bvtree_subset(c = win_not_part_of_btwn2@c,
                                                  coi1 = 4,
                                                  coi2 = 1)
  testthat::expect_equal(out, c(3,4,5))
  
  #......................
  # test nested right side
  #......................
  # host dist
  coi <- c(4, 1)
  # conn and time
  t <- c(-1,8,7,6, 5)
  c <- c(-1,0,1,2, 3)
  # update z
  tmp <- t
  tmp[tmp == -1] <- NA
  z <- order(tmp) - 1
  # make tree
  win_not_part_of_btwn3 <- new("bvtree", c = c, t = t, z = z)
  out <- polySimIBD:::get_withinIBD_bvtree_subset(c = win_not_part_of_btwn3@c,
                                                  coi1 = 4,
                                                  coi2 = 1)
  testthat::expect_equal(out, c(1:5))
  
  })
