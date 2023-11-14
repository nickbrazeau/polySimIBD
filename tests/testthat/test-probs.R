testthat::test_that("probability functions work", {
  num3 <- rdirichlet(shape = rep(1,3)) 
  num1 <- rbetabinom(n = 1, k = 10, alpha = 1, beta = 1)
  # tests out
  testthat::expect_length(num3, 3)
  testthat::expect_length(num1, 1)
})

