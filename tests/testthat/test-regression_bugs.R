# Unit Tests for Previously Identified Bugs (Regression Tests)
# Author: Claude Code (Code Review AI)
# Human reviewed for accuracy 
# Date: 2025-11-12
# Purpose: Ensure previously identified bugs don't reappear


# ============================================================================
# Test for Underflow Handling in Recombination Probability
# ============================================================================

testthat::test_that("odd_prob calculation handles extreme values", {
  # Test with very large distance that might cause underflow
  pos_extreme <- c(0, 1e10)  # Very large distance

  swf <- polySimIBD::sim_swf(
    pos = pos_extreme,
    N = 3,
    m = 0.5,
    mean_coi = 2,
    migr_mat = 1,
    rho = 1e-5,
    tlim = 2
  )

  testthat::expect_s3_class(swf, "swfsim")
})


testthat::test_that("odd_prob calculation handles very small rho", {
  # Test with very small recombination rate
  swf <- polySimIBD::sim_swf(
    pos = c(0, 1000, 2000),
    N = 3,
    m = 0.5,
    mean_coi = 2,
    migr_mat = 1,
    rho = 1e-15,  # Extremely small
    tlim = 2
  )

  testthat::expect_s3_class(swf, "swfsim")
})


# ============================================================================
# Test for COI = 1 Edge Case in get_within_ibd
# ============================================================================

testthat::test_that("get_within_ibd properly rejects COI = 1", {
  # Create a scenario where we can control COI
  set.seed(42)

  swf <- polySimIBD::sim_swf(
    pos = c(0, 1000),
    N = 10,
    m = 0.5,
    mean_coi = 0.5,  # Low mean_coi to get some COI = 1
    migr_mat = 1,
    rho = 1e-2,
    tlim = 2
  )

  # Find a host with COI = 1
  host_coi1 <- which(swf$coi == 1)[1]

  if (is.na(host_coi1)) {
    testthat::skip("No host with COI = 1 in this sample")
  }

  # Should error with appropriate message
  testthat::expect_error(
    polySimIBD::get_within_ibd(swf = swf, host_index = host_coi1),
    "COI is 1"
  )
})


# ============================================================================
# Test for Consistency Between ARG and IBD Calculations
# ============================================================================

testthat::test_that("ARG-based IBD matches direct IBD calculation", {
  # This is an integration test to ensure consistency

  set.seed(5555)
  swf <- polySimIBD::sim_swf(
    pos = c(0, 1000, 2000, 3000),
    N = c(3, 3),
    m = c(0.5, 0.5),
    mean_coi = c(2, 2),
    migr_mat = matrix(c(0.9, 0.1, 0.1, 0.9), nrow = 2),
    rho = 1e-2,
    tlim = 3
  )

  # Calculate IBD using get_bvibd
  weights <- c(min(swf$pos),  diff(swf$pos))
  weights <- weights / sum(weights)
  bvibd_direct <- polySimIBD::get_bvibd(
    swf = swf,
    host_index = c(1, 4),
    weight_loci = weights
  )

  # Calculate using ARG and manual computation
  arg <- polySimIBD::get_arg(swf = swf, host_index = c(1, 4))
  coi <- swf$coi[c(1, 4)]

  # Manual calculation
  ibd_manual <- numeric(length(arg))
  for (i in 1:length(arg)) {  # Check all loci
    bvtree <- arg[[i]]@c
    # Check if any haplotype from host 2 coalesces with host 1
    ibd_manual[i] <- any(bvtree[(coi[1] + 1):sum(coi)] %in% (1:coi[1] - 1))
  }

  bvibd_manual <- sum(ibd_manual * weights) / sum(weights)

  # Should be very close (allowing for numerical precision)
  testthat::expect_equal(bvibd_direct, bvibd_manual, tolerance = 1e-10)
})


# ============================================================================
# Test for Memory/Performance with Large Simulations
# ============================================================================

testthat::test_that("large simulation completes without error", {
  # This is more of a stress test

  swf <- polySimIBD::sim_swf(
    pos = seq(0, 5000, by = 100),  # 51 loci
    N = 20,
    m = 0.5,
    mean_coi = 5,
    migr_mat = 1,
    rho = 1e-3,
    tlim = 10
  )

  testthat::expect_s3_class(swf, "swfsim")
  testthat::expect_equal(length(swf$coi), 20)

  # Get ARG for subset to test memory handling
  arg <- polySimIBD::get_arg(swf = swf, host_index = 1:5)
  testthat::expect_s3_class(arg, "argraph")
})


# ============================================================================
# Test for Deterministic Behavior with Set Seed
# ============================================================================

testthat::test_that("simulations are reproducible with set seed", {
  # Set seed and run simulation
  set.seed(12345)
  swf1 <- polySimIBD::sim_swf(
    pos = c(0, 1000, 2000),
    N = 5,
    m = 0.5,
    mean_coi = 3,
    migr_mat = 1,
    rho = 1e-2,
    tlim = 3
  )

  # Reset seed and run again
  set.seed(12345)
  swf2 <- polySimIBD::sim_swf(
    pos = c(0, 1000, 2000),
    N = 5,
    m = 0.5,
    mean_coi = 3,
    migr_mat = 1,
    rho = 1e-2,
    tlim = 3
  )

  # Should produce identical results
  testthat::expect_identical(swf1$coi, swf2$coi)
  testthat::expect_identical(swf1$recomb, swf2$recomb)
})


# ============================================================================
# Test for Proper Error Messages
# ============================================================================

testthat::test_that("error messages are informative", {
  # Test that error messages guide users properly

  # Non-increasing positions
  testthat::expect_error(
    polySimIBD::sim_swf(
      pos = c(0, 2000, 1000),
      N = 5,
      m = 0.5,
      mean_coi = 2,
      rho = 1e-2,
      tlim = 2
    ),
    "increasing"
  )

  # Mismatched dimensions
  testthat::expect_error(
    polySimIBD::sim_swf(
      pos = c(0, 1000),
      N = c(5, 5),
      m = c(0.5, 0.5),
      mean_coi = c(2, 2),
      migr_mat = matrix(c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
                        nrow = 3),  # 3x3 matrix but only 2 demes
      rho = 1e-2,
      tlim = 2
    ),
    "dimensions"
  )
})
