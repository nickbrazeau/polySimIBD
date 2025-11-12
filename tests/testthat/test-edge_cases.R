# Unit Tests for Edge Cases and Input Validation
# Author: Claude Code (Code Review AI)
# Human reviewed for accuracy
# Date: 2025-11-12
# Purpose: Test edge cases, error handling, and boundary conditions


# ============================================================================
# Tests for Migration Matrix Edge Cases
# ============================================================================

testthat::test_that("migration matrix with extreme probabilities works", {
  # Test with very high migration probability
  migr_mat_high <- matrix(c(0.01, 0.99,
                            0.99, 0.01),
                          nrow = 2, byrow = TRUE)

  swf <- polySimIBD::sim_swf(
    pos = c(0, 1000, 2000),
    N = c(5, 5),
    m = c(0.5, 0.5),
    mean_coi = c(2, 2),
    migr_mat = migr_mat_high,
    rho = 1e-3,
    tlim = 3
  )

  testthat::expect_s3_class(swf, "swfsim")
  testthat::expect_equal(length(swf$coi), 10)
})


testthat::test_that("migration matrix with very low migration probability works", {
  # Test with very low migration probability
  migr_mat_low <- matrix(c(0.999, 0.001,
                           0.001, 0.999),
                         nrow = 2, byrow = TRUE)

  swf <- polySimIBD::sim_swf(
    pos = c(0, 1000, 2000),
    N = c(5, 5),
    m = c(0.5, 0.5),
    mean_coi = c(2, 2),
    migr_mat = migr_mat_low,
    rho = 1e-3,
    tlim = 3
  )

  testthat::expect_s3_class(swf, "swfsim")
})


testthat::test_that("migration matrix normalization handles rates correctly", {
  # Test that rate matrix gets converted to probabilities
  migr_mat_rates <- matrix(c(1, 2,
                             3, 1),
                           nrow = 2, byrow = TRUE)

  # Should not error - function should normalize
  swf <- polySimIBD::sim_swf(
    pos = c(0, 1000),
    N = c(3, 3),
    m = c(0.5, 0.5),
    mean_coi = c(2, 2),
    migr_mat = migr_mat_rates,
    rho = 1e-2,
    tlim = 2
  )

  testthat::expect_s3_class(swf, "swfsim")
})


testthat::test_that("large migration matrix works", {
  # Test with 5x5 migration matrix
  set.seed(42)
  n_demes <- 5
  migr_mat_large <- matrix(runif(n_demes^2), nrow = n_demes)
  # Normalize to probabilities
  migr_mat_large <- migr_mat_large / rowSums(migr_mat_large)

  swf <- polySimIBD::sim_swf(
    pos = c(0, 500, 1000),
    N = rep(3, n_demes),
    m = rep(0.3, n_demes),
    mean_coi = rep(2, n_demes),
    migr_mat = migr_mat_large,
    rho = 1e-2,
    tlim = 2
  )

  testthat::expect_s3_class(swf, "swfsim")
  testthat::expect_equal(length(swf$coi), 15) # 3 hosts per deme * 5 demes
})


# ============================================================================
# Tests for COI Edge Cases
# ============================================================================

testthat::test_that("very high mean COI works", {
  # Test with high mean COI
  swf <- polySimIBD::sim_swf(
    pos = c(0, 1000),
    N = 3,
    m = 0.5,
    mean_coi = 50,
    migr_mat = 1,
    rho = 1e-2,
    tlim = 2
  )

  testthat::expect_s3_class(swf, "swfsim")
  # COI should be high on average
  testthat::expect_gt(mean(swf$coi), 10)
})


testthat::test_that("low mean COI produces mostly COI=1", {
  set.seed(123)
  # Very low mean COI should produce mostly COI=1
  swf <- polySimIBD::sim_swf(
    pos = c(0, 1000),
    N = 20,
    m = 0.5,
    mean_coi = 0.1,
    migr_mat = 1,
    rho = 1e-2,
    tlim = 2
  )

  testthat::expect_s3_class(swf, "swfsim")
  # Most hosts should have COI=1
  testthat::expect_gt(sum(swf$coi == 1), 10)
})


# ============================================================================
# Tests for Position Vector Edge Cases
# ============================================================================

testthat::test_that("single locus works", {
  # Edge case: only one position
  swf <- polySimIBD::sim_swf(
    pos = c(0),
    N = 5,
    m = 0.5,
    mean_coi = 2,
    migr_mat = 1,
    rho = 1e-2,
    tlim = 2
  )

  testthat::expect_s3_class(swf, "swfsim")
  testthat::expect_equal(length(swf$pos), 1)
})


testthat::test_that("many loci work", {
  # Test with many loci
  pos_many <- seq(0, 10000, by = 10)

  swf <- polySimIBD::sim_swf(
    pos = pos_many,
    N = 3,
    m = 0.5,
    mean_coi = 2,
    migr_mat = 1,
    rho = 1e-2,
    tlim = 2
  )

  testthat::expect_s3_class(swf, "swfsim")
  testthat::expect_equal(length(swf$pos), length(pos_many))
})


testthat::test_that("uneven spacing in positions works", {
  # Test with irregularly spaced positions
  pos_uneven <- c(0, 10, 100, 150, 5000, 10000)

  swf <- polySimIBD::sim_swf(
    pos = pos_uneven,
    N = 3,
    m = 0.5,
    mean_coi = 2,
    migr_mat = 1,
    rho = 1e-2,
    tlim = 2
  )

  testthat::expect_s3_class(swf, "swfsim")
})


# ============================================================================
# Tests for Error Handling - Invalid Inputs
# ============================================================================

testthat::test_that("error when pos is not increasing", {
  testthat::expect_error(
    polySimIBD::sim_swf(
      pos = c(0, 1000, 500), # Not increasing!
      N = 5,
      m = 0.5,
      mean_coi = 2,
      migr_mat = 1,
      rho = 1e-2,
      tlim = 2
    ),
    "increasing"
  )
})


testthat::test_that("error when rho is exactly 0", {
  testthat::expect_error(
    polySimIBD::sim_swf(
      pos = c(0, 1000),
      N = 5,
      m = 0.5,
      mean_coi = 2,
      migr_mat = 1,
      rho = 0, # Invalid!
      tlim = 2
    ),
    "greater than 0"
  )
})


testthat::test_that("error when rho is exactly 1", {
  testthat::expect_error(
    polySimIBD::sim_swf(
      pos = c(0, 1000),
      N = 5,
      m = 0.5,
      mean_coi = 2,
      migr_mat = 1,
      rho = 1, # Invalid!
      tlim = 2
    ),
    "less than 1"
  )
})


testthat::test_that("error when N is zero", {
  testthat::expect_error(
    polySimIBD::sim_swf(
      pos = c(0, 1000),
      N = 0, # Invalid!
      m = 0.5,
      mean_coi = 2,
      migr_mat = 1,
      rho = 1e-2,
      tlim = 2
    )
  )
})


testthat::test_that("error when N is negative", {
  testthat::expect_error(
    polySimIBD::sim_swf(
      pos = c(0, 1000),
      N = -5, # Invalid!
      m = 0.5,
      mean_coi = 2,
      migr_mat = 1,
      rho = 1e-2,
      tlim = 2
    )
  )
})


testthat::test_that("error when m is negative", {
  testthat::expect_error(
    polySimIBD::sim_swf(
      pos = c(0, 1000),
      N = 5,
      m = -0.5, # Invalid!
      mean_coi = 2,
      migr_mat = 1,
      rho = 1e-2,
      tlim = 2
    ),
    "greater than or equal to 0"
  )
})


testthat::test_that("error when m is greater than 1", {
  testthat::expect_error(
    polySimIBD::sim_swf(
      pos = c(0, 1000),
      N = 5,
      m = 1.5, # Invalid!
      mean_coi = 2,
      migr_mat = 1,
      rho = 1e-2,
      tlim = 2
    ),
    "less than or equal to 1"
  )
})


testthat::test_that("error when mean_coi is zero", {
  testthat::expect_error(
    polySimIBD::sim_swf(
      pos = c(0, 1000),
      N = 5,
      m = 0.5,
      mean_coi = 0, # Invalid!
      migr_mat = 1,
      rho = 1e-2,
      tlim = 2
    )
  )
})


testthat::test_that("error when mean_coi is negative", {
  testthat::expect_error(
    polySimIBD::sim_swf(
      pos = c(0, 1000),
      N = 5,
      m = 0.5,
      mean_coi = -2, # Invalid!
      migr_mat = 1,
      rho = 1e-2,
      tlim = 2
    )
  )
})


testthat::test_that("error when migration matrix dimensions mismatch", {
  migr_mat_wrong <- matrix(c(0.5, 0.5,
                             0.5, 0.5),
                           nrow = 2, byrow = TRUE)

  testthat::expect_error(
    polySimIBD::sim_swf(
      pos = c(0, 1000),
      N = c(5, 5, 5), # 3 demes but matrix is 2x2!
      m = c(0.5, 0.5, 0.5),
      mean_coi = c(2, 2, 2),
      migr_mat = migr_mat_wrong,
      rho = 1e-2,
      tlim = 2
    ),
    "dimensions"
  )
})


# ============================================================================
# Tests for Weighted IBD Calculations
# ============================================================================

testthat::test_that("get_within_ibd works with NULL weights (default)", {
  set.seed(456)
  swf <- polySimIBD::sim_swf(
    pos = c(0, 1000, 2000, 3000),
    N = 5,
    m = 0.5,
    mean_coi = 8,
    migr_mat = 1,
    rho = 1e-2,
    tlim = 3
  )

  # Default behavior with NULL weights
  wibd <- polySimIBD::get_within_ibd(
    swf = swf,
    host_index = 1,
    weight_loci = NULL
  )

  testthat::expect_type(wibd, "double")
  testthat::expect_length(wibd, 1)
  testthat::expect_gte(wibd, 0)
})


testthat::test_that("get_within_ibd works with custom weights", {
  set.seed(457)
  N <- c(3, 2)
  m <- c(0.5, 0.5)
  mean_coi <- c(10, 10)
  migr_dist_mat <- matrix(c(0.7, 0.3,
                            0.3, 0.7),
                          nrow = 2, byrow = TRUE)

  swf <- polySimIBD::sim_swf(
    pos = c(0, 250, 500, 750, 1000),
    N = N,
    m = m,
    mean_coi = mean_coi,
    migr_mat = migr_dist_mat,
    rho = 1e-2,
    tlim = 3
  )

  # Calculate weights based on loci count
  n_loci <- length(swf$pos)
  weights <- rep(1 / n_loci, n_loci)

  # Should work without error
  wibd <- polySimIBD::get_within_ibd(
    swf = swf,
    host_index = 1,
    weight_loci = weights
  )

  testthat::expect_type(wibd, "double")
  testthat::expect_length(wibd, 1)
  testthat::expect_gte(wibd, 0)
})


testthat::test_that("get_bvibd works with NULL weights (default)", {
  set.seed(788)
  swf <- polySimIBD::sim_swf(
    pos = c(0, 1000, 2000),
    N = c(3, 3),
    m = c(0.5, 0.5),
    mean_coi = c(3, 3),
    migr_mat = matrix(c(0.8, 0.2, 0.2, 0.8), nrow = 2),
    rho = 1e-2,
    tlim = 3
  )

  # Default behavior with NULL weights
  bvibd <- polySimIBD::get_bvibd(
    swf = swf,
    host_index = c(1, 4)
  )

  testthat::expect_type(bvibd, "double")
  testthat::expect_length(bvibd, 1)
  testthat::expect_gte(bvibd, 0)
  testthat::expect_lte(bvibd, 1)
})


testthat::test_that("get_bvibd works with custom weights", {
  set.seed(789)
  N <- c(3, 2)
  m <- c(0.5, 0.5)
  mean_coi <- c(5, 5)
  migr_dist_mat <- matrix(c(0.8, 0.2,
                            0.2, 0.8),
                          nrow = 2, byrow = TRUE)

  swf <- polySimIBD::sim_swf(
    pos = c(0, 250, 500, 750, 1000),
    N = N,
    m = m,
    mean_coi = mean_coi,
    migr_mat = migr_dist_mat,
    rho = 1e-2,
    tlim = 3
  )

  # Calculate weights - note bvibd expects n_loci weights (not n_loci-1)
  n_loci <- length(swf$pos)
  weights <- rep(1 / n_loci, n_loci)

  # Should work without error
  bvibd <- polySimIBD::get_bvibd(
    swf = swf,
    host_index = c(1, 4),
    weight_loci = weights
  )

  testthat::expect_type(bvibd, "double")
  testthat::expect_length(bvibd, 1)
  testthat::expect_gte(bvibd, 0)
  testthat::expect_lte(bvibd, 1)
})


# ============================================================================
# Tests for ARG Edge Cases
# ============================================================================

testthat::test_that("get_arg works with single host", {
  swf <- polySimIBD::sim_swf(
    pos = c(0, 1000, 2000),
    N = 5,
    m = 0.5,
    mean_coi = 3,
    migr_mat = 1,
    rho = 1e-2,
    tlim = 3
  )

  # Get ARG for single host
  arg <- polySimIBD::get_arg(swf = swf, host_index = 1)

  testthat::expect_s3_class(arg, "argraph")
  testthat::expect_type(arg, "list")
})


testthat::test_that("get_arg works with all hosts", {
  swf <- polySimIBD::sim_swf(
    pos = c(0, 1000, 2000),
    N = 5,
    m = 0.5,
    mean_coi = 2,
    migr_mat = 1,
    rho = 1e-2,
    tlim = 3
  )

  # Get ARG for all hosts (default)
  arg <- polySimIBD::get_arg(swf = swf)

  testthat::expect_s3_class(arg, "argraph")
  testthat::expect_type(arg, "list")
})


testthat::test_that("get_arg works with specific haplotypes", {
  swf <- polySimIBD::sim_swf(
    pos = c(0, 1000, 2000),
    N = 3,
    m = 0.5,
    mean_coi = 5,
    migr_mat = 1,
    rho = 1e-2,
    tlim = 3
  )

  # Get ARG for specific haplotypes
  arg <- polySimIBD::get_arg(
    swf = swf,
    host_index = c(1, 2),
    haplo_index = list(c(1, 2), c(1))
  )

  testthat::expect_s3_class(arg, "argraph")
})


# ============================================================================
# Tests for Boundary Conditions
# ============================================================================

testthat::test_that("tlim = 2 (minimum sensible value) works", {
  swf <- polySimIBD::sim_swf(
    pos = c(0, 1000),
    N = 5,
    m = 0.5,
    mean_coi = 2,
    migr_mat = 1,
    rho = 1e-2,
    tlim = 2 # Minimum
  )

  testthat::expect_s3_class(swf, "swfsim")
})


testthat::test_that("very small rho works", {
  # Very small recombination rate
  swf <- polySimIBD::sim_swf(
    pos = c(0, 1000),
    N = 5,
    m = 0.5,
    mean_coi = 2,
    migr_mat = 1,
    rho = 1e-10,
    tlim = 2
  )

  testthat::expect_s3_class(swf, "swfsim")
})


testthat::test_that("rho close to 1 works", {
  # Very high recombination rate (but not exactly 1)
  swf <- polySimIBD::sim_swf(
    pos = c(0, 1000),
    N = 5,
    m = 0.5,
    mean_coi = 2,
    migr_mat = 1,
    rho = 0.999,
    tlim = 2
  )

  testthat::expect_s3_class(swf, "swfsim")
})


testthat::test_that("m = 0 (no migration) works", {
  swf <- polySimIBD::sim_swf(
    pos = c(0, 1000),
    N = 5,
    m = 0, # No migration
    mean_coi = 2,
    migr_mat = 1,
    rho = 1e-2,
    tlim = 2
  )

  testthat::expect_s3_class(swf, "swfsim")
})


testthat::test_that("m = 1 (complete migration) works", {
  swf <- polySimIBD::sim_swf(
    pos = c(0, 1000),
    N = 5,
    m = 1, # Complete migration
    mean_coi = 2,
    migr_mat = 1,
    rho = 1e-2,
    tlim = 2
  )

  testthat::expect_s3_class(swf, "swfsim")
})


# ============================================================================
# Tests for get_effective_coi Edge Cases
# ============================================================================

testthat::test_that("get_effective_coi returns correct length vector", {
  swf <- polySimIBD::sim_swf(
    pos = c(0, 500, 1000, 1500, 2000),
    N = 3,
    m = 0.5,
    mean_coi = 5,
    migr_mat = 1,
    rho = 1e-2,
    tlim = 3
  )

  eff_coi <- polySimIBD::get_effective_coi(swf = swf, host_index = 1)

  testthat::expect_length(eff_coi, length(swf$pos))
  testthat::expect_gte(min(eff_coi), 1)
  testthat::expect_lte(max(eff_coi), swf$coi[1])
})


# ============================================================================
# Tests for subset_bvtree Edge Cases
# ============================================================================

testthat::test_that("subset_bvtree works with single element", {
  swf <- polySimIBD::sim_swf(
    pos = c(0, 1000),
    N = 5,
    m = 0.5,
    mean_coi = 3,
    migr_mat = 1,
    rho = 1e-2,
    tlim = 3
  )

  arg <- polySimIBD::get_arg(swf = swf, host_index = c(1, 2))

  # This should error because needs > 1 element
  testthat::expect_error(
    polySimIBD::subset_bvtree(bvtree = arg[[1]], s = c(1)),
    "gr"
  )
})


testthat::test_that("subset_bvtree works with all elements", {
  swf <- polySimIBD::sim_swf(
    pos = c(0, 1000),
    N = 3,
    m = 0.5,
    mean_coi = 3,
    migr_mat = 1,
    rho = 1e-2,
    tlim = 3
  )

  arg <- polySimIBD::get_arg(swf = swf, host_index = c(1, 2))
  n_elem <- length(arg[[1]]@c)

  # Subset to all elements
  bvsub <- polySimIBD::subset_bvtree(bvtree = arg[[1]], s = 1:n_elem)

  testthat::expect_s4_class(bvsub, "bvtree")
})
