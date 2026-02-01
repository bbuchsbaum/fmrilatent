library(testthat)

# Tests for voxel_subset_to_gsp function in R/graph_bridge.R
#
# NOTE: There appears to be a bug in the source code. The function documentation
# states it accepts "LogicalNeuroVol or 3D logical array", but the implementation
# on line 19 uses `which(as.logical(mask_arr), arr.ind = TRUE)`. The as.logical()
# call flattens plain R arrays to vectors, causing arr.ind = TRUE to return a
# simple integer vector instead of a coordinate matrix. This results in the error:
# "incorrect number of dimensions"
#
# The function works correctly with LogicalNeuroVol objects because neuroim2's
# as.logical method for LogicalNeuroVol preserves array structure.
#
# Suggested fix: Change line 19 from:
#   idx_all <- which(as.logical(mask_arr), arr.ind = TRUE)
# to:
#   idx_all <- which(mask_arr, arr.ind = TRUE)
#
# or ensure mask_arr is coerced to logical array while preserving dimensions:
#   storage.mode(mask_arr) <- "logical"
#   idx_all <- which(mask_arr, arr.ind = TRUE)

# -----------------------------------------------------------------------------
# Tests with LogicalNeuroVol (the currently working case)
# -----------------------------------------------------------------------------

test_that("voxel_subset_to_gsp returns gsp_graph object with LogicalNeuroVol", {
  skip_if_not_installed("rgsp")
  skip_if_not_installed("neuroim2")

  mask_arr <- array(TRUE, dim = c(3, 3, 3))
  spc <- neuroim2::NeuroSpace(dim(mask_arr))
  mask_vol <- neuroim2::LogicalNeuroVol(mask_arr, spc)
  voxel_indices <- 1:10

  result <- voxel_subset_to_gsp(mask_vol, voxel_indices, k_neighbors = 3L)

  expect_s3_class(result, "gsp_graph")
})

test_that("voxel_subset_to_gsp works with all TRUE LogicalNeuroVol mask", {
  skip_if_not_installed("rgsp")
  skip_if_not_installed("neuroim2")

  mask_arr <- array(TRUE, dim = c(2, 2, 2))
  spc <- neuroim2::NeuroSpace(dim(mask_arr))
  mask_vol <- neuroim2::LogicalNeuroVol(mask_arr, spc)
  n_voxels <- sum(mask_arr)
  voxel_indices <- seq_len(n_voxels)

  result <- voxel_subset_to_gsp(mask_vol, voxel_indices, k_neighbors = 3L)

  expect_s3_class(result, "gsp_graph")
})

test_that("voxel_subset_to_gsp works with partial LogicalNeuroVol mask", {
  skip_if_not_installed("rgsp")
  skip_if_not_installed("neuroim2")

  mask_arr <- array(FALSE, dim = c(3, 3, 3))
  mask_arr[1:2, 1:2, 1:2] <- TRUE  # 8 voxels are TRUE
  spc <- neuroim2::NeuroSpace(dim(mask_arr))
  mask_vol <- neuroim2::LogicalNeuroVol(mask_arr, spc)

  n_true <- sum(mask_arr)
  voxel_indices <- seq_len(n_true)

  result <- voxel_subset_to_gsp(mask_vol, voxel_indices, k_neighbors = 2L)

  expect_s3_class(result, "gsp_graph")
})

test_that("voxel_subset_to_gsp works with subset of voxel indices", {
  skip_if_not_installed("rgsp")
  skip_if_not_installed("neuroim2")

  mask_arr <- array(TRUE, dim = c(3, 3, 3))
  spc <- neuroim2::NeuroSpace(dim(mask_arr))
  mask_vol <- neuroim2::LogicalNeuroVol(mask_arr, spc)
  # Select only a subset of voxels
  voxel_indices <- c(1, 3, 5, 7, 9)

  result <- voxel_subset_to_gsp(mask_vol, voxel_indices, k_neighbors = 2L)

  expect_s3_class(result, "gsp_graph")
})

test_that("voxel_subset_to_gsp respects different k_neighbors values", {
  skip_if_not_installed("rgsp")
  skip_if_not_installed("neuroim2")

  mask_arr <- array(TRUE, dim = c(3, 3, 2))
  spc <- neuroim2::NeuroSpace(dim(mask_arr))
  mask_vol <- neuroim2::LogicalNeuroVol(mask_arr, spc)
  voxel_indices <- seq_len(sum(mask_arr))

  # Test with different k values
  result_k2 <- voxel_subset_to_gsp(mask_vol, voxel_indices, k_neighbors = 2L)
  result_k4 <- voxel_subset_to_gsp(mask_vol, voxel_indices, k_neighbors = 4L)
  result_k6 <- voxel_subset_to_gsp(mask_vol, voxel_indices, k_neighbors = 6L)

  expect_s3_class(result_k2, "gsp_graph")
  expect_s3_class(result_k4, "gsp_graph")
  expect_s3_class(result_k6, "gsp_graph")
})

test_that("voxel_subset_to_gsp default k_neighbors is 6", {
  skip_if_not_installed("rgsp")
  skip_if_not_installed("neuroim2")

  mask_arr <- array(TRUE, dim = c(3, 3, 3))
  spc <- neuroim2::NeuroSpace(dim(mask_arr))
  mask_vol <- neuroim2::LogicalNeuroVol(mask_arr, spc)
  voxel_indices <- seq_len(sum(mask_arr))

  # Call without specifying k_neighbors
  result <- voxel_subset_to_gsp(mask_vol, voxel_indices)

  expect_s3_class(result, "gsp_graph")
})

test_that("voxel_subset_to_gsp works with 2D-like mask (z=1)", {
  skip_if_not_installed("rgsp")
  skip_if_not_installed("neuroim2")

  # A "flat" 3D array (effectively 2D)
  mask_arr <- array(TRUE, dim = c(4, 4, 1))
  spc <- neuroim2::NeuroSpace(dim(mask_arr))
  mask_vol <- neuroim2::LogicalNeuroVol(mask_arr, spc)
  voxel_indices <- seq_len(sum(mask_arr))

  result <- voxel_subset_to_gsp(mask_vol, voxel_indices, k_neighbors = 4L)

  expect_s3_class(result, "gsp_graph")
})

test_that("voxel_subset_to_gsp errors with single voxel (rgsp requires >=2)", {
  skip_if_not_installed("rgsp")
  skip_if_not_installed("neuroim2")

  mask_arr <- array(TRUE, dim = c(2, 2, 2))
  spc <- neuroim2::NeuroSpace(dim(mask_arr))
  mask_vol <- neuroim2::LogicalNeuroVol(mask_arr, spc)
  voxel_indices <- 1L

  # rgsp::graph_knn requires at least 2 coordinates
  expect_error(
    voxel_subset_to_gsp(mask_vol, voxel_indices, k_neighbors = 1L),
    "at least 2"
  )
})

test_that("voxel_subset_to_gsp works with two voxels (minimum for rgsp)", {
  skip_if_not_installed("rgsp")
  skip_if_not_installed("neuroim2")

  mask_arr <- array(TRUE, dim = c(2, 2, 2))
  spc <- neuroim2::NeuroSpace(dim(mask_arr))
  mask_vol <- neuroim2::LogicalNeuroVol(mask_arr, spc)
  voxel_indices <- 1:2

  # Two voxels should work
  result <- voxel_subset_to_gsp(mask_vol, voxel_indices, k_neighbors = 1L)

  expect_s3_class(result, "gsp_graph")
})

test_that("voxel_subset_to_gsp works with non-contiguous voxel indices", {
  skip_if_not_installed("rgsp")
  skip_if_not_installed("neuroim2")

  mask_arr <- array(TRUE, dim = c(4, 4, 4))
  spc <- neuroim2::NeuroSpace(dim(mask_arr))
  mask_vol <- neuroim2::LogicalNeuroVol(mask_arr, spc)
  # Select non-contiguous voxel indices
  voxel_indices <- c(1, 10, 20, 30, 40, 50)

  result <- voxel_subset_to_gsp(mask_vol, voxel_indices, k_neighbors = 3L)

  expect_s3_class(result, "gsp_graph")
})

test_that("voxel_subset_to_gsp works with large partial mask", {
  skip_if_not_installed("rgsp")
  skip_if_not_installed("neuroim2")

  mask_arr <- array(FALSE, dim = c(4, 4, 4))
  mask_arr[2:3, 2:3, 2:3] <- TRUE  # 8 TRUE voxels in center
  spc <- neuroim2::NeuroSpace(dim(mask_arr))
  mask_vol <- neuroim2::LogicalNeuroVol(mask_arr, spc)

  n_true <- sum(mask_arr)
  voxel_indices <- seq_len(n_true)

  result <- voxel_subset_to_gsp(mask_vol, voxel_indices, k_neighbors = 4L)

  expect_s3_class(result, "gsp_graph")
})

# -----------------------------------------------------------------------------
# Tests for error conditions
# -----------------------------------------------------------------------------

test_that("voxel_subset_to_gsp errors when rgsp not available", {
  # Mock the situation where rgsp is not installed
  # We can't truly test this if rgsp is installed, so we skip if it is
  skip_if(requireNamespace("rgsp", quietly = TRUE),
          "rgsp is installed, cannot test missing package error")

  mask <- array(TRUE, dim = c(2, 2, 2))
  voxel_indices <- 1:4

  expect_error(
    voxel_subset_to_gsp(mask, voxel_indices),
    "rgsp not installed"
  )
})

test_that("voxel_subset_to_gsp errors with NULL mask", {
  skip_if_not_installed("rgsp")

  expect_error(
    voxel_subset_to_gsp(NULL, 1:4),
    "mask must be array-like or LogicalNeuroVol"
  )
})

test_that("voxel_subset_to_gsp errors with function as mask", {
  skip_if_not_installed("rgsp")

  # A function cannot be converted to array
  invalid_mask <- function(x) x

  expect_error(
    voxel_subset_to_gsp(invalid_mask, 1:4),
    "mask must be array-like or LogicalNeuroVol"
  )
})

test_that("voxel_subset_to_gsp errors with out-of-bounds voxel indices", {
  skip_if_not_installed("rgsp")
  skip_if_not_installed("neuroim2")

  mask_arr <- array(TRUE, dim = c(2, 2, 2))  # 8 TRUE voxels
  spc <- neuroim2::NeuroSpace(dim(mask_arr))
  mask_vol <- neuroim2::LogicalNeuroVol(mask_arr, spc)
  # voxel_indices that exceed the number of TRUE voxels
  voxel_indices <- c(1, 2, 100)

  expect_error(voxel_subset_to_gsp(mask_vol, voxel_indices))
})

test_that("voxel_subset_to_gsp errors with mixed positive/negative voxel indices", {
  skip_if_not_installed("rgsp")
  skip_if_not_installed("neuroim2")

  mask_arr <- array(TRUE, dim = c(2, 2, 2))
  spc <- neuroim2::NeuroSpace(dim(mask_arr))
  mask_vol <- neuroim2::LogicalNeuroVol(mask_arr, spc)

  # R subsetting does not allow mixing positive and negative indices
  # (except when positive are 0)
  expect_error(voxel_subset_to_gsp(mask_vol, c(-1, 2, 3), k_neighbors = 1L))
})

test_that("voxel_subset_to_gsp handles zero index (ignored by R)", {
  skip_if_not_installed("rgsp")
  skip_if_not_installed("neuroim2")

  mask_arr <- array(TRUE, dim = c(2, 2, 2))
  spc <- neuroim2::NeuroSpace(dim(mask_arr))
  mask_vol <- neuroim2::LogicalNeuroVol(mask_arr, spc)

  # R subsetting with 0 ignores that index
  # c(0, 1, 2) would select indices 1 and 2 only
  result <- voxel_subset_to_gsp(mask_vol, c(0, 1, 2), k_neighbors = 1L)
  expect_s3_class(result, "gsp_graph")
})

# -----------------------------------------------------------------------------
# Tests documenting the plain array bug (expected to fail until bug is fixed)
# These tests document the discrepancy between documentation and behavior
# -----------------------------------------------------------------------------

test_that("voxel_subset_to_gsp with plain 3D logical array works", {
  skip_if_not_installed("rgsp")

  mask <- array(TRUE, dim = c(3, 3, 3))
  voxel_indices <- 1:10

  result <- voxel_subset_to_gsp(mask, voxel_indices, k_neighbors = 3L)

  expect_s3_class(result, "gsp_graph")
})
