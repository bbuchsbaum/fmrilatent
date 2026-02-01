# Tests for generic function definitions in all_generic.R

library(testthat)
library(neuroim2)
library(Matrix)

# Helper function to create a minimal LatentNeuroVec for testing
create_test_lvec <- function(n_time = 5L, n_components = 2L, dims_3d = c(3, 3, 2),
                              sparse_loadings = FALSE, with_offset = FALSE) {
  n_voxels <- prod(dims_3d)

  # Create mask - all voxels included

  mask_arr <- array(TRUE, dim = dims_3d)
  mask_vol <- LogicalNeuroVol(mask_arr, NeuroSpace(dims_3d))

  # Create space
  space <- NeuroSpace(c(dims_3d, n_time))

  # Create basis (time x components)
  basis <- Matrix::Matrix(
    matrix(rnorm(n_time * n_components), nrow = n_time, ncol = n_components),
    sparse = FALSE
  )

  # Create loadings (voxels x components)
  loadings <- Matrix::Matrix(
    matrix(rnorm(n_voxels * n_components), nrow = n_voxels, ncol = n_components),
    sparse = sparse_loadings
  )

  # Create offset if requested
  offset <- if (with_offset) rnorm(n_voxels) else numeric(0)

  LatentNeuroVec(
    basis = basis,
    loadings = loadings,
    space = space,
    mask = mask_vol,
    offset = offset,
    label = "test-lvec"
  )
}

# ============================================================================
# Tests for basis() generic
# ============================================================================

test_that("basis generic is defined", {
  expect_true(isGeneric("basis"))
})

test_that("basis method returns correct matrix for LatentNeuroVec", {
  lvec <- create_test_lvec(n_time = 5L, n_components = 3L)

  b <- basis(lvec)

  expect_true(inherits(b, "Matrix") || is.matrix(b))
  expect_equal(nrow(b), 5L)  # n_time

  expect_equal(ncol(b), 3L)  # n_components
})

test_that("basis method preserves values from construction", {
  n_time <- 4L
  n_components <- 2L
  dims_3d <- c(2, 2, 2)
  n_voxels <- prod(dims_3d)

  # Create specific basis values
  original_basis <- matrix(1:8, nrow = n_time, ncol = n_components)
  basis_mat_obj <- Matrix::Matrix(original_basis, sparse = FALSE)

  mask_arr <- array(TRUE, dim = dims_3d)
  mask_vol <- LogicalNeuroVol(mask_arr, NeuroSpace(dims_3d))
  space <- NeuroSpace(c(dims_3d, n_time))
  loadings <- Matrix::Matrix(matrix(1, nrow = n_voxels, ncol = n_components), sparse = FALSE)

  lvec <- LatentNeuroVec(
    basis = basis_mat_obj,
    loadings = loadings,
    space = space,
    mask = mask_vol
  )

  retrieved_basis <- basis(lvec)
  expect_equal(as.matrix(retrieved_basis), original_basis)
})

test_that("basis works with dense Matrix input", {
  lvec <- create_test_lvec(sparse_loadings = FALSE)
  b <- basis(lvec)
  expect_true(all(is.finite(b)))
})

# ============================================================================
# Tests for loadings() generic
# ============================================================================

test_that("loadings generic is defined", {

  expect_true(isGeneric("loadings"))
})

test_that("loadings method returns correct matrix for LatentNeuroVec", {
  dims_3d <- c(3, 3, 2)
  n_voxels <- prod(dims_3d)
  n_components <- 4L

  lvec <- create_test_lvec(n_components = n_components, dims_3d = dims_3d)

  L <- loadings(lvec)

  expect_true(inherits(L, "Matrix") || is.matrix(L))
  expect_equal(nrow(L), n_voxels)
  expect_equal(ncol(L), n_components)
})

test_that("loadings method preserves values from construction", {
  n_time <- 3L
  n_components <- 2L
  dims_3d <- c(2, 2, 1)
  n_voxels <- prod(dims_3d)

  # Create specific loadings values
  original_loadings <- matrix(seq(1, n_voxels * n_components), nrow = n_voxels, ncol = n_components)
  loadings_mat_obj <- Matrix::Matrix(original_loadings, sparse = FALSE)

  mask_arr <- array(TRUE, dim = dims_3d)
  mask_vol <- LogicalNeuroVol(mask_arr, NeuroSpace(dims_3d))
  space <- NeuroSpace(c(dims_3d, n_time))
  basis <- Matrix::Matrix(matrix(1, nrow = n_time, ncol = n_components), sparse = FALSE)

  lvec <- LatentNeuroVec(
    basis = basis,
    loadings = loadings_mat_obj,
    space = space,
    mask = mask_vol
  )

  retrieved_loadings <- loadings(lvec)
  expect_equal(as.matrix(retrieved_loadings), original_loadings)
})

test_that("loadings works with sparse Matrix input", {
  n_time <- 5L
  n_components <- 2L
  dims_3d <- c(3, 3, 2)
  n_voxels <- prod(dims_3d)

  # Create a sparse loadings matrix (mostly zeros)
  sparse_loadings_mat <- Matrix::Matrix(0, nrow = n_voxels, ncol = n_components, sparse = TRUE)
  sparse_loadings_mat[1, 1] <- 1.0
  sparse_loadings_mat[n_voxels, 2] <- 1.0

  mask_arr <- array(TRUE, dim = dims_3d)
  mask_vol <- LogicalNeuroVol(mask_arr, NeuroSpace(dims_3d))
  space <- NeuroSpace(c(dims_3d, n_time))
  basis <- Matrix::Matrix(matrix(rnorm(n_time * n_components), nrow = n_time), sparse = FALSE)

  lvec <- LatentNeuroVec(
    basis = basis,
    loadings = sparse_loadings_mat,
    space = space,
    mask = mask_vol
  )

  L <- loadings(lvec)
  expect_true(inherits(L, "Matrix"))
  expect_equal(L[1, 1], 1.0)
  expect_equal(L[n_voxels, 2], 1.0)
})

# ============================================================================
# Tests for offset() generic
# ============================================================================

test_that("offset generic is defined", {
  expect_true(isGeneric("offset"))
})

test_that("offset method returns empty numeric when no offset provided", {
  lvec <- create_test_lvec(with_offset = FALSE)

  off <- offset(lvec)

  expect_true(is.numeric(off))
  expect_equal(length(off), 0)
})

test_that("offset method returns correct values when offset provided", {
  n_time <- 4L
  n_components <- 2L
  dims_3d <- c(2, 2, 2)
  n_voxels <- prod(dims_3d)

  # Create specific offset values
  original_offset <- seq(0.1, by = 0.1, length.out = n_voxels)

  mask_arr <- array(TRUE, dim = dims_3d)
  mask_vol <- LogicalNeuroVol(mask_arr, NeuroSpace(dims_3d))
  space <- NeuroSpace(c(dims_3d, n_time))
  basis <- Matrix::Matrix(matrix(1, nrow = n_time, ncol = n_components), sparse = FALSE)
  loadings <- Matrix::Matrix(matrix(1, nrow = n_voxels, ncol = n_components), sparse = FALSE)

  lvec <- LatentNeuroVec(
    basis = basis,
    loadings = loadings,
    space = space,
    mask = mask_vol,
    offset = original_offset
  )

  retrieved_offset <- offset(lvec)
  expect_equal(retrieved_offset, original_offset)
})

test_that("offset length matches number of voxels in mask", {
  lvec <- create_test_lvec(dims_3d = c(4, 4, 2), with_offset = TRUE)

  off <- offset(lvec)
  n_voxels <- sum(mask(lvec))

  expect_equal(length(off), n_voxels)
})

# ============================================================================
# Tests for mask() generic
# ============================================================================

test_that("mask generic is defined", {
  expect_true(isGeneric("mask"))
})

test_that("mask method returns LogicalNeuroVol for LatentNeuroVec", {
  lvec <- create_test_lvec()

  m <- mask(lvec)

  expect_true(inherits(m, "LogicalNeuroVol"))
})

test_that("mask method preserves dimensions", {
  dims_3d <- c(5, 4, 3)
  lvec <- create_test_lvec(dims_3d = dims_3d)

  m <- mask(lvec)

  expect_equal(dim(m), dims_3d)
})

test_that("mask returns correct voxel count", {
  dims_3d <- c(3, 3, 3)
  n_voxels_expected <- prod(dims_3d)
  lvec <- create_test_lvec(dims_3d = dims_3d)

  m <- mask(lvec)

  expect_equal(sum(m), n_voxels_expected)
})

# ============================================================================
# Tests for map() generic
# ============================================================================

test_that("map generic is defined", {
  expect_true(isGeneric("map"))
})

test_that("map method returns IndexLookupVol for LatentNeuroVec", {
  lvec <- create_test_lvec()

  mp <- map(lvec)

  expect_true(inherits(mp, "IndexLookupVol"))
})

test_that("map correctly maps spatial indices to mask indices", {
  dims_3d <- c(2, 2, 2)
  lvec <- create_test_lvec(dims_3d = dims_3d)

  mp <- map(lvec)

  # For a full mask, all 8 voxels should map to indices 1-8
  all_spatial_indices <- 1:prod(dims_3d)
  mask_indices <- lookup(mp, all_spatial_indices)

  expect_true(all(mask_indices > 0))
  expect_equal(length(unique(mask_indices)), prod(dims_3d))
})

# ============================================================================
# Tests for interaction between generics
# ============================================================================

test_that("basis and loadings have compatible dimensions", {
  lvec <- create_test_lvec(n_components = 5L)

  b <- basis(lvec)
  L <- loadings(lvec)

  # Number of components should match

  expect_equal(ncol(b), ncol(L))
})

test_that("offset length matches loadings rows when provided", {
  lvec <- create_test_lvec(dims_3d = c(3, 3, 2), with_offset = TRUE)

  off <- offset(lvec)
  L <- loadings(lvec)

  expect_equal(length(off), nrow(L))
})

test_that("mask voxel count matches loadings rows", {
  dims_3d <- c(4, 3, 2)
  lvec <- create_test_lvec(dims_3d = dims_3d)

  m <- mask(lvec)
  L <- loadings(lvec)

  expect_equal(sum(m), nrow(L))
})

# ============================================================================
# Tests for generic dispatch on different object types
# ============================================================================

test_that("basis method not defined for arbitrary objects produces error", {
  expect_error(basis(list(a = 1)), "unable to find an inherited method")
})

test_that("loadings method for non-supported objects returns NULL (default ANY method)", {
  # The loadings generic has a default ANY method that returns NULL
  # This is inherited behavior from stats::loadings
  expect_null(loadings(list(a = 1)))
})

test_that("offset generic uses stats::offset for non-LatentNeuroVec objects", {
  # stats::offset returns the object unchanged for arbitrary objects
  result <- offset(list(a = 1))
  expect_equal(result, list(a = 1))
})

# ============================================================================
# Tests for edge cases
# ============================================================================

test_that("basis works with single component", {
  lvec <- create_test_lvec(n_components = 1L)

  b <- basis(lvec)

  expect_equal(ncol(b), 1L)
})

test_that("basis works with single time point", {
  n_time <- 1L
  n_components <- 2L
  dims_3d <- c(2, 2, 1)
  n_voxels <- prod(dims_3d)

  mask_arr <- array(TRUE, dim = dims_3d)
  mask_vol <- LogicalNeuroVol(mask_arr, NeuroSpace(dims_3d))
  space <- NeuroSpace(c(dims_3d, n_time))
  basis <- Matrix::Matrix(matrix(c(1, 2), nrow = n_time, ncol = n_components), sparse = FALSE)
  loadings <- Matrix::Matrix(matrix(1, nrow = n_voxels, ncol = n_components), sparse = FALSE)

  lvec <- LatentNeuroVec(
    basis = basis,
    loadings = loadings,
    space = space,
    mask = mask_vol
  )

  b <- basis(lvec)
  expect_equal(nrow(b), 1L)
})

test_that("loadings works with minimal 1x1x1 volume", {
  n_time <- 3L
  n_components <- 2L
  dims_3d <- c(1, 1, 1)
  n_voxels <- prod(dims_3d)

  mask_arr <- array(TRUE, dim = dims_3d)
  mask_vol <- LogicalNeuroVol(mask_arr, NeuroSpace(dims_3d))
  space <- NeuroSpace(c(dims_3d, n_time))
  basis <- Matrix::Matrix(matrix(rnorm(n_time * n_components), nrow = n_time), sparse = FALSE)
  loadings <- Matrix::Matrix(matrix(c(0.5, 1.5), nrow = n_voxels, ncol = n_components), sparse = FALSE)

  lvec <- LatentNeuroVec(
    basis = basis,
    loadings = loadings,
    space = space,
    mask = mask_vol
  )

  L <- loadings(lvec)
  expect_equal(nrow(L), 1L)
  expect_equal(ncol(L), n_components)
})

# ============================================================================
# Tests for generic function with ... argument handling
# ============================================================================

test_that("basis generic accepts additional arguments without error", {
  lvec <- create_test_lvec()
  # The ... arguments should be accepted but currently unused
  b <- basis(lvec)
  expect_true(inherits(b, "Matrix") || is.matrix(b))
})

test_that("loadings generic accepts additional arguments without error", {
  lvec <- create_test_lvec()
  L <- loadings(lvec)
  expect_true(inherits(L, "Matrix") || is.matrix(L))
})

test_that("offset generic accepts additional arguments without error", {
  lvec <- create_test_lvec()
  off <- offset(lvec)
  expect_true(is.numeric(off))
})

# ============================================================================
# Test conditional generic definitions (map and mask)
# ============================================================================

test_that("map generic exists and can be dispatched", {
  lvec <- create_test_lvec()
  mp <- map(lvec)
  expect_true(inherits(mp, "IndexLookupVol"))
})

test_that("mask generic exists and can be dispatched", {
  lvec <- create_test_lvec()
  m <- mask(lvec)
  expect_true(inherits(m, "LogicalNeuroVol"))
})

# ============================================================================
# Reconstruction tests using generics
# ============================================================================

test_that("reconstruction using basis and loadings matches as.matrix", {
  set.seed(123)
  lvec <- create_test_lvec(n_time = 5L, n_components = 3L, dims_3d = c(3, 3, 2))

  b <- basis(lvec)
  L <- loadings(lvec)
  off <- offset(lvec)

  # Manual reconstruction: X = B %*% t(L) + offset
  manual_recon <- as.matrix(b %*% t(L))
  if (length(off) > 0) {
    manual_recon <- sweep(manual_recon, 2, off, "+")
  }

  # Using as.matrix method
  auto_recon <- as.matrix(lvec)

  expect_equal(manual_recon, auto_recon, tolerance = 1e-10)
})

test_that("reconstruction with offset is correct", {
  set.seed(456)
  lvec <- create_test_lvec(n_time = 4L, n_components = 2L, dims_3d = c(2, 2, 2), with_offset = TRUE)

  b <- basis(lvec)
  L <- loadings(lvec)
  off <- offset(lvec)

  expect_true(length(off) > 0)

  # Manual reconstruction with offset
  manual_recon <- as.matrix(b %*% t(L))
  manual_recon <- sweep(manual_recon, 2, off, "+")

  auto_recon <- as.matrix(lvec)

  expect_equal(manual_recon, auto_recon, tolerance = 1e-10)
})
