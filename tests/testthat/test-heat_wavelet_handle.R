# Tests for heat_wavelet_handle.R
# Tests the heat_wavelet_loadings_handle() function which creates LoadingsHandle
# objects for shared heat-wavelet spatial dictionaries.

library(testthat)

# All tests require rgsp package for heat wavelet computations
skip_if_not_installed("rgsp")

# Helper function to create a minimal ClusterReduction for testing
create_test_reduction <- function(mask_dim = c(2, 2, 2)) {
  mask <- array(TRUE, dim = mask_dim)
  map <- seq_len(sum(mask))
  make_cluster_reduction(mask, map)
}

# --- Constructor Tests ---

test_that("heat_wavelet_loadings_handle returns a LoadingsHandle object", {

  reduction <- create_test_reduction()
  spec <- basis_heat_wavelet(scales = 1, order = 10, threshold = 0)

  handle <- heat_wavelet_loadings_handle(reduction, spec)

  expect_s4_class(handle, "LoadingsHandle")
})

test_that("heat_wavelet_loadings_handle stores correct dimensions", {
  mask_dim <- c(2, 2, 2)
  n_vox <- prod(mask_dim)
  scales <- c(1, 2)

  reduction <- create_test_reduction(mask_dim)
  spec <- basis_heat_wavelet(scales = scales, order = 10, threshold = 0)

  handle <- heat_wavelet_loadings_handle(reduction, spec)

  # Dimensions should be (n_voxels, n_voxels * n_scales)
  expect_equal(handle@dim[1], n_vox)
  expect_equal(handle@dim[2], n_vox * length(scales))
})

test_that("heat_wavelet_loadings_handle sets kind to 'lifted'", {
  reduction <- create_test_reduction()
  spec <- basis_heat_wavelet(scales = 1, order = 8, threshold = 0)

  handle <- heat_wavelet_loadings_handle(reduction, spec)

  expect_equal(handle@kind, "lifted")
})

test_that("heat_wavelet_loadings_handle stores the default label", {
  reduction <- create_test_reduction()
  spec <- basis_heat_wavelet(scales = 1, order = 8, threshold = 0)

  handle <- heat_wavelet_loadings_handle(reduction, spec)

  expect_equal(handle@label, "heat-wavelet")
})

test_that("heat_wavelet_loadings_handle accepts custom label", {
  reduction <- create_test_reduction()
  spec <- basis_heat_wavelet(scales = 1, order = 8, threshold = 0)

  custom_label <- "my-custom-wavelet-basis"
  handle <- heat_wavelet_loadings_handle(reduction, spec, label = custom_label)

  expect_equal(handle@label, custom_label)
})

test_that("heat_wavelet_loadings_handle stores spec with correct family", {
  reduction <- create_test_reduction()
  spec <- basis_heat_wavelet(scales = c(1, 2), order = 10, threshold = 1e-6)

  handle <- heat_wavelet_loadings_handle(reduction, spec)

  expect_equal(handle@spec$family, "heat_wavelet")
  expect_identical(handle@spec$reduction, reduction)
  expect_identical(handle@spec$basis_spec, spec)
  expect_null(handle@spec$data)
})

# --- ID Handling Tests ---

test_that("heat_wavelet_loadings_handle generates random ID when not provided", {
  reduction <- create_test_reduction()
  spec <- basis_heat_wavelet(scales = 1, order = 8, threshold = 0)

  handle <- heat_wavelet_loadings_handle(reduction, spec)

  expect_true(nchar(handle@id) > 0)
  expect_true(grepl("^heat-wavelet-", handle@id))
})

test_that("heat_wavelet_loadings_handle uses provided ID", {
  reduction <- create_test_reduction()
  spec <- basis_heat_wavelet(scales = 1, order = 8, threshold = 0)

  custom_id <- "my-stable-test-id"
  handle <- heat_wavelet_loadings_handle(reduction, spec, id = custom_id)

  expect_equal(handle@id, custom_id)
})

test_that("heat_wavelet_loadings_handle generates unique IDs for different calls", {
  reduction <- create_test_reduction()
  spec <- basis_heat_wavelet(scales = 1, order = 8, threshold = 0)

  handle1 <- heat_wavelet_loadings_handle(reduction, spec)
  handle2 <- heat_wavelet_loadings_handle(reduction, spec)

  expect_false(handle1@id == handle2@id)
})

# --- Registry Tests ---

test_that("heat_wavelet_loadings_handle registers matrix in loadings registry", {
  reduction <- create_test_reduction()
  spec <- basis_heat_wavelet(scales = 1, order = 8, threshold = 0)

  handle <- heat_wavelet_loadings_handle(reduction, spec)

  expect_true(fmrilatent:::.latent_has_matrix(handle@id, type = "loadings"))
})

test_that("registered matrix can be retrieved via loadings_mat", {
  reduction <- create_test_reduction()
  spec <- basis_heat_wavelet(scales = 1, order = 8, threshold = 0)

  handle <- heat_wavelet_loadings_handle(reduction, spec)

  # loadings_mat should retrieve the registered matrix
  mat <- loadings_mat(handle)

  expect_true(inherits(mat, "Matrix") || is.matrix(mat))
  expect_equal(dim(mat), handle@dim)
})

test_that("registered matrix matches dimensions from handle", {
  mask_dim <- c(3, 3, 2)
  n_vox <- prod(mask_dim)
  scales <- c(1, 2, 4)

  reduction <- create_test_reduction(mask_dim)
  spec <- basis_heat_wavelet(scales = scales, order = 15, threshold = 0)

  handle <- heat_wavelet_loadings_handle(reduction, spec)
  mat <- loadings_mat(handle)

  expect_equal(nrow(mat), n_vox)
  expect_equal(ncol(mat), n_vox * length(scales))
})

# --- Data Parameter Tests ---

test_that("heat_wavelet_loadings_handle accepts data parameter", {
  reduction <- create_test_reduction()
  spec <- basis_heat_wavelet(scales = 1, order = 8, threshold = 0)

  # data is passed to lift() but often ignored for heat wavelets
  # This test verifies it doesn't cause errors
  handle <- heat_wavelet_loadings_handle(reduction, spec, data = NULL)

  expect_s4_class(handle, "LoadingsHandle")
  expect_null(handle@spec$data)
})

# --- Materialization Tests ---

test_that("loadings_mat materializes correct sparse structure", {
  reduction <- create_test_reduction()
  spec <- basis_heat_wavelet(scales = 1, order = 8, threshold = 1e-6)

  handle <- heat_wavelet_loadings_handle(reduction, spec)
  mat <- loadings_mat(handle)

  # Heat wavelet loadings should be sparse
  expect_s4_class(mat, "dgCMatrix")
})

test_that("loadings_mat with subsetting returns correct dimensions", {
  mask_dim <- c(2, 2, 2)
  n_vox <- prod(mask_dim)

  reduction <- create_test_reduction(mask_dim)
  spec <- basis_heat_wavelet(scales = c(1, 2), order = 10, threshold = 0)

  handle <- heat_wavelet_loadings_handle(reduction, spec)

  # Subset rows
  i_subset <- c(1, 3, 5)
  mat_subset <- loadings_mat(handle, i = i_subset)
  expect_equal(nrow(mat_subset), length(i_subset))

  # Subset columns
  j_subset <- c(2, 4)
  mat_subset_j <- loadings_mat(handle, j = j_subset)
  expect_equal(ncol(mat_subset_j), length(j_subset))

  # Subset both
  mat_both <- loadings_mat(handle, i = i_subset, j = j_subset)
  expect_equal(dim(mat_both), c(length(i_subset), length(j_subset)))
})

# --- Integration with LatentNeuroVec ---

test_that("LoadingsHandle from heat_wavelet integrates with LatentNeuroVec", {
  mask_dim <- c(2, 2, 2)
  n_vox <- prod(mask_dim)
  n_time <- 5L

  reduction <- create_test_reduction(mask_dim)
  spec <- basis_heat_wavelet(scales = 1, order = 8, threshold = 0)

  loadings_handle <- heat_wavelet_loadings_handle(reduction, spec)
  k <- loadings_handle@dim[2]

  # Create compatible basis
  basis <- Matrix::Matrix(matrix(rnorm(n_time * k), nrow = n_time, ncol = k))

  # Create mask and space
  mask_arr <- array(TRUE, dim = mask_dim)
  mask_vol <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(mask_dim))
  space <- neuroim2::NeuroSpace(c(mask_dim, n_time))

  # Construct LatentNeuroVec with handle
  lvec <- LatentNeuroVec(
    basis = basis,
    loadings = loadings_handle,
    space = space,
    mask = mask_vol,
    offset = numeric(0),
    label = "heat-wavelet-test"
  )

  expect_s4_class(lvec, "LatentNeuroVec")
  expect_equal(dim(lvec), c(mask_dim, n_time))
})

# --- Edge Cases ---

test_that("heat_wavelet_loadings_handle works with single scale", {
  reduction <- create_test_reduction()
  spec <- basis_heat_wavelet(scales = 2, order = 8, threshold = 0)

  handle <- heat_wavelet_loadings_handle(reduction, spec)
  mat <- loadings_mat(handle)

  n_vox <- prod(c(2, 2, 2))
  expect_equal(nrow(mat), n_vox)
  # Single scale means n_vox columns (one set of wavelets)
  expect_equal(ncol(mat), n_vox)
})

test_that("heat_wavelet_loadings_handle works with multiple scales", {
  reduction <- create_test_reduction()
  scales <- c(1, 2, 4, 8)
  spec <- basis_heat_wavelet(scales = scales, order = 15, threshold = 0)

  handle <- heat_wavelet_loadings_handle(reduction, spec)
  mat <- loadings_mat(handle)

  n_vox <- prod(c(2, 2, 2))
  expect_equal(nrow(mat), n_vox)
  expect_equal(ncol(mat), n_vox * length(scales))
})

test_that("heat_wavelet_loadings_handle respects threshold parameter", {
  reduction <- create_test_reduction()

  # Very loose threshold - more zeros
  spec_loose <- basis_heat_wavelet(scales = 1, order = 8, threshold = 0.5)
  handle_loose <- heat_wavelet_loadings_handle(reduction, spec_loose)
  mat_loose <- loadings_mat(handle_loose)


  # No threshold - no zeros introduced
  spec_none <- basis_heat_wavelet(scales = 1, order = 8, threshold = 0)
  handle_none <- heat_wavelet_loadings_handle(reduction, spec_none)
  mat_none <- loadings_mat(handle_none)

  # Loose threshold should have >= zeros than no threshold
  nnz_loose <- Matrix::nnzero(mat_loose)
  nnz_none <- Matrix::nnzero(mat_none)

  expect_true(nnz_loose <= nnz_none)
})

# --- Validity Tests ---

test_that("LoadingsHandle dim slot has correct type", {
  reduction <- create_test_reduction()
  spec <- basis_heat_wavelet(scales = 1, order = 8, threshold = 0)

  handle <- heat_wavelet_loadings_handle(reduction, spec)

  expect_type(handle@dim, "integer")
  expect_length(handle@dim, 2)
})

test_that("LoadingsHandle slots are properly initialized", {
  reduction <- create_test_reduction()
  spec <- basis_heat_wavelet(scales = 1, order = 8, threshold = 0)
  custom_id <- "validity-test-id"
  custom_label <- "validity-test-label"

  handle <- heat_wavelet_loadings_handle(reduction, spec, id = custom_id, label = custom_label)

  # Check all slots are non-empty/properly set

expect_true(nchar(handle@id) > 0)
  expect_true(length(handle@dim) == 2)
  expect_true(nchar(handle@kind) > 0)
  expect_true(is.list(handle@spec))
  expect_true(nchar(handle@label) > 0)
})
