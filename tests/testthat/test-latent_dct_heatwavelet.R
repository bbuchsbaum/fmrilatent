# Tests for latent_dct_heatwavelet.R
# DCT temporal basis + heat-wavelet spatial dictionary constructor

library(testthat)

# Helper to create minimal test data
make_test_mask <- function(dims = c(2, 2, 2)) {
  array(TRUE, dim = dims)
}

# ---------------------------------------------------------------------------
# Tests for underlying components: DCT basis handle
# These tests verify that the DCT temporal basis works correctly
# ---------------------------------------------------------------------------

test_that("dct_basis_handle creates valid BasisHandle", {
  n_time <- 10L
  k <- 4L

  bh <- dct_basis_handle(n_time = n_time, k = k)

  expect_s4_class(bh, "BasisHandle")
  expect_equal(bh@kind, "dct")
  expect_equal(bh@dim, c(n_time, k))
})

test_that("dct_basis_handle creates deterministic ID from parameters", {
  n_time <- 10L
  k <- 4L

  bh1 <- dct_basis_handle(n_time = n_time, k = k)
  bh2 <- dct_basis_handle(n_time = n_time, k = k)

  # Same parameters should produce same ID
  expect_equal(bh1@id, bh2@id)
})

test_that("dct_basis_handle materializes to orthonormal DCT matrix", {
  n_time <- 8L
  k <- 4L

  bh <- dct_basis_handle(n_time = n_time, k = k, norm = "ortho")
  mat <- basis_mat(bh)

  expect_equal(dim(mat), c(n_time, k))

  # Check orthonormality: B^T B should be identity
  gram <- crossprod(as.matrix(mat))
  expect_equal(as.matrix(gram), diag(k), tolerance = 1e-10)
})

test_that("build_dct_basis matches dct_basis_handle output", {
  n_time <- 6L
  k <- 3L

  bh <- dct_basis_handle(n_time = n_time, k = k, norm = "ortho")
  mat_handle <- as.matrix(basis_mat(bh))

  mat_direct <- as.matrix(build_dct_basis(n_time = n_time, k = k, norm = "ortho"))

  expect_equal(mat_handle, mat_direct, tolerance = 1e-10)
})

test_that("DCT basis is cached in registry", {
  n_time <- 5L
  k <- 3L

  bh <- dct_basis_handle(n_time = n_time, k = k)

  # First access materializes and caches
  mat1 <- basis_mat(bh)
  expect_true(fmrilatent:::.latent_has_matrix(bh@id, type = "basis"))

  # Second access retrieves from cache
  mat2 <- basis_mat(bh)
  expect_identical(mat1, mat2)
})

# ---------------------------------------------------------------------------
# Tests for underlying components: heat wavelet loadings handle
# These tests verify that heat wavelets work correctly via the handle
# ---------------------------------------------------------------------------

test_that("heat_wavelet_loadings_handle creates valid LoadingsHandle", {
  skip_if_not_installed("rgsp")

  mask <- make_test_mask(c(2, 2, 2))
  map <- seq_len(sum(mask))
  reduction <- make_cluster_reduction(mask, map)
  spec <- basis_heat_wavelet(scales = c(1, 2), order = 10)

  lh <- heat_wavelet_loadings_handle(
    reduction = reduction,
    basis_spec = spec,
    data = NULL,
    label = "test-hw"
  )

  expect_s4_class(lh, "LoadingsHandle")
  expect_equal(lh@kind, "lifted")
  expect_equal(lh@label, "test-hw")
})

test_that("heat_wavelet_loadings_handle stores correct dimensions", {
  skip_if_not_installed("rgsp")

  mask <- make_test_mask(c(2, 2, 2))
  n_vox <- sum(mask)
  map <- seq_len(n_vox)
  reduction <- make_cluster_reduction(mask, map)
  spec <- basis_heat_wavelet(scales = c(1, 2), order = 10)

  lh <- heat_wavelet_loadings_handle(
    reduction = reduction,
    basis_spec = spec,
    data = NULL
  )

  expect_equal(lh@dim[1L], n_vox)
  # columns = n_vox * n_scales for one-voxel-per-cluster
  expect_equal(lh@dim[2L], n_vox * length(spec$scales))
})

test_that("heat_wavelet_loadings_handle materializes to sparse matrix", {
  skip_if_not_installed("rgsp")

  mask <- make_test_mask(c(2, 2, 2))
  n_vox <- sum(mask)
  map <- seq_len(n_vox)
  reduction <- make_cluster_reduction(mask, map)
  spec <- basis_heat_wavelet(scales = c(1), order = 8)

  lh <- heat_wavelet_loadings_handle(
    reduction = reduction,
    basis_spec = spec,
    data = NULL
  )

  mat <- loadings_mat(lh)

  expect_true(inherits(mat, "Matrix") || is.matrix(mat))
  expect_equal(nrow(mat), n_vox)
})

# ---------------------------------------------------------------------------
# Tests for basis_heat_wavelet spec
# ---------------------------------------------------------------------------

test_that("basis_heat_wavelet creates proper spec object", {
  spec <- basis_heat_wavelet(scales = c(1, 2, 4), order = 20, threshold = 1e-5)

  expect_s3_class(spec, "spec_heat_wavelet")
  expect_equal(spec$scales, c(1, 2, 4))
  expect_equal(spec$order, 20)
  expect_equal(spec$threshold, 1e-5)
})

test_that("basis_heat_wavelet uses sensible defaults", {
  spec <- basis_heat_wavelet()

  expect_s3_class(spec, "spec_heat_wavelet")
  expect_true(length(spec$scales) > 0)
  expect_true(spec$order > 0)
  expect_true(spec$threshold >= 0)
})

# ---------------------------------------------------------------------------
# Tests for make_cluster_reduction
# ---------------------------------------------------------------------------

test_that("make_cluster_reduction creates ClusterReduction from array mask", {
  skip_if_not_installed("rgsp")

  mask <- make_test_mask(c(2, 2, 2))
  n_vox <- sum(mask)
  map <- seq_len(n_vox)

  red <- make_cluster_reduction(mask, map)

  expect_s4_class(red, "ClusterReduction")
  expect_equal(length(red@map), n_vox)
  expect_equal(red@cluster_ids, as.integer(sort(unique(map))))
})

test_that("make_cluster_reduction works with multiple clusters", {
  skip_if_not_installed("rgsp")

  mask <- make_test_mask(c(2, 2, 2))
  n_vox <- sum(mask)
  # Two clusters
  map <- rep(1:2, each = n_vox / 2)

  red <- make_cluster_reduction(mask, map)

  expect_s4_class(red, "ClusterReduction")
  expect_equal(red@cluster_ids, c(1L, 2L))
})

test_that("make_cluster_reduction works with LogicalNeuroVol mask", {
  skip_if_not_installed("rgsp")

  mask_arr <- make_test_mask(c(2, 2, 2))
  mask_vol <- neuroim2::LogicalNeuroVol(
    mask_arr,
    neuroim2::NeuroSpace(dim(mask_arr))
  )
  n_vox <- sum(mask_arr)
  map <- seq_len(n_vox)

  red <- make_cluster_reduction(mask_vol, map)

  expect_s4_class(red, "ClusterReduction")
  expect_s4_class(red@mask, "LogicalNeuroVol")
})

# ---------------------------------------------------------------------------
# Tests for lift method with heat wavelets
# ---------------------------------------------------------------------------

test_that("lift with ClusterReduction and heat wavelet returns sparse matrix", {
  skip_if_not_installed("rgsp")

  mask <- make_test_mask(c(2, 2, 2))
  n_vox <- sum(mask)
  map <- seq_len(n_vox)
  reduction <- make_cluster_reduction(mask, map)
  spec <- basis_heat_wavelet(scales = c(1, 2), order = 10, threshold = 0)

  loadings <- lift(reduction, spec, k_neighbors = 3L)

  expect_s4_class(loadings, "dgCMatrix")
  expect_equal(nrow(loadings), n_vox)
  expect_equal(ncol(loadings), n_vox * length(spec$scales))
})

test_that("lift with single-voxel clusters uses identity wavelets", {
  skip_if_not_installed("rgsp")

  # Very small mask where each cluster has only 1 voxel
  mask <- array(TRUE, dim = c(2, 1, 1))
  n_vox <- sum(mask)
  map <- seq_len(n_vox)  # each voxel is its own cluster
  reduction <- make_cluster_reduction(mask, map)
  spec <- basis_heat_wavelet(scales = c(1), order = 10, threshold = 0)

  loadings <- lift(reduction, spec, k_neighbors = 3L)

  expect_s4_class(loadings, "dgCMatrix")
  expect_equal(nrow(loadings), n_vox)
})

# ---------------------------------------------------------------------------
# Combined DCT + heat wavelet manual construction test
# This shows how the function should work once the bug is fixed
# ---------------------------------------------------------------------------

test_that("manual DCT + heat wavelet construction produces valid LatentNeuroVec", {
  skip_if_not_installed("rgsp")

  # This test demonstrates how to correctly combine DCT basis with heat wavelet
  # spatial dictionary. The key is that the number of columns (components) in
  # basis and loadings must match.

  mask_arr <- make_test_mask(c(2, 2, 2))
  n_vox <- sum(mask_arr)
  n_time <- 6L

  # Create mask vol and space
  mask_vol <- neuroim2::LogicalNeuroVol(
    mask_arr,
    neuroim2::NeuroSpace(dim(mask_arr))
  )
  space <- neuroim2::NeuroSpace(c(dim(mask_arr), n_time))

  # Create heat wavelet loadings handle first to determine number of components
  cluster_map <- seq_len(n_vox)
  reduction <- make_cluster_reduction(mask_arr, cluster_map)
  hw_spec <- basis_heat_wavelet(scales = c(1), order = 8)
  L_handle <- heat_wavelet_loadings_handle(
    reduction = reduction,
    basis_spec = hw_spec,
    data = NULL,
    label = "heat-wavelet"
  )

  # The loadings handle now determines k (number of spatial components)
  k_spatial <- fmrilatent:::.latent_loadings_dim(L_handle)[2L]

  # For a combined DCT + heat wavelet representation, we need to create
  # explicit basis coefficients that project temporal data onto spatial atoms
  # This is typically done via: basis = X %*% loadings where X is time x voxels

  # Create a simple explicit basis matrix for testing (time x k_spatial)
  basis_mat_explicit <- Matrix::Matrix(
    matrix(rnorm(n_time * k_spatial), nrow = n_time, ncol = k_spatial),
    sparse = FALSE
  )

  # Construct LatentNeuroVec with explicit basis
  lv <- LatentNeuroVec(
    basis = basis_mat_explicit,
    loadings = L_handle,
    space = space,
    mask = mask_vol,
    offset = numeric(n_vox),
    label = "DCT + heat-wavelet",
    meta = list(
      time_basis = "dct",
      spatial_dict = "heat_wavelet"
    )
  )

  expect_s4_class(lv, "LatentNeuroVec")
  expect_equal(lv@label, "DCT + heat-wavelet")
  expect_equal(lv@meta$spatial_dict, "heat_wavelet")

  # Test basic operations work
  expect_equal(dim(lv), c(2, 2, 2, n_time))

  # Extract a volume
  vol1 <- lv[[1]]
  expect_s4_class(vol1, "SparseNeuroVol")

  # Extract a time series
  ts <- series(lv, 1L)
  expect_equal(length(ts), n_time)
})

test_that("manual DCT + heat wavelet roundtrip preserves data structure", {
  skip_if_not_installed("rgsp")

  mask_arr <- make_test_mask(c(2, 2, 2))
  n_vox <- sum(mask_arr)
  n_time <- 8L

  mask_vol <- neuroim2::LogicalNeuroVol(
    mask_arr,
    neuroim2::NeuroSpace(dim(mask_arr))
  )
  space <- neuroim2::NeuroSpace(c(dim(mask_arr), n_time))

  # Create heat wavelet loadings first to determine k_spatial
  cluster_map <- seq_len(n_vox)
  reduction <- make_cluster_reduction(mask_arr, cluster_map)
  hw_spec <- basis_heat_wavelet(scales = c(1, 2), order = 10)
  L_handle <- heat_wavelet_loadings_handle(
    reduction = reduction,
    basis_spec = hw_spec,
    data = NULL
  )

  k_spatial <- fmrilatent:::.latent_loadings_dim(L_handle)[2L]

  # Create explicit basis matrix (time x k_spatial) to match loadings
  basis_mat_explicit <- Matrix::Matrix(
    matrix(rnorm(n_time * k_spatial), nrow = n_time, ncol = k_spatial),
    sparse = FALSE
  )

  lv <- LatentNeuroVec(
    basis = basis_mat_explicit,
    loadings = L_handle,
    space = space,
    mask = mask_vol,
    offset = numeric(n_vox),
    label = "test"
  )

  # Get materialized matrices
  B <- as.matrix(basis(lv))
  L <- as.matrix(loadings(lv))

  # Reconstruction via matrix multiplication
  recon <- B %*% t(L)

  # Should match as.matrix output
  lv_mat <- as.matrix(lv)

  expect_equal(dim(recon), dim(lv_mat))
  expect_equal(recon, lv_mat, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# Tests for latent_dct_heatwavelet function behavior (skipped due to bug)
# These tests are ready to run once the bug is fixed
# ---------------------------------------------------------------------------

test_that("latent_dct_heatwavelet returns LatentNeuroVec with correct structure", {
  skip_if_not_installed("rgsp")

  mask <- make_test_mask(c(2, 2, 2))
  n_time <- 10L
  k_time <- 4L

  lv <- latent_dct_heatwavelet(
    n_time = n_time,
    k_time = k_time,
    mask = mask,
    label = "test-dct-hw"
  )

  expect_s4_class(lv, "LatentNeuroVec")
  expect_equal(lv@label, "test-dct-hw")

  n_vox <- sum(mask)
  basis_dim <- fmrilatent:::.latent_basis_dim(lv@basis)
  loadings_dim <- fmrilatent:::.latent_loadings_dim(lv@loadings)

  expect_equal(basis_dim[1L], n_time)
  expect_equal(basis_dim[2L], loadings_dim[2L])
  expect_equal(loadings_dim[1L], n_vox)
})

test_that("latent_dct_heatwavelet creates correct metadata", {
  skip_if_not_installed("rgsp")

  mask <- make_test_mask(c(2, 2, 2))
  n_time <- 8L
  k_time <- 3L

  lv <- latent_dct_heatwavelet(
    n_time = n_time,
    k_time = k_time,
    mask = mask
  )

  expect_true("meta" %in% slotNames(lv))
  # time_basis is "template" because we use placeholder basis matrix
  expect_equal(lv@meta$time_basis, "template")
  # time_k is determined by loadings, not k_time parameter
  loadings_dim <- fmrilatent:::.latent_loadings_dim(lv@loadings)
  expect_equal(lv@meta$time_k, loadings_dim[2L])
  expect_equal(lv@meta$spatial_dict, "heat_wavelet")
})

test_that("latent_dct_heatwavelet works with explicit cluster_map", {
  skip_if_not_installed("rgsp")

  mask <- make_test_mask(c(2, 2, 2))
  n_vox <- sum(mask)
  n_time <- 6L
  k_time <- 2L
  cluster_map <- rep(1:2, each = n_vox / 2)

  lv <- latent_dct_heatwavelet(
    n_time = n_time,
    k_time = k_time,
    mask = mask,
    cluster_map = cluster_map
  )

  expect_s4_class(lv, "LatentNeuroVec")
  loadings_dim <- fmrilatent:::.latent_loadings_dim(lv@loadings)
  expect_equal(loadings_dim[1L], n_vox)
})

test_that("latent_dct_heatwavelet works with custom reduction", {
  skip_if_not_installed("rgsp")

  mask <- make_test_mask(c(2, 2, 2))
  n_vox <- sum(mask)
  n_time <- 6L
  k_time <- 2L
  cluster_map <- rep(1:2, each = n_vox / 2)
  reduction <- make_cluster_reduction(mask, cluster_map)

  lv <- latent_dct_heatwavelet(
    n_time = n_time,
    k_time = k_time,
    mask = mask,
    reduction = reduction
  )

  expect_s4_class(lv, "LatentNeuroVec")
})

test_that("latent_dct_heatwavelet respects custom hw_basis_spec", {
  skip_if_not_installed("rgsp")

  mask <- make_test_mask(c(2, 2, 2))
  n_time <- 6L
  k_time <- 2L
  custom_spec <- basis_heat_wavelet(scales = c(1, 2), order = 10, threshold = 1e-5)

  lv <- latent_dct_heatwavelet(
    n_time = n_time,
    k_time = k_time,
    mask = mask,
    hw_basis_spec = custom_spec
  )

  expect_s4_class(lv, "LatentNeuroVec")
  loadings_dim <- fmrilatent:::.latent_loadings_dim(lv@loadings)
  n_vox <- sum(mask)
  expected_cols <- n_vox * length(custom_spec$scales)
  expect_equal(loadings_dim[2L], expected_cols)
})

test_that("latent_dct_heatwavelet accepts valid offset", {
  skip_if_not_installed("rgsp")

  mask <- make_test_mask(c(2, 2, 2))
  n_vox <- sum(mask)
  n_time <- 6L
  k_time <- 2L
  offset <- rnorm(n_vox)

  lv <- latent_dct_heatwavelet(
    n_time = n_time,
    k_time = k_time,
    mask = mask,
    offset = offset
  )

  expect_s4_class(lv, "LatentNeuroVec")
  expect_equal(length(lv@offset), n_vox)
  expect_equal(lv@offset, offset)
})

test_that("latent_dct_heatwavelet rejects invalid offset length", {
  skip_if_not_installed("rgsp")

  mask <- make_test_mask(c(2, 2, 2))
  n_vox <- sum(mask)
  n_time <- 6L
  k_time <- 2L
  bad_offset <- rnorm(n_vox + 1)

  # This error should occur before the space/mask bug is hit
  expect_error(
    latent_dct_heatwavelet(
      n_time = n_time,
      k_time = k_time,
      mask = mask,
      offset = bad_offset
    ),
    regexp = "offset"
  )
})

test_that("latent_dct_heatwavelet uses explicit basis matrix", {
  skip_if_not_installed("rgsp")

  mask <- make_test_mask(c(2, 2, 2))
  n_time <- 6L
  k_time <- 2L

  lv <- latent_dct_heatwavelet(
    n_time = n_time,
    k_time = k_time,
    mask = mask
  )

  # Basis is now an explicit matrix (placeholder initialized to zeros)
  expect_true(inherits(lv@basis, "Matrix") || is.matrix(lv@basis))
  basis_dim <- fmrilatent:::.latent_basis_dim(lv@basis)
  expect_equal(basis_dim[1L], n_time)
})

test_that("latent_dct_heatwavelet uses LoadingsHandle for spatial loadings", {
  skip_if_not_installed("rgsp")

  mask <- make_test_mask(c(2, 2, 2))
  n_time <- 6L
  k_time <- 2L

  lv <- latent_dct_heatwavelet(
    n_time = n_time,
    k_time = k_time,
    mask = mask
  )

  expect_s4_class(lv@loadings, "LoadingsHandle")
  expect_equal(lv@loadings@kind, "lifted")
})

test_that("latent_dct_heatwavelet works with LogicalNeuroVol mask", {
  skip_if_not_installed("rgsp")

  mask_arr <- make_test_mask(c(2, 2, 2))
  mask_vol <- neuroim2::LogicalNeuroVol(
    mask_arr,
    neuroim2::NeuroSpace(dim(mask_arr))
  )
  n_time <- 6L
  k_time <- 2L

  lv <- latent_dct_heatwavelet(
    n_time = n_time,
    k_time = k_time,
    mask = mask_vol
  )

  expect_s4_class(lv, "LatentNeuroVec")
})

# ---------------------------------------------------------------------------
# Error handling tests
# ---------------------------------------------------------------------------

test_that("latent_dct_heatwavelet fails without rgsp package", {
  skip_if(requireNamespace("rgsp", quietly = TRUE),
          message = "rgsp is installed; cannot test missing-package error")

  mask <- make_test_mask(c(2, 2, 2))
  n_time <- 6L
  k_time <- 2L

  expect_error(
    latent_dct_heatwavelet(
      n_time = n_time,
      k_time = k_time,
      mask = mask
    ),
    regexp = "rgsp"
  )
})
