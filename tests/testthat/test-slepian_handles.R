# Tests for R/slepian_handles.R
# Covers slepian_temporal_handle and slepian_spatial_loadings_handle

# =============================================================================
# slepian_temporal_handle tests
# =============================================================================

test_that("slepian_temporal_handle creates BasisHandle with correct structure", {
  n_time <- 20L
  tr <- 2.0
  bandwidth <- 0.1

  handle <- slepian_temporal_handle(n_time = n_time, tr = tr, bandwidth = bandwidth)

  expect_s4_class(handle, "BasisHandle")
  expect_equal(handle@kind, "slepian_temporal")
  expect_equal(handle@dim[1], n_time)
  expect_true(handle@dim[2] > 0)
  expect_true(nchar(handle@id) > 0)
  expect_true(nchar(handle@label) > 0)
})

test_that("slepian_temporal_handle computes k from bandwidth when not provided", {
  n_time <- 50L
  tr <- 1.5
  bandwidth <- 0.05

  # Expected k = floor(2 * NW) - 1 where NW = n_time * bandwidth * tr
  NW <- n_time * bandwidth * tr
  expected_k <- floor(2 * NW) - 1L

  handle <- slepian_temporal_handle(n_time = n_time, tr = tr, bandwidth = bandwidth)

  expect_equal(handle@dim[2], expected_k)
  expect_equal(handle@spec$k, expected_k)
})

test_that("slepian_temporal_handle uses explicit k when provided", {
  n_time <- 30L
  tr <- 2.0
  bandwidth <- 0.1
  k <- 5L

  handle <- slepian_temporal_handle(n_time = n_time, tr = tr, bandwidth = bandwidth, k = k)

  expect_equal(handle@dim[2], k)
  expect_equal(handle@spec$k, k)
})

test_that("slepian_temporal_handle stores all spec parameters correctly", {
  n_time <- 25L
  tr <- 1.0
  bandwidth <- 0.08
  k <- 4L
  backend <- "dense"

  handle <- slepian_temporal_handle(
    n_time = n_time,
    tr = tr,
    bandwidth = bandwidth,
    k = k,
    backend = backend
  )

  expect_equal(handle@spec$n_time, n_time)
  expect_equal(handle@spec$tr, tr)
  expect_equal(handle@spec$bandwidth, bandwidth)
  expect_equal(handle@spec$k, k)
  expect_equal(handle@spec$backend, backend)
})

test_that("slepian_temporal_handle uses custom id when provided", {
  custom_id <- "my-custom-slepian-id"
  handle <- slepian_temporal_handle(n_time = 10L, tr = 2.0, bandwidth = 0.1, id = custom_id)

  expect_equal(handle@id, custom_id)
})

test_that("slepian_temporal_handle uses custom label when provided", {
  custom_label <- "My Custom Slepian Basis"
  handle <- slepian_temporal_handle(n_time = 10L, tr = 2.0, bandwidth = 0.1, label = custom_label)

  expect_equal(handle@label, custom_label)
})

test_that("slepian_temporal_handle generates unique ids by default", {
  handle1 <- slepian_temporal_handle(n_time = 10L, tr = 2.0, bandwidth = 0.1)
  handle2 <- slepian_temporal_handle(n_time = 10L, tr = 2.0, bandwidth = 0.2)
  handle3 <- slepian_temporal_handle(n_time = 20L, tr = 2.0, bandwidth = 0.1)

  # Different bandwidths should give different ids
  expect_false(handle1@id == handle2@id)
  # Different n_time should give different ids
  expect_false(handle1@id == handle3@id)

  # Same parameters should give same auto-generated id
  handle4 <- slepian_temporal_handle(n_time = 10L, tr = 2.0, bandwidth = 0.1)
  expect_equal(handle1@id, handle4@id)
})

test_that("slepian_temporal_handle auto id changes with tr", {
  handle1 <- slepian_temporal_handle(n_time = 20L, tr = 1.0, bandwidth = 0.1, k = 3L)
  handle2 <- slepian_temporal_handle(n_time = 20L, tr = 2.0, bandwidth = 0.1, k = 3L)

  expect_false(handle1@id == handle2@id)
})

test_that("slepian_temporal_handle defaults to tridiag backend", {
  handle <- slepian_temporal_handle(n_time = 15L, tr = 2.0, bandwidth = 0.1)

  expect_equal(handle@spec$backend, "tridiag")
})

test_that("slepian_temporal_handle accepts dense backend", {
  handle <- slepian_temporal_handle(n_time = 10L, tr = 2.0, bandwidth = 0.1, backend = "dense")

  expect_equal(handle@spec$backend, "dense")
})

test_that("slepian_temporal_handle converts n_time to integer", {
  # Pass numeric, should be converted to integer
  handle <- slepian_temporal_handle(n_time = 15.0, tr = 2.0, bandwidth = 0.1)

  expect_true(is.integer(handle@dim[1]))
  expect_equal(handle@dim[1], 15L)
  expect_true(is.integer(handle@spec$n_time))
})

test_that("slepian_temporal_handle converts k to integer", {
  handle <- slepian_temporal_handle(n_time = 20L, tr = 2.0, bandwidth = 0.1, k = 3.0)

  expect_true(is.integer(handle@spec$k))
  expect_equal(handle@spec$k, 3L)
})

test_that("slepian_temporal_handle label contains parameters", {
  n_time <- 30L
  bandwidth <- 0.05

  handle <- slepian_temporal_handle(n_time = n_time, tr = 2.0, bandwidth = bandwidth)

  # Auto-generated label should contain n_time, k, and bandwidth info
  expect_true(grepl("30", handle@label))
  expect_true(grepl(as.character(handle@spec$k), handle@label))
  expect_true(grepl("0.05", handle@label))
})

test_that("slepian_temporal_handle works with small bandwidth", {
  # Small bandwidth means fewer components
  handle <- slepian_temporal_handle(n_time = 100L, tr = 2.0, bandwidth = 0.01)

  expect_s4_class(handle, "BasisHandle")
  expect_true(handle@dim[2] >= 1)  # At least one component
})

test_that("slepian_temporal_handle works with large bandwidth", {
  handle <- slepian_temporal_handle(n_time = 50L, tr = 2.0, bandwidth = 0.2)

  expect_s4_class(handle, "BasisHandle")
  expect_true(handle@dim[2] > 0)
})

# =============================================================================
# basis_mat materialization for slepian_temporal BasisHandle
# =============================================================================

test_that("basis_mat materializes slepian_temporal handle correctly", {
  n_time <- 12L
  tr <- 2.0
  bandwidth <- 0.1

  handle <- slepian_temporal_handle(n_time = n_time, tr = tr, bandwidth = bandwidth)
  mat <- basis_mat(handle)

  expect_true(is.matrix(mat) || inherits(mat, "Matrix"))
  expect_equal(nrow(mat), n_time)
  expect_equal(ncol(mat), handle@dim[2])
})

test_that("basis_mat for slepian handle produces orthonormal columns", {
  handle <- slepian_temporal_handle(n_time = 20L, tr = 2.0, bandwidth = 0.1)
  mat <- as.matrix(basis_mat(handle))

  # Check orthonormality
  crossprod_mat <- crossprod(mat)
  expect_equal(crossprod_mat, diag(ncol(mat)), tolerance = 1e-10)
})

test_that("basis_mat with subsetting returns correct dimensions", {
  handle <- slepian_temporal_handle(n_time = 20L, tr = 2.0, bandwidth = 0.1, k = 5L)

  # Subset rows
  mat_sub_rows <- basis_mat(handle, i = 1:10)
  expect_equal(nrow(mat_sub_rows), 10)
  expect_equal(ncol(mat_sub_rows), 5)

  # Subset columns
  mat_sub_cols <- basis_mat(handle, j = 1:3)
  expect_equal(nrow(mat_sub_cols), 20)
  expect_equal(ncol(mat_sub_cols), 3)

  # Subset both
  mat_sub_both <- basis_mat(handle, i = 5:15, j = 2:4)
  expect_equal(nrow(mat_sub_both), 11)
  expect_equal(ncol(mat_sub_both), 3)
})

test_that("basis_mat caches materialized slepian handle in registry", {
  # Create a fresh handle with unique id to avoid conflicts
  unique_id <- paste0("test-slepian-cache-", format(Sys.time(), "%Y%m%d%H%M%S"), "-", sample(1000, 1))
  handle <- slepian_temporal_handle(n_time = 15L, tr = 2.0, bandwidth = 0.1, id = unique_id)

  # First call should materialize
  mat1 <- basis_mat(handle)

  # Should be in registry now
  expect_true(fmrilatent:::`.latent_has_matrix`(handle@id, type = "basis"))

  # Second call should return same object
  mat2 <- basis_mat(handle)
  expect_identical(mat1, mat2)
})

test_that("slepian_temporal_handle with tridiag backend matches dpss_time_basis", {
  n_time <- 15L
  tr <- 2.0
  bandwidth <- 0.08
  k <- 4L

  handle <- slepian_temporal_handle(
    n_time = n_time,
    tr = tr,
    bandwidth = bandwidth,
    k = k,
    backend = "tridiag"
  )

  mat_handle <- as.matrix(basis_mat(handle))
  mat_direct <- dpss_time_basis(n_time, tr = tr, bandwidth = bandwidth, k = k, backend = "tridiag")

  expect_equal(mat_handle, mat_direct, tolerance = 1e-10)
})

test_that("slepian_temporal_handle with dense backend matches dpss_time_basis", {
  n_time <- 12L
  tr <- 1.5
  bandwidth <- 0.1
  k <- 3L

  # Use unique id to avoid cache collision with tridiag backend test
  # Note: The auto-generated id does not include backend, which can cause

  # cache collisions if the same parameters are used with different backends.
  unique_id <- paste0("test-dense-", format(Sys.time(), "%H%M%S"), "-", sample(10000, 1))

  handle <- slepian_temporal_handle(
    n_time = n_time,
    tr = tr,
    bandwidth = bandwidth,
    k = k,
    backend = "dense",
    id = unique_id
  )

  mat_handle <- as.matrix(basis_mat(handle))
  mat_direct <- dpss_time_basis(n_time, tr = tr, bandwidth = bandwidth, k = k, backend = "dense")

  expect_equal(mat_handle, mat_direct, tolerance = 1e-10)
})

# =============================================================================
# slepian_spatial_loadings_handle tests
# =============================================================================

test_that("slepian_spatial_loadings_handle creates LoadingsHandle with correct structure", {
  skip_if_not_installed("RSpectra")

  # Create a simple cluster reduction
  mask <- array(TRUE, dim = c(3, 3, 2))
  map <- seq_len(sum(mask))
  red <- make_cluster_reduction(mask, map)
  spec <- basis_slepian(k = 2)

  handle <- slepian_spatial_loadings_handle(
    reduction = red,
    basis_spec = spec,
    label = "test-spatial"
  )

  expect_s4_class(handle, "LoadingsHandle")
  expect_equal(handle@kind, "slepian_spatial")
  expect_equal(handle@label, "test-spatial")
  expect_true(nchar(handle@id) > 0)
})

test_that("slepian_spatial_loadings_handle reuses deterministic id when inputs match", {
  skip_if_not_installed("RSpectra")

  mask <- array(TRUE, dim = c(2, 2, 2))
  map <- seq_len(sum(mask))
  red <- make_cluster_reduction(mask, map)
  spec <- basis_slepian(k = 1)

  handle1 <- slepian_spatial_loadings_handle(reduction = red, basis_spec = spec)
  handle2 <- slepian_spatial_loadings_handle(reduction = red, basis_spec = spec)

  expect_equal(handle1@id, handle2@id)
})

test_that("slepian_spatial_loadings_handle uses custom id when provided", {
  skip_if_not_installed("RSpectra")

  mask <- array(TRUE, dim = c(2, 2, 2))
  map <- seq_len(sum(mask))
  red <- make_cluster_reduction(mask, map)
  spec <- basis_slepian(k = 1)
  custom_id <- "my-spatial-slepian-id"

  handle <- slepian_spatial_loadings_handle(
    reduction = red,
    basis_spec = spec,
    id = custom_id
  )

  expect_equal(handle@id, custom_id)
})

test_that("slepian_spatial_loadings_handle stores spec parameters correctly", {
  skip_if_not_installed("RSpectra")

  mask <- array(TRUE, dim = c(2, 2, 2))
  map <- seq_len(sum(mask))
  red <- make_cluster_reduction(mask, map)
  spec <- basis_slepian(k = 2)

  handle <- slepian_spatial_loadings_handle(reduction = red, basis_spec = spec)

  expect_equal(handle@spec$family, "slepian_spatial")
  expect_identical(handle@spec$basis_spec, spec)
  expect_s4_class(handle@spec$reduction, "ClusterReduction")
})

test_that("slepian_spatial_loadings_handle dim matches lifted matrix", {
  skip_if_not_installed("RSpectra")

  mask <- array(TRUE, dim = c(3, 3, 2))
  n_vox <- sum(mask)
  map <- seq_len(n_vox)
  red <- make_cluster_reduction(mask, map)
  spec <- basis_slepian(k = 3)

  handle <- slepian_spatial_loadings_handle(reduction = red, basis_spec = spec)

  # dim should be (n_voxels, k_components)
  expect_equal(handle@dim[1], n_vox)
  expect_true(handle@dim[2] > 0)
})

test_that("slepian_spatial_loadings_handle registers matrix in loadings registry", {
  skip_if_not_installed("RSpectra")

  mask <- array(TRUE, dim = c(2, 2, 2))
  map <- seq_len(sum(mask))
  red <- make_cluster_reduction(mask, map)
  spec <- basis_slepian(k = 2)
  unique_id <- paste0("test-spatial-reg-", format(Sys.time(), "%Y%m%d%H%M%S"), "-", sample(1000, 1))

  handle <- slepian_spatial_loadings_handle(
    reduction = red,
    basis_spec = spec,
    id = unique_id
  )

  # Matrix should be registered during construction
  expect_true(fmrilatent:::`.latent_has_matrix`(handle@id, type = "loadings"))
})

test_that("loadings_mat retrieves correct matrix for slepian_spatial handle", {
  skip_if_not_installed("RSpectra")

  mask <- array(TRUE, dim = c(2, 2, 2))
  map <- seq_len(sum(mask))
  red <- make_cluster_reduction(mask, map)
  spec <- basis_slepian(k = 2)

  handle <- slepian_spatial_loadings_handle(reduction = red, basis_spec = spec)

  mat <- loadings_mat(handle)

  expect_true(inherits(mat, "Matrix") || is.matrix(mat))
  expect_equal(nrow(mat), handle@dim[1])
  expect_equal(ncol(mat), handle@dim[2])
})

test_that("loadings_mat with subsetting returns correct dimensions", {
  skip_if_not_installed("RSpectra")

  mask <- array(TRUE, dim = c(3, 3, 2))
  n_vox <- sum(mask)
  map <- seq_len(n_vox)
  red <- make_cluster_reduction(mask, map)
  spec <- basis_slepian(k = 4)

  handle <- slepian_spatial_loadings_handle(reduction = red, basis_spec = spec)

  # Subset rows (voxels)
  mat_sub_rows <- loadings_mat(handle, i = 1:5)
  expect_equal(nrow(mat_sub_rows), 5)

  # Subset columns (components)
  mat_sub_cols <- loadings_mat(handle, j = 1:2)
  expect_equal(ncol(mat_sub_cols), 2)

  # Subset both
  mat_sub_both <- loadings_mat(handle, i = 1:8, j = 2:3)
  expect_equal(nrow(mat_sub_both), 8)
  expect_equal(ncol(mat_sub_both), 2)
})

# =============================================================================
# Edge cases and integration tests
# =============================================================================

test_that("slepian_temporal_handle works with minimum viable n_time", {
  # Very small n_time
  handle <- slepian_temporal_handle(n_time = 3L, tr = 2.0, bandwidth = 0.1, k = 1L)

  expect_s4_class(handle, "BasisHandle")
  expect_equal(handle@dim[1], 3L)
  expect_equal(handle@dim[2], 1L)

  mat <- basis_mat(handle)
  expect_equal(nrow(mat), 3)
  expect_equal(ncol(mat), 1)
})

test_that("slepian_temporal_handle works with various tr values", {
  for (tr in c(0.5, 1.0, 2.0, 3.0)) {
    handle <- slepian_temporal_handle(n_time = 20L, tr = tr, bandwidth = 0.1)
    expect_s4_class(handle, "BasisHandle")
    expect_true(handle@dim[2] >= 1)
  }
})

test_that("slepian handles integrate with LatentNeuroVec", {
  # Test using slepian_temporal_handle with LatentNeuroVec
  mask_arr <- array(TRUE, dim = c(2, 2, 1))
  mask_vol <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(2, 2, 1)))
  n_time <- 10L
  k <- 3L

  # Create slepian temporal handle
  b_handle <- slepian_temporal_handle(n_time = n_time, tr = 2.0, bandwidth = 0.1, k = k)

  # Create explicit loadings
  load_mat <- Matrix::Matrix(
    matrix(rnorm(sum(mask_arr) * k), nrow = sum(mask_arr), ncol = k),
    sparse = FALSE
  )

  space <- neuroim2::NeuroSpace(c(2, 2, 1, n_time))

  lvec <- LatentNeuroVec(
    basis = b_handle,
    loadings = load_mat,
    space = space,
    mask = mask_vol,
    offset = numeric(0),
    label = "slepian-test"
  )

  expect_s4_class(lvec, "LatentNeuroVec")
  expect_equal(dim(lvec)[4], n_time)

  # Verify we can materialize and compute
  arr <- as.array(lvec)
  expect_equal(dim(arr), c(2, 2, 1, n_time))
})

test_that("multiple slepian handles with same parameters share cached basis", {
  # Two handles with same auto-generated id should share the cached matrix
  n_time <- 15L
  tr <- 2.0
  bandwidth <- 0.1

  handle1 <- slepian_temporal_handle(n_time = n_time, tr = tr, bandwidth = bandwidth)
  mat1 <- basis_mat(handle1)

  handle2 <- slepian_temporal_handle(n_time = n_time, tr = tr, bandwidth = bandwidth)
  mat2 <- basis_mat(handle2)

  # Same id means same cache entry
  expect_equal(handle1@id, handle2@id)
  expect_identical(mat1, mat2)
})

test_that("basis_slepian creates correct spec structure", {
  spec <- basis_slepian(k = 5, type = "laplacian")

  expect_equal(class(spec), "spec_slepian")
  expect_equal(spec$k, 5)
  expect_equal(spec$type, "laplacian")
})

test_that("basis_slepian uses default parameters", {
  spec <- basis_slepian()

  expect_equal(spec$k, 3)
  expect_equal(spec$type, "laplacian")
})

# =============================================================================
# Dimension helper tests for handles
# =============================================================================

test_that(".latent_basis_dim works with BasisHandle", {
  handle <- slepian_temporal_handle(n_time = 20L, tr = 2.0, bandwidth = 0.1, k = 5L)

  dims <- fmrilatent:::`.latent_basis_dim`(handle)

  expect_equal(dims, c(20L, 5L))
})

test_that(".latent_loadings_dim works with LoadingsHandle", {
  skip_if_not_installed("RSpectra")

  mask <- array(TRUE, dim = c(3, 3, 2))
  map <- seq_len(sum(mask))
  red <- make_cluster_reduction(mask, map)
  spec <- basis_slepian(k = 3)

  handle <- slepian_spatial_loadings_handle(reduction = red, basis_spec = spec)

  dims <- fmrilatent:::`.latent_loadings_dim`(handle)

  expect_equal(dims[1], sum(mask))
  expect_equal(dims[2], handle@dim[2])
})

# =============================================================================
# Edge case: auto-generated ids reflect backend-specific cache keys
# =============================================================================

test_that("handles with same params but different backends get different auto-generated ids", {
  h1 <- slepian_temporal_handle(n_time = 10L, tr = 2.0, bandwidth = 0.1, backend = "tridiag")
  h2 <- slepian_temporal_handle(n_time = 10L, tr = 2.0, bandwidth = 0.1, backend = "dense")

  expect_false(h1@id == h2@id)

  # Different spec backends
  expect_equal(h1@spec$backend, "tridiag")
  expect_equal(h2@spec$backend, "dense")
})

test_that("explicit unique ids prevent cache collisions between backends", {
  # When using unique explicit ids, each handle gets its own cache entry
  unique_id1 <- paste0("backend-test-tridiag-", sample(100000, 1))
  unique_id2 <- paste0("backend-test-dense-", sample(100000, 1))

  h1 <- slepian_temporal_handle(n_time = 10L, tr = 2.0, bandwidth = 0.1, k = 2L,
                                 backend = "tridiag", id = unique_id1)
  h2 <- slepian_temporal_handle(n_time = 10L, tr = 2.0, bandwidth = 0.1, k = 2L,
                                 backend = "dense", id = unique_id2)

  mat1 <- as.matrix(basis_mat(h1))
  mat2 <- as.matrix(basis_mat(h2))

  # Different backends may produce numerically different results
  # (though both are valid orthonormal bases)
  expect_equal(dim(mat1), dim(mat2))
  expect_equal(dim(mat1), c(10L, 2L))

  # Each should match its respective direct call
  mat1_direct <- dpss_time_basis(10L, tr = 2.0, bandwidth = 0.1, k = 2L, backend = "tridiag")
  mat2_direct <- dpss_time_basis(10L, tr = 2.0, bandwidth = 0.1, k = 2L, backend = "dense")

  expect_equal(mat1, mat1_direct, tolerance = 1e-10)
  expect_equal(mat2, mat2_direct, tolerance = 1e-10)
})
