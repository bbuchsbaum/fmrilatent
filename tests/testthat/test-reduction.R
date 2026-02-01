library(testthat)

# =============================================================================
# Tests for R/reduction.R - Graph reduction scaffolds and basis specifications
# =============================================================================

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

# Helper to create a minimal LogicalNeuroVol for testing
make_test_mask_vol <- function(dims = c(3, 3, 3), all_true = TRUE) {
  if (all_true) {
    mask_arr <- array(TRUE, dim = dims)
  } else {
    # Create a partial mask with some FALSE values
    mask_arr <- array(FALSE, dim = dims)
    mask_arr[1:2, 1:2, 1:2] <- TRUE
  }
  spc <- neuroim2::NeuroSpace(dims)
  neuroim2::LogicalNeuroVol(mask_arr, spc)
}

# -----------------------------------------------------------------------------
# Tests for basis specification constructors
# -----------------------------------------------------------------------------

test_that("basis_slepian creates spec_slepian object with defaults", {
  spec <- basis_slepian()

  expect_s3_class(spec, "spec_slepian")
  expect_equal(spec$k, 3)
  expect_equal(spec$type, "laplacian")
})

test_that("basis_slepian creates spec_slepian with custom parameters", {
  spec <- basis_slepian(k = 5, type = "adjacency")

  expect_s3_class(spec, "spec_slepian")
  expect_equal(spec$k, 5)
  expect_equal(spec$type, "adjacency")
})

test_that("basis_slepian handles edge case k = 1", {
  spec <- basis_slepian(k = 1)

  expect_s3_class(spec, "spec_slepian")
  expect_equal(spec$k, 1)
})

test_that("basis_slepian handles large k values", {
  spec <- basis_slepian(k = 100)

  expect_s3_class(spec, "spec_slepian")
  expect_equal(spec$k, 100)
})

test_that("basis_pca creates spec_pca object with defaults", {
  spec <- basis_pca()

  expect_s3_class(spec, "spec_pca")
  expect_equal(spec$k, 3)
  expect_false(spec$whiten)
})

test_that("basis_pca creates spec_pca with custom parameters", {
  spec <- basis_pca(k = 10, whiten = TRUE)

  expect_s3_class(spec, "spec_pca")
  expect_equal(spec$k, 10)
  expect_true(spec$whiten)
})

test_that("basis_pca handles k = 1", {
  spec <- basis_pca(k = 1, whiten = FALSE)

  expect_s3_class(spec, "spec_pca")
  expect_equal(spec$k, 1)
})

test_that("basis_flat creates spec_flat object with empty list", {
  spec <- basis_flat()

  expect_s3_class(spec, "spec_flat")
  expect_length(spec, 0)
})

# -----------------------------------------------------------------------------
# Tests for GraphReduction base class
# -----------------------------------------------------------------------------

test_that("GraphReduction class can be instantiated with valid slots", {
  mask_vol <- make_test_mask_vol()

  red <- new("GraphReduction", mask = mask_vol, info = list())

  expect_s4_class(red, "GraphReduction")
  expect_s4_class(red@mask, "LogicalNeuroVol")
  expect_type(red@info, "list")
})

test_that("GraphReduction stores metadata in info slot", {
  mask_vol <- make_test_mask_vol()
  meta <- list(source = "test", version = "1.0")

  red <- new("GraphReduction", mask = mask_vol, info = meta)

  expect_equal(red@info$source, "test")
  expect_equal(red@info$version, "1.0")
})

test_that("GraphReduction default lift method throws informative error", {
  mask_vol <- make_test_mask_vol()
  red <- new("GraphReduction", mask = mask_vol, info = list())
  spec <- basis_pca()

  expect_error(
    lift(red, spec),
    "No lift\\(\\) implementation"
  )
})

test_that("GraphReduction default lift method works with any basis_spec",
  {
  mask_vol <- make_test_mask_vol()
  red <- new("GraphReduction", mask = mask_vol, info = list())

  # Test with various spec types
  expect_error(lift(red, basis_slepian()), "No lift\\(\\) implementation")
  expect_error(lift(red, basis_pca()), "No lift\\(\\) implementation")
  expect_error(lift(red, basis_flat()), "No lift\\(\\) implementation")
})

# -----------------------------------------------------------------------------
# Tests for ClusterReduction class
# -----------------------------------------------------------------------------

test_that("ClusterReduction class can be instantiated with valid slots", {
  mask_vol <- make_test_mask_vol()
  n_vox <- sum(as.array(mask_vol))
  map <- seq_len(n_vox)
  ids <- unique(map)

  red <- new("ClusterReduction",
             mask = mask_vol,
             info = list(),
             map = as.integer(map),
             cluster_ids = as.integer(ids))

  expect_s4_class(red, "ClusterReduction")
  expect_s4_class(red, "GraphReduction")  # Inherits from GraphReduction
  expect_equal(length(red@map), n_vox)
  expect_equal(red@cluster_ids, seq_len(n_vox))
})

test_that("ClusterReduction stores cluster map correctly", {
  mask_vol <- make_test_mask_vol(dims = c(2, 2, 2))
  # Create clusters: voxels 1-4 in cluster 1, voxels 5-8 in cluster 2
  map <- as.integer(c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L))
  ids <- as.integer(c(1L, 2L))

  red <- new("ClusterReduction",
             mask = mask_vol,
             info = list(),
             map = map,
             cluster_ids = ids)

  expect_equal(red@map, map)
  expect_equal(red@cluster_ids, ids)
  expect_equal(sum(red@map == 1L), 4L)
  expect_equal(sum(red@map == 2L), 4L)
})

test_that("ClusterReduction works with single cluster (all voxels)", {
  mask_vol <- make_test_mask_vol(dims = c(2, 2, 2))
  n_vox <- 8L
  map <- rep(1L, n_vox)
  ids <- 1L

  red <- new("ClusterReduction",
             mask = mask_vol,
             info = list(),
             map = as.integer(map),
             cluster_ids = as.integer(ids))

  expect_s4_class(red, "ClusterReduction")
  expect_equal(length(unique(red@map)), 1L)
  expect_equal(red@cluster_ids, 1L)
})

test_that("ClusterReduction works with each voxel as its own cluster", {
  mask_vol <- make_test_mask_vol(dims = c(2, 2, 2))
  n_vox <- 8L
  map <- seq_len(n_vox)
  ids <- seq_len(n_vox)

  red <- new("ClusterReduction",
             mask = mask_vol,
             info = list(),
             map = as.integer(map),
             cluster_ids = as.integer(ids))

  expect_equal(length(unique(red@map)), n_vox)
  expect_equal(red@cluster_ids, seq_len(n_vox))
})

test_that("ClusterReduction inherits from GraphReduction", {
  mask_vol <- make_test_mask_vol()
  n_vox <- sum(as.array(mask_vol))

  red <- new("ClusterReduction",
             mask = mask_vol,
             info = list(test = TRUE),
             map = as.integer(seq_len(n_vox)),
             cluster_ids = as.integer(seq_len(n_vox)))

  expect_true(is(red, "GraphReduction"))
  expect_s4_class(red@mask, "LogicalNeuroVol")
  expect_true(red@info$test)
})

test_that("make_cluster_reduction helper creates valid ClusterReduction", {
  mask <- array(TRUE, dim = c(2, 2, 2))
  map <- seq_len(sum(mask))

  red <- make_cluster_reduction(mask, map)

  expect_s4_class(red, "ClusterReduction")
  expect_equal(length(red@map), 8L)
  expect_equal(red@cluster_ids, seq_len(8L))
})

test_that("make_cluster_reduction works with LogicalNeuroVol input", {
  mask_vol <- make_test_mask_vol(dims = c(2, 2, 2))
  map <- seq_len(8L)

  red <- make_cluster_reduction(mask_vol, map)

  expect_s4_class(red, "ClusterReduction")
  expect_s4_class(red@mask, "LogicalNeuroVol")
})

test_that("make_cluster_reduction handles non-contiguous cluster IDs", {
  mask <- array(TRUE, dim = c(2, 2, 1))
  # Non-contiguous IDs: 1, 5, 10
  map <- c(1L, 1L, 5L, 10L)

  red <- make_cluster_reduction(mask, map)

  expect_equal(sort(red@cluster_ids), c(1L, 5L, 10L))
})

test_that("make_cluster_reduction converts numeric map to integer", {
  mask <- array(TRUE, dim = c(2, 2, 1))
  map <- c(1.0, 1.0, 2.0, 2.0)  # Numeric, not integer

  red <- make_cluster_reduction(mask, map)

  expect_type(red@map, "integer")
  expect_type(red@cluster_ids, "integer")
})

# -----------------------------------------------------------------------------
# Tests for CoarsenedReduction class
# -----------------------------------------------------------------------------

test_that("CoarsenedReduction class can be instantiated with valid slots", {
  mask_vol <- make_test_mask_vol(dims = c(2, 2, 2))
  n_fine <- 8L
  n_coarse <- 4L

  # Create a simple prolongation matrix (fine x coarse)
  P_mat <- Matrix::sparseMatrix(
    i = 1:8,
    j = rep(1:4, each = 2),
    x = rep(1, 8),
    dims = c(n_fine, n_coarse)
  )

  # Create coarse adjacency
  coarse_adj <- Matrix::sparseMatrix(
    i = c(1, 2, 3, 4, 1, 2, 3, 4),
    j = c(2, 1, 4, 3, 3, 4, 1, 2),
    x = rep(1, 8),
    dims = c(n_coarse, n_coarse)
  )

  red <- new("CoarsenedReduction",
             mask = mask_vol,
             info = list(),
             P_matrix = P_mat,
             coarse_adj = coarse_adj)

  expect_s4_class(red, "CoarsenedReduction")
  expect_s4_class(red, "GraphReduction")
  expect_s4_class(red@P_matrix, "dgCMatrix")
  expect_s4_class(red@coarse_adj, "dgCMatrix")
})

test_that("CoarsenedReduction inherits from GraphReduction", {
  mask_vol <- make_test_mask_vol(dims = c(2, 2, 2))

  P_mat <- Matrix::sparseMatrix(
    i = 1:8, j = rep(1:2, each = 4), x = rep(1, 8),
    dims = c(8, 2)
  )
  coarse_adj <- Matrix::sparseMatrix(
    i = c(1, 2), j = c(2, 1), x = c(1, 1), dims = c(2, 2)
  )

  red <- new("CoarsenedReduction",
             mask = mask_vol,
             info = list(coarsening_type = "test"),
             P_matrix = P_mat,
             coarse_adj = coarse_adj)

  expect_true(is(red, "GraphReduction"))
  expect_equal(red@info$coarsening_type, "test")
})

test_that("CoarsenedReduction P_matrix dimensions are correct", {
  mask_vol <- make_test_mask_vol(dims = c(2, 2, 2))
  n_fine <- 8L
  n_coarse <- 2L

  P_mat <- Matrix::sparseMatrix(
    i = 1:8, j = rep(1:2, each = 4), x = rep(0.5, 8),
    dims = c(n_fine, n_coarse)
  )
  coarse_adj <- Matrix::sparseMatrix(i = 1:2, j = 2:1, x = c(1, 1), dims = c(2, 2))

  red <- new("CoarsenedReduction",
             mask = mask_vol,
             info = list(),
             P_matrix = P_mat,
             coarse_adj = coarse_adj)

  expect_equal(nrow(red@P_matrix), n_fine)
  expect_equal(ncol(red@P_matrix), n_coarse)
  expect_equal(nrow(red@coarse_adj), n_coarse)
  expect_equal(ncol(red@coarse_adj), n_coarse)
})

test_that("CoarsenedReduction handles identity coarsening (no reduction)", {
  mask_vol <- make_test_mask_vol(dims = c(2, 2, 1))
  n <- 4L

  # Identity prolongation (fine == coarse)
  P_mat <- Matrix::sparseMatrix(
    i = 1:n, j = 1:n, x = rep(1, n), dims = c(n, n)
  )
  # Simple adjacency (explicitly make it non-symmetric to get dgCMatrix)
  coarse_adj <- Matrix::sparseMatrix(
    i = c(1, 2, 3, 2, 3, 4), j = c(2, 3, 4, 1, 2, 3), x = rep(1, 6),
    dims = c(n, n)
  )

  red <- new("CoarsenedReduction",
             mask = mask_vol,
             info = list(),
             P_matrix = P_mat,
             coarse_adj = coarse_adj)

  expect_equal(nrow(red@P_matrix), ncol(red@P_matrix))
})

test_that("CoarsenedReduction works with empty coarse_adj", {
  mask_vol <- make_test_mask_vol(dims = c(2, 2, 1))

  P_mat <- Matrix::sparseMatrix(
    i = 1:4, j = c(1, 1, 2, 2), x = rep(1, 4), dims = c(4, 2)
  )
  # Empty adjacency
  coarse_adj <- Matrix::sparseMatrix(
    i = integer(0), j = integer(0), x = numeric(0), dims = c(2, 2)
  )

  red <- new("CoarsenedReduction",
             mask = mask_vol,
             info = list(),
             P_matrix = P_mat,
             coarse_adj = coarse_adj)

  expect_equal(Matrix::nnzero(red@coarse_adj), 0L)
})

# -----------------------------------------------------------------------------
# Tests for lift generic and default method
# -----------------------------------------------------------------------------

test_that("lift generic is defined", {
  expect_true(isGeneric("lift"))
})

test_that("lift default method provides informative error message", {
  mask_vol <- make_test_mask_vol()
  red <- new("GraphReduction", mask = mask_vol, info = list())
  spec <- basis_slepian()

  err <- expect_error(lift(red, spec))
  expect_match(err$message, "external method")
})

test_that("lift with data = NULL works for default method (still errors)", {
  mask_vol <- make_test_mask_vol()
  red <- new("GraphReduction", mask = mask_vol, info = list())
  spec <- basis_pca()

  expect_error(lift(red, spec, data = NULL), "No lift\\(\\) implementation")
})

test_that("lift with additional ... arguments passes through", {
  mask_vol <- make_test_mask_vol()
  red <- new("GraphReduction", mask = mask_vol, info = list())
  spec <- basis_flat()

  # Extra arguments should not change the error behavior
  expect_error(
    lift(red, spec, data = NULL, extra_arg = 42),
    "No lift\\(\\) implementation"
  )
})

# -----------------------------------------------------------------------------
# Tests for ClusterReduction lift methods (requires rgsp/RSpectra)
# -----------------------------------------------------------------------------

test_that("lift for ClusterReduction with spec_slepian works", {
  skip_if_not_installed("rgsp")
  skip_if_not_installed("RSpectra")

  mask <- array(TRUE, dim = c(2, 2, 2))
  map <- seq_len(sum(mask))
  red <- make_cluster_reduction(mask, map)
  spec <- basis_slepian(k = 2)

  L <- lift(red, spec, k_neighbors = 3L)

  expect_s4_class(L, "dgCMatrix")
  expect_equal(nrow(L), length(map))
  expect_gt(ncol(L), 0)
})

test_that("lift for ClusterReduction handles single-voxel cluster", {
  skip_if_not_installed("rgsp")
  skip_if_not_installed("RSpectra")

  mask <- array(TRUE, dim = c(2, 2, 2))
  # All voxels in separate clusters (8 clusters, 1 voxel each)
  map <- seq_len(sum(mask))
  red <- make_cluster_reduction(mask, map)
  spec <- basis_slepian(k = 1)

  L <- lift(red, spec, k_neighbors = 3L)

  expect_s4_class(L, "dgCMatrix")
  expect_equal(nrow(L), 8L)
})

test_that("lift for ClusterReduction with multiple clusters", {
  skip_if_not_installed("rgsp")
  skip_if_not_installed("RSpectra")

  mask <- array(TRUE, dim = c(2, 2, 2))
  # 2 clusters of 4 voxels each
  map <- rep(c(1L, 2L), each = 4)
  red <- make_cluster_reduction(mask, map)
  spec <- basis_slepian(k = 2)

  L <- lift(red, spec, k_neighbors = 2L)

  expect_s4_class(L, "dgCMatrix")
  expect_equal(nrow(L), 8L)
})

# -----------------------------------------------------------------------------
# Tests for edge cases
# -----------------------------------------------------------------------------

test_that("basis_slepian with k = 0 creates valid object", {
  # Edge case: zero components requested
  spec <- basis_slepian(k = 0)

  expect_s3_class(spec, "spec_slepian")
  expect_equal(spec$k, 0)
})

test_that("basis_pca with k = 0 creates valid object", {
  spec <- basis_pca(k = 0)

  expect_s3_class(spec, "spec_pca")
  expect_equal(spec$k, 0)
})

test_that("ClusterReduction with partial mask works", {
  mask_vol <- make_test_mask_vol(dims = c(3, 3, 3), all_true = FALSE)
  n_vox <- sum(as.array(mask_vol))
  map <- seq_len(n_vox)

  red <- new("ClusterReduction",
             mask = mask_vol,
             info = list(),
             map = as.integer(map),
             cluster_ids = as.integer(unique(map)))

  expect_s4_class(red, "ClusterReduction")
  expect_equal(length(red@map), n_vox)
})

test_that("GraphReduction with empty info list works", {
  mask_vol <- make_test_mask_vol()

  red <- new("GraphReduction", mask = mask_vol, info = list())

  expect_length(red@info, 0)
})

test_that("ClusterReduction with complex info metadata", {
  mask_vol <- make_test_mask_vol(dims = c(2, 2, 2))
  map <- seq_len(8L)
  info <- list(
    atlas_name = "test_atlas",
    n_parcels = 8L,
    resolution = 2.0,
    nested = list(a = 1, b = 2)
  )

  red <- new("ClusterReduction",
             mask = mask_vol,
             info = info,
             map = as.integer(map),
             cluster_ids = as.integer(unique(map)))

  expect_equal(red@info$atlas_name, "test_atlas")
  expect_equal(red@info$n_parcels, 8L)
  expect_equal(red@info$nested$a, 1)
})

# -----------------------------------------------------------------------------
# Tests for type coercion and validation
# -----------------------------------------------------------------------------

test_that("make_cluster_reduction coerces array mask to LogicalNeuroVol", {
  mask <- array(TRUE, dim = c(2, 2, 2))
  map <- seq_len(8L)

  red <- make_cluster_reduction(mask, map)

  expect_s4_class(red@mask, "LogicalNeuroVol")
})

test_that("make_cluster_reduction handles numeric array mask", {
  mask <- array(1, dim = c(2, 2, 2))  # Numeric 1s instead of TRUE
  map <- seq_len(8L)

  # The function converts to LogicalNeuroVol which should coerce to logical
  red <- make_cluster_reduction(mask, map)

  expect_s4_class(red@mask, "LogicalNeuroVol")
})

test_that("ClusterReduction cluster_ids are sorted and unique", {
  mask <- array(TRUE, dim = c(2, 2, 2))
  # Unsorted, with duplicates in different positions
  map <- c(3L, 1L, 2L, 1L, 3L, 2L, 1L, 3L)

  red <- make_cluster_reduction(mask, map)

  expect_equal(red@cluster_ids, c(1L, 2L, 3L))  # Sorted unique
})

# -----------------------------------------------------------------------------
# Integration tests: reduction -> lift workflow
# -----------------------------------------------------------------------------

test_that("full workflow: create ClusterReduction and lift slepian basis", {
  skip_if_not_installed("rgsp")
  skip_if_not_installed("RSpectra")

  # Setup
  mask <- array(TRUE, dim = c(3, 3, 2))
  n_vox <- sum(mask)
  map <- seq_len(n_vox)

  # Create reduction
  red <- make_cluster_reduction(mask, map)
  expect_s4_class(red, "ClusterReduction")

  # Lift with slepian basis
  spec <- basis_slepian(k = 2)
  L <- lift(red, spec, k_neighbors = 4L)

  expect_s4_class(L, "Matrix")
  expect_equal(nrow(L), n_vox)
})

test_that("full workflow: create ClusterReduction with atlas-like parcellation", {
  skip_if_not_installed("rgsp")
  skip_if_not_installed("RSpectra")

  # Simulate an atlas with 4 parcels
  mask <- array(TRUE, dim = c(4, 4, 2))
  n_vox <- sum(mask)
  # 4 parcels, 8 voxels each
  map <- rep(1:4, each = 8)

  red <- make_cluster_reduction(mask, map)
  expect_equal(length(red@cluster_ids), 4L)

  spec <- basis_slepian(k = 3)
  L <- lift(red, spec, k_neighbors = 3L)

  expect_s4_class(L, "Matrix")
  expect_equal(nrow(L), n_vox)
  # Each of 4 parcels should contribute up to k=3 components
  expect_lte(ncol(L), 4 * 3)  # At most 12 columns
})
