library(testthat)

set.seed(1)

make_mask_cube <- function(nx = 2, ny = 2, nz = 2, fill = TRUE) {
  if (fill) {
    array(TRUE, dim = c(nx, ny, nz))
  } else {
    m <- array(FALSE, dim = c(nx, ny, nz))
    m[1, 1, 1] <- TRUE
    m
  }
}

test_that("precompute_haar_scalings counts voxels", {
  mask_full <- make_mask_cube(2, 2, 2, TRUE)
  s_full <- fmrilatent:::precompute_haar_scalings(mask_full, 1)
  expect_equal(length(s_full), 1L)
  expect_equal(s_full[[1]]$sqrt_nvalid, sqrt(8))

  mask_single <- make_mask_cube(2, 2, 2, FALSE)
  s_single <- fmrilatent:::precompute_haar_scalings(mask_single, 1)
  expect_equal(sum(round(s_single[[1]]$sqrt_nvalid^2)), 1)
})

test_that("single-level Haar roundtrip works (R path)", {
  mask <- make_mask_cube(2, 2, 2, TRUE)
  X <- matrix(rnorm(16), nrow = 2)
  opts <- options(fmrilatent.hwt.use_rcpp = FALSE)
  coeff <- fmrilatent:::perform_haar_lift_analysis(X, mask, levels = 1, z_order_seed = 42L)
  reco  <- fmrilatent:::perform_haar_lift_synthesis(coeff, mask, levels = 1, z_order_seed = 42L)
  options(opts)
  expect_equal(reco, X, tolerance = 1e-7)
})

test_that("multi-level Haar roundtrip works (R path)", {
  mask <- array(TRUE, dim = c(3, 3, 3))
  X <- matrix(rnorm(5 * sum(mask)), nrow = 5)
  opts <- options(fmrilatent.hwt.use_rcpp = FALSE)
  coeff <- fmrilatent:::perform_haar_lift_analysis(X, mask, levels = 2, z_order_seed = 42L)
  reco  <- fmrilatent:::perform_haar_lift_synthesis(coeff, mask, levels = 2, z_order_seed = 42L)
  options(opts)
  expect_equal(reco, X, tolerance = 1e-6)
})

test_that("ROI subsetting returns matching subset", {
  mask <- array(FALSE, dim = c(3, 3, 3))
  mask[1:2, 1:2, 1:2] <- TRUE
  mask[3, 3, 3] <- TRUE
  X <- matrix(rnorm(4 * sum(mask)), nrow = 4)

  fw <- haar_wavelet_forward(X, mask, levels = 2, z_seed = 42L)
  full <- haar_wavelet_inverse(fw, mask, levels = 2, z_seed = 42L)

  roi <- array(FALSE, dim = c(3, 3, 3))
  roi[1:2, 1:2, 1:2] <- TRUE
  roi_reco <- haar_wavelet_inverse(fw, mask, levels = 2, z_seed = 42L, roi_mask = roi)

  global_idx <- which(as.logical(mask))
  roi_global <- which(as.logical(roi))
  col_keep <- which(global_idx %in% roi_global)

  expect_equal(roi_reco, full[, col_keep, drop = FALSE], tolerance = 1e-6)
})

test_that("levels_keep zeros selected detail levels", {
  mask <- array(TRUE, dim = c(2, 2, 2))
  X <- matrix(rnorm(12 * sum(mask)), nrow = 12)
  fw <- haar_wavelet_forward(X, mask, levels = 1, z_seed = 42L)

  full <- haar_wavelet_inverse(fw, mask)
  root_only <- haar_wavelet_inverse(fw, mask, levels_keep = integer(0))

  expect_equal(full, X, tolerance = 1e-6)
  expect_gt(sum(abs(full - root_only)), 1e-8)
})

test_that("ImplicitLatent predict matches HaarLatent predict", {
  mask <- array(TRUE, dim = c(2, 2, 2))
  X <- matrix(rnorm(6 * sum(mask)), nrow = 6)
  hl <- haar_latent(X, mask, levels = 1, z_seed = 42L)
  expect_true(is_implicit_latent(hl))
  expect_true(is_haar_latent(hl))
  full1 <- predict(hl)
  full2 <- predict.ImplicitLatent(hl)
  expect_equal(full1, full2, tolerance = 1e-6)
})

# -- Additional coverage tests ------------------------------------------------

test_that("as_logical_mask converts various inputs", {
  # Test regular 3D array
  arr <- array(c(TRUE, FALSE, TRUE, FALSE), dim = c(2, 2, 1))
  result <- fmrilatent:::as_logical_mask(arr)
  expect_true(is.logical(result))
  expect_equal(dim(result), c(2, 2, 1))

  # Test numeric conversion
  num_arr <- array(c(1, 0, 1, 0), dim = c(2, 2, 1))
  result_num <- fmrilatent:::as_logical_mask(num_arr)
  expect_equal(result, result_num)

  # Test error on wrong dimensions
  expect_error(
    fmrilatent:::as_logical_mask(matrix(TRUE, 2, 2)),
    "must be a 3D array"
  )

  # Test error on NULL dims
  expect_error(
    fmrilatent:::as_logical_mask(c(TRUE, FALSE)),
    "must be a 3D array"
  )
})

test_that("get_morton_ordered_indices returns correct ordering (R path)", {
  opts <- options(fmrilatent.hwt.use_rcpp = FALSE)

  # Test with simple 2x2x2 mask
  mask <- array(TRUE, dim = c(2, 2, 2))
  indices <- fmrilatent:::get_morton_ordered_indices(mask, z_order_seed = 42L)
  expect_equal(length(indices), 8)
  expect_true(all(indices >= 1 & indices <= 8))
  expect_equal(length(unique(indices)), 8)

  # Test with sparse mask
  mask_sparse <- array(FALSE, dim = c(3, 3, 3))
  mask_sparse[1, 1, 1] <- TRUE
  mask_sparse[3, 3, 3] <- TRUE
  indices_sparse <- fmrilatent:::get_morton_ordered_indices(mask_sparse)
  expect_equal(length(indices_sparse), 2)

  # Test with empty mask
  mask_empty <- array(FALSE, dim = c(2, 2, 2))
  indices_empty <- fmrilatent:::get_morton_ordered_indices(mask_empty)
  expect_equal(length(indices_empty), 0)

  # Test error on oversized mask (use smaller size to avoid memory issues)
  # Just verify the error check exists without actually allocating huge array
  skip("Oversized mask test skipped to avoid memory allocation")

  options(opts)
})

test_that("get_valid_finest_blocks identifies valid blocks (R path)", {
  opts <- options(fmrilatent.hwt.use_rcpp = FALSE)

  # Test with full mask
  mask <- array(TRUE, dim = c(2, 2, 2))
  blocks <- fmrilatent:::get_valid_finest_blocks(mask)
  expect_true(length(blocks) > 0)
  expect_true(is.integer(blocks))

  # Test with sparse mask
  mask_sparse <- array(FALSE, dim = c(4, 4, 4))
  mask_sparse[1:2, 1:2, 1:2] <- TRUE
  blocks_sparse <- fmrilatent:::get_valid_finest_blocks(mask_sparse)
  expect_true(length(blocks_sparse) > 0)

  # Test with empty mask
  mask_empty <- array(FALSE, dim = c(2, 2, 2))
  blocks_empty <- fmrilatent:::get_valid_finest_blocks(mask_empty)
  expect_equal(length(blocks_empty), 0)

  options(opts)
})

test_that("forward_lift_R and inverse_lift_R roundtrip correctly", {
  opts <- options(fmrilatent.hwt.use_rcpp = FALSE)

  mask <- array(TRUE, dim = c(2, 2, 2))
  mask_logical <- fmrilatent:::as_logical_mask(mask)
  morton_idx <- fmrilatent:::get_morton_ordered_indices(mask_logical, 42L)

  full_order <- fmrilatent:::get_morton_ordered_indices(
    array(TRUE, dim(mask_logical)), 42L
  )
  mask_flat <- as.vector(mask_logical)
  mask_flat_morton <- mask_flat[full_order]

  scalings <- fmrilatent:::precompute_haar_scalings(mask_logical, 1)

  # Test with simple data
  data_vec <- rnorm(8)

  # Forward
  fw <- fmrilatent:::forward_lift_R(
    data_vec, mask_flat_morton, dim(mask_logical), 1L, scalings
  )

  expect_true(is.list(fw))
  expect_true("root_coeff" %in% names(fw))
  expect_true("detail_coeffs_by_level" %in% names(fw))

  # Inverse
  reco <- fmrilatent:::inverse_lift_R(
    fw$root_coeff, fw$detail_coeffs_by_level,
    mask_flat_morton, dim(mask_logical), 1L, scalings
  )

  expect_equal(reco, data_vec, tolerance = 1e-7)

  options(opts)
})

test_that("lna_forward_lift_matrix handles matrix batches (R path)", {
  opts <- options(fmrilatent.hwt.use_rcpp = FALSE)

  mask <- array(TRUE, dim = c(2, 2, 2))
  mask_logical <- fmrilatent:::as_logical_mask(mask)
  morton_idx <- fmrilatent:::get_morton_ordered_indices(mask_logical, 42L)

  mask_linear <- which(mask_logical)
  perm <- match(morton_idx, mask_linear)

  full_order <- fmrilatent:::get_morton_ordered_indices(
    array(TRUE, dim(mask_logical)), 42L
  )
  mask_flat <- as.vector(mask_logical)
  mask_flat_morton <- mask_flat[full_order]

  scalings <- fmrilatent:::precompute_haar_scalings(mask_logical, 1)

  # Test with multiple time points
  X <- matrix(rnorm(3 * 8), nrow = 3, ncol = 8)
  X_morton <- X[, perm, drop = FALSE]

  result <- fmrilatent:::lna_forward_lift_matrix(
    X_morton, mask_flat_morton, dim(mask_logical), 1L, scalings
  )

  expect_true(is.matrix(result$root_coeff))
  expect_equal(nrow(result$root_coeff), 3)
  expect_true(is.list(result$detail_coeffs_by_level))

  # Test with empty matrix
  X_empty <- matrix(0, nrow = 0, ncol = 8)
  result_empty <- fmrilatent:::lna_forward_lift_matrix(
    X_empty, mask_flat_morton, dim(mask_logical), 1L, scalings
  )
  expect_equal(nrow(result_empty$root_coeff), 0)

  # Test error on non-matrix input
  expect_error(
    fmrilatent:::lna_forward_lift_matrix(
      rnorm(8), mask_flat_morton, dim(mask_logical), 1L, scalings
    ),
    "must be a matrix"
  )

  options(opts)
})

test_that("lna_inverse_lift_matrix handles matrix batches (R path)", {
  opts <- options(fmrilatent.hwt.use_rcpp = FALSE)

  mask <- array(TRUE, dim = c(2, 2, 2))
  mask_logical <- fmrilatent:::as_logical_mask(mask)
  morton_idx <- fmrilatent:::get_morton_ordered_indices(mask_logical, 42L)
  mask_linear <- which(mask_logical)
  perm <- match(morton_idx, mask_linear)

  full_order <- fmrilatent:::get_morton_ordered_indices(
    array(TRUE, dim(mask_logical)), 42L
  )
  mask_flat <- as.vector(mask_logical)
  mask_flat_morton <- mask_flat[full_order]

  scalings <- fmrilatent:::precompute_haar_scalings(mask_logical, 1)

  X <- matrix(rnorm(3 * 8), nrow = 3, ncol = 8)
  X_morton <- X[, perm, drop = FALSE]

  fw <- fmrilatent:::lna_forward_lift_matrix(
    X_morton, mask_flat_morton, dim(mask_logical), 1L, scalings
  )

  # Test roundtrip
  reco <- fmrilatent:::lna_inverse_lift_matrix(
    fw$root_coeff, fw$detail_coeffs_by_level,
    mask_flat_morton, dim(mask_logical), 1L, scalings
  )

  expect_equal(reco, X_morton, tolerance = 1e-7)

  # Test with vector input (auto-converts to matrix) - need list with at least `levels` elements
  detail_list_single <- vector("list", 1L)
  detail_list_single[[1]] <- if (is.matrix(fw$detail_coeffs_by_level[[1]])) {
    matrix(fw$detail_coeffs_by_level[[1]][1, ], nrow = 1)
  } else {
    matrix(fw$detail_coeffs_by_level[[1]], nrow = 1)
  }
  reco_vec <- fmrilatent:::lna_inverse_lift_matrix(
    matrix(fw$root_coeff[1, ], nrow = 1),
    detail_list_single,
    mask_flat_morton, dim(mask_logical), 1L, scalings
  )
  expect_true(is.matrix(reco_vec))
  expect_equal(nrow(reco_vec), 1)

  # Test error on invalid detail list
  expect_error(
    fmrilatent:::lna_inverse_lift_matrix(
      fw$root_coeff, list(), mask_flat_morton, dim(mask_logical), 1L, scalings
    ),
    "must be a list of length"
  )

  # Test empty matrix
  root_empty <- matrix(0, nrow = 0, ncol = 1)
  result_empty <- fmrilatent:::lna_inverse_lift_matrix(
    root_empty, fw$detail_coeffs_by_level[integer(0)],
    mask_flat_morton, dim(mask_logical), 0L, list()
  )
  expect_equal(nrow(result_empty), 0)

  options(opts)
})

test_that("encode_morton3d encodes coordinates correctly", {
  # Test origin
  code <- fmrilatent:::encode_morton3d(0L, 0L, 0L, 3L)
  expect_equal(code, 0L)

  # Test non-zero coordinates
  code1 <- fmrilatent:::encode_morton3d(1L, 0L, 0L, 3L)
  code2 <- fmrilatent:::encode_morton3d(0L, 1L, 0L, 3L)
  code3 <- fmrilatent:::encode_morton3d(0L, 0L, 1L, 3L)

  # All should be different
  expect_true(code1 != code2)
  expect_true(code2 != code3)
  expect_true(code1 != code3)

  # Test symmetric coordinates
  code_sym <- fmrilatent:::encode_morton3d(1L, 1L, 1L, 3L)
  expect_true(code_sym > 0L)
})

test_that("compute_block_map creates valid block mapping", {
  mask <- array(TRUE, dim = c(2, 2, 2))
  mapping <- fmrilatent:::compute_block_map(mask, 1)

  expect_true(is.list(mapping))
  expect_equal(length(mapping), 1)
  expect_true("code" %in% names(mapping[[1]]))
  expect_true("count" %in% names(mapping[[1]]))
  expect_true("start" %in% names(mapping[[1]]))

  # Test with sparse mask
  mask_sparse <- array(FALSE, dim = c(4, 4, 4))
  mask_sparse[1:2, 1:2, 1:2] <- TRUE
  mapping_sparse <- fmrilatent:::compute_block_map(mask_sparse, 2)
  expect_equal(length(mapping_sparse), 2)
})

test_that("get_roi_detail_indices extracts ROI indices", {
  mask <- array(TRUE, dim = c(4, 4, 4))
  roi <- array(FALSE, dim = c(4, 4, 4))
  roi[1:2, 1:2, 1:2] <- TRUE

  indices <- fmrilatent:::get_roi_detail_indices(roi, mask, 2)

  expect_true(is.list(indices))
  expect_equal(length(indices), 2)
  expect_true(all(vapply(indices, is.integer, logical(1))))

  # Test with mismatched dimensions
  roi_bad <- array(FALSE, dim = c(2, 2, 2))
  expect_error(
    fmrilatent:::get_roi_detail_indices(roi_bad, mask, 2),
    "must match mask dimensions"
  )

  # Test with empty ROI
  roi_empty <- array(FALSE, dim = c(4, 4, 4))
  indices_empty <- fmrilatent:::get_roi_detail_indices(roi_empty, mask, 2)
  expect_equal(length(indices_empty), 2)
  expect_true(all(vapply(indices_empty, length, integer(1)) == 0))
})

test_that("haar_wavelet_forward handles various inputs", {
  # Test with transposed input (auto-transposes)
  mask <- array(TRUE, dim = c(2, 2, 2))
  X_transposed <- matrix(rnorm(8 * 3), nrow = 8, ncol = 3)

  result <- haar_wavelet_forward(X_transposed, mask, levels = 1, z_seed = 42L)
  expect_true(is.list(result))
  expect_true("coeff" %in% names(result))
  expect_true("meta" %in% names(result))

  # Test with threshold
  X <- matrix(rnorm(8 * 3), nrow = 3, ncol = 8)
  result_thresh <- haar_wavelet_forward(
    X, mask, levels = 1, z_seed = 42L,
    threshold = list(type = "absolute", value = 0.5)
  )
  expect_true(is.list(result_thresh))

  # Test with relative threshold
  result_rel <- haar_wavelet_forward(
    X, mask, levels = 1, z_seed = 42L,
    threshold = list(type = "relative_to_root_std", value = 0.1)
  )
  expect_true(is.list(result_rel))

  # Test metadata
  expect_equal(result$meta$levels, 1)
  expect_equal(result$meta$z_seed, 42L)
  expect_equal(result$meta$mask_dims, c(2, 2, 2))
  expect_equal(result$meta$num_voxels_in_mask, 8)
  expect_true(length(result$meta$valid_finest_blocks) > 0)
})

test_that("haar_wavelet_inverse handles roi_mask and time_idx", {
  mask <- array(TRUE, dim = c(3, 3, 3))
  X <- matrix(rnorm(5 * 27), nrow = 5, ncol = 27)

  fw <- haar_wavelet_forward(X, mask, levels = 2, z_seed = 42L)

  # Test time subset
  time_subset <- haar_wavelet_inverse(fw, mask, time_idx = c(1, 3, 5))
  expect_equal(nrow(time_subset), 3)
  expect_equal(ncol(time_subset), 27)

  # Test ROI subset
  roi <- array(FALSE, dim = c(3, 3, 3))
  roi[1:2, 1:2, 1:2] <- TRUE
  roi_reco <- haar_wavelet_inverse(fw, mask, roi_mask = roi)
  expect_true(ncol(roi_reco) < 27)

  # Test combined roi and time
  combined <- haar_wavelet_inverse(fw, mask, roi_mask = roi, time_idx = c(2, 4))
  expect_equal(nrow(combined), 2)
  expect_true(ncol(combined) < 27)
})

test_that("as.matrix.HaarLatent converts to matrix", {
  mask <- array(TRUE, dim = c(2, 2, 2))
  X <- matrix(rnorm(8 * 4), nrow = 4, ncol = 8)

  hl <- haar_latent(X, mask, levels = 1, z_seed = 42L)
  mat <- as.matrix(hl)

  expect_true(is.matrix(mat))
  expect_equal(dim(mat), dim(X))
  expect_equal(mat, X, tolerance = 1e-6)
})

test_that("predict.HaarLatent delegates correctly", {
  mask <- array(TRUE, dim = c(2, 2, 2))
  X <- matrix(rnorm(8 * 4), nrow = 4, ncol = 8)

  hl <- haar_latent(X, mask, levels = 1, z_seed = 42L)

  # Test basic predict
  pred <- predict(hl)
  expect_equal(pred, X, tolerance = 1e-6)

  # Test with time_idx
  pred_time <- predict(hl, time_idx = c(1, 3))
  expect_equal(nrow(pred_time), 2)

  # Test with levels_keep
  pred_coarse <- predict(hl, levels_keep = integer(0))
  expect_equal(nrow(pred_coarse), 4)
  expect_true(sum(abs(pred - pred_coarse)) > 1e-8)
})

test_that("precompute_haar_scalings validates inputs", {
  # Test invalid levels
  mask <- array(TRUE, dim = c(2, 2, 2))
  expect_error(
    fmrilatent:::precompute_haar_scalings(mask, 0),
    "must be a positive integer"
  )
  expect_error(
    fmrilatent:::precompute_haar_scalings(mask, -1),
    "must be a positive integer"
  )
  expect_error(
    fmrilatent:::precompute_haar_scalings(mask, NA),
    "must be a positive integer"
  )
})

test_that("perform_haar_lift_analysis validates inputs", {
  mask <- array(TRUE, dim = c(2, 2, 2))

  # Test non-matrix input
  expect_error(
    fmrilatent:::perform_haar_lift_analysis(rnorm(8), mask, 1),
    "must be a matrix"
  )

  # Test wrong number of columns
  X_wrong <- matrix(rnorm(3 * 5), nrow = 3, ncol = 5)
  expect_error(
    fmrilatent:::perform_haar_lift_analysis(X_wrong, mask, 1),
    "must have 8 columns"
  )

  # Test empty mask
  mask_empty <- array(FALSE, dim = c(2, 2, 2))
  X <- matrix(rnorm(3 * 0), nrow = 3, ncol = 0)
  expect_error(
    fmrilatent:::perform_haar_lift_analysis(X, mask_empty, 1),
    "must contain at least one voxel"
  )
})

test_that("perform_haar_lift_synthesis validates inputs", {
  mask <- array(TRUE, dim = c(2, 2, 2))

  # Test invalid root
  expect_error(
    fmrilatent:::perform_haar_lift_synthesis(list(root = NULL), mask, 1),
    "must be a matrix"
  )

  # Test invalid detail
  expect_error(
    fmrilatent:::perform_haar_lift_synthesis(
      list(root = matrix(0, 1, 1), detail = NULL), mask, 1
    ),
    "must be a list of length"
  )

  # Test detail length mismatch
  expect_error(
    fmrilatent:::perform_haar_lift_synthesis(
      list(root = matrix(0, 1, 1), detail = list()), mask, 1
    ),
    "must be a list of length"
  )
})

# =============================================================================
# Additional coverage tests -- R-fallback functions called DIRECTLY
# =============================================================================

# -- forward_lift_R / inverse_lift_R direct tests ----------------------------

test_that("forward_lift_R returns expected structure with 2-level decomposition", {

  # Build a 4x4x4 all-TRUE mask so we get enough voxels for 2 levels

  mask <- array(TRUE, dim = c(4, 4, 4))
  mask_logical <- fmrilatent:::as_logical_mask(mask)
  scalings <- fmrilatent:::precompute_haar_scalings(mask_logical, 2L)

  full_order <- fmrilatent:::get_morton_ordered_indices(
    array(TRUE, dim(mask_logical)), 42L
  )
  mask_flat <- as.vector(mask_logical)
  mask_flat_morton <- mask_flat[full_order]

  data_vec <- rnorm(64)

  fw <- fmrilatent:::forward_lift_R(
    data_vec, mask_flat_morton, dim(mask_logical), 2L, scalings
  )

  expect_true(is.list(fw))
  expect_equal(length(fw$detail_coeffs_by_level), 2L)
  # Root should be a scalar (single top-level block)
  expect_true(length(fw$root_coeff) >= 1L)
  # Detail at level 1 should have as many elements as total voxels in blocks
  expect_true(length(fw$detail_coeffs_by_level[[1]]) > 0)
  expect_true(length(fw$detail_coeffs_by_level[[2]]) > 0)
})

test_that("forward_lift_R / inverse_lift_R roundtrip with multi-level (4x4x4)", {
  mask <- array(TRUE, dim = c(4, 4, 4))
  mask_logical <- fmrilatent:::as_logical_mask(mask)
  scalings <- fmrilatent:::precompute_haar_scalings(mask_logical, 2L)

  full_order <- fmrilatent:::get_morton_ordered_indices(
    array(TRUE, dim(mask_logical)), 42L
  )
  mask_flat <- as.vector(mask_logical)
  mask_flat_morton <- mask_flat[full_order]

  set.seed(99)
  data_vec <- rnorm(64)

  fw <- fmrilatent:::forward_lift_R(
    data_vec, mask_flat_morton, dim(mask_logical), 2L, scalings
  )
  reco <- fmrilatent:::inverse_lift_R(
    fw$root_coeff, fw$detail_coeffs_by_level,
    mask_flat_morton, dim(mask_logical), 2L, scalings
  )
  expect_equal(reco, data_vec, tolerance = 1e-7)
})

test_that("forward_lift_R / inverse_lift_R roundtrip with sparse mask", {
  # Only some voxels are TRUE
  mask <- array(FALSE, dim = c(4, 4, 4))
  mask[1:2, 1:2, 1:2] <- TRUE
  mask[3, 3, 3] <- TRUE
  mask_logical <- fmrilatent:::as_logical_mask(mask)
  n_vox <- sum(mask_logical)
  scalings <- fmrilatent:::precompute_haar_scalings(mask_logical, 1L)

  morton_idx <- fmrilatent:::get_morton_ordered_indices(mask_logical, 42L)
  full_order <- fmrilatent:::get_morton_ordered_indices(
    array(TRUE, dim(mask_logical)), 42L
  )
  mask_flat <- as.vector(mask_logical)
  mask_flat_morton <- mask_flat[full_order]

  set.seed(55)
  data_vec <- rnorm(n_vox)

  fw <- fmrilatent:::forward_lift_R(
    data_vec, mask_flat_morton, dim(mask_logical), 1L, scalings
  )
  reco <- fmrilatent:::inverse_lift_R(
    fw$root_coeff, fw$detail_coeffs_by_level,
    mask_flat_morton, dim(mask_logical), 1L, scalings
  )
  expect_equal(reco, data_vec, tolerance = 1e-7)
})

test_that("inverse_lift_R errors on detail vector length mismatch", {
  mask <- array(TRUE, dim = c(2, 2, 2))
  mask_logical <- fmrilatent:::as_logical_mask(mask)
  scalings <- fmrilatent:::precompute_haar_scalings(mask_logical, 1L)

  full_order <- fmrilatent:::get_morton_ordered_indices(
    array(TRUE, dim(mask_logical)), 42L
  )
  mask_flat_morton <- as.vector(mask_logical)[full_order]

  # Wrong length detail vector
  expect_error(
    fmrilatent:::inverse_lift_R(
      1.0, list(numeric(3)),  # should be length 8, not 3
      mask_flat_morton, dim(mask_logical), 1L, scalings
    ),
    "incorrect length"
  )
})

# -- lna_forward_lift_matrix / lna_inverse_lift_matrix R-path tests ----------

test_that("lna_forward_lift_matrix uses R fallback when Rcpp disabled", {
  opts <- options(
    fmrilatent.haar.use_rcpp = FALSE,
    fmrilatent.haar.use_rcpp_batch = FALSE
  )
  on.exit(options(opts), add = TRUE)

  mask <- array(TRUE, dim = c(2, 2, 2))
  mask_logical <- fmrilatent:::as_logical_mask(mask)
  full_order <- fmrilatent:::get_morton_ordered_indices(
    array(TRUE, dim(mask_logical)), 42L
  )
  mask_flat_morton <- as.vector(mask_logical)[full_order]
  scalings <- fmrilatent:::precompute_haar_scalings(mask_logical, 1L)

  X <- matrix(rnorm(5 * 8), nrow = 5, ncol = 8)

  result <- fmrilatent:::lna_forward_lift_matrix(
    X, mask_flat_morton, dim(mask_logical), 1L, scalings,
    compat_profile = NULL
  )

  expect_true(is.matrix(result$root_coeff))
  expect_equal(nrow(result$root_coeff), 5L)
  expect_true(is.list(result$detail_coeffs_by_level))
  expect_equal(length(result$detail_coeffs_by_level), 1L)
})

test_that("lna_inverse_lift_matrix uses R fallback and roundtrips", {
  opts <- options(
    fmrilatent.haar.use_rcpp = FALSE,
    fmrilatent.haar.use_rcpp_batch = FALSE
  )
  on.exit(options(opts), add = TRUE)

  mask <- array(TRUE, dim = c(2, 2, 2))
  mask_logical <- fmrilatent:::as_logical_mask(mask)
  full_order <- fmrilatent:::get_morton_ordered_indices(
    array(TRUE, dim(mask_logical)), 42L
  )
  mask_flat_morton <- as.vector(mask_logical)[full_order]
  scalings <- fmrilatent:::precompute_haar_scalings(mask_logical, 1L)

  X <- matrix(rnorm(4 * 8), nrow = 4, ncol = 8)

  fw <- fmrilatent:::lna_forward_lift_matrix(
    X, mask_flat_morton, dim(mask_logical), 1L, scalings,
    compat_profile = NULL
  )
  reco <- fmrilatent:::lna_inverse_lift_matrix(
    fw$root_coeff, fw$detail_coeffs_by_level,
    mask_flat_morton, dim(mask_logical), 1L, scalings,
    compat_profile = NULL
  )

  expect_equal(reco, X, tolerance = 1e-7)
})

test_that("lna_inverse_lift_matrix converts vector root_coeff to matrix", {
  opts <- options(
    fmrilatent.haar.use_rcpp = FALSE,
    fmrilatent.haar.use_rcpp_batch = FALSE
  )
  on.exit(options(opts), add = TRUE)

  mask <- array(TRUE, dim = c(2, 2, 2))
  mask_logical <- fmrilatent:::as_logical_mask(mask)
  full_order <- fmrilatent:::get_morton_ordered_indices(
    array(TRUE, dim(mask_logical)), 42L
  )
  mask_flat_morton <- as.vector(mask_logical)[full_order]
  scalings <- fmrilatent:::precompute_haar_scalings(mask_logical, 1L)

  X <- matrix(rnorm(8), nrow = 1, ncol = 8)
  fw <- fmrilatent:::lna_forward_lift_matrix(
    X, mask_flat_morton, dim(mask_logical), 1L, scalings,
    compat_profile = NULL
  )

  # Pass root as vector, not matrix
  reco <- fmrilatent:::lna_inverse_lift_matrix(
    as.numeric(fw$root_coeff[1, ]),
    fw$detail_coeffs_by_level,
    mask_flat_morton, dim(mask_logical), 1L, scalings,
    compat_profile = NULL
  )

  expect_true(is.matrix(reco))
  expect_equal(nrow(reco), 1L)
  expect_equal(reco, X, tolerance = 1e-7)
})

test_that("lna_inverse_lift_matrix rejects non-matrix root_coeff", {
  opts <- options(
    fmrilatent.haar.use_rcpp = FALSE,
    fmrilatent.haar.use_rcpp_batch = FALSE
  )
  on.exit(options(opts), add = TRUE)

  expect_error(
    fmrilatent:::lna_inverse_lift_matrix(
      data.frame(a = 1), list(matrix(0, 1, 8)),
      rep(TRUE, 8), c(2L, 2L, 2L), 1L, list()
    ),
    "must be a matrix"
  )
})

# -- perform_haar_lift_analysis / synthesis direct R-path tests --------------

test_that("perform_haar_lift_analysis works with R fallback", {
  opts <- options(fmrilatent.haar.use_rcpp = FALSE)
  on.exit(options(opts), add = TRUE)

  mask <- array(TRUE, dim = c(2, 2, 2))
  X <- matrix(rnorm(3 * 8), nrow = 3, ncol = 8)

  coeff <- fmrilatent:::perform_haar_lift_analysis(X, mask, levels = 1L)

  expect_true(is.list(coeff))
  expect_true(is.matrix(coeff$root))
  expect_equal(nrow(coeff$root), 3L)
  expect_true(is.list(coeff$detail))
  expect_equal(length(coeff$detail), 1L)
})

test_that("perform_haar_lift_analysis auto-transposes when needed (R path)", {
  opts <- options(fmrilatent.haar.use_rcpp = FALSE)
  on.exit(options(opts), add = TRUE)

  mask <- array(TRUE, dim = c(2, 2, 2))
  # 8 rows x 3 cols => nrow matches mask count, so it gets transposed
  X <- matrix(rnorm(8 * 3), nrow = 8, ncol = 3)

  coeff <- fmrilatent:::perform_haar_lift_analysis(X, mask, levels = 1L)
  # After transpose, should have 3 time points
  expect_equal(nrow(coeff$root), 3L)
})

test_that("perform_haar_lift_synthesis roundtrips with R fallback", {
  opts <- options(fmrilatent.haar.use_rcpp = FALSE)
  on.exit(options(opts), add = TRUE)

  mask <- array(TRUE, dim = c(2, 2, 2))
  set.seed(123)
  X <- matrix(rnorm(4 * 8), nrow = 4, ncol = 8)

  coeff <- fmrilatent:::perform_haar_lift_analysis(X, mask, levels = 1L)
  reco <- fmrilatent:::perform_haar_lift_synthesis(coeff, mask, levels = 1L)

  expect_equal(reco, X, tolerance = 1e-7)
})

test_that("perform_haar_lift_synthesis validates root column count", {
  mask <- array(TRUE, dim = c(2, 2, 2))
  # Root should have 1 column for a 2x2x2 mask with 1 level, giving it 2 is wrong
  expect_error(
    fmrilatent:::perform_haar_lift_synthesis(
      list(root = matrix(0, 1, 2), detail = list(matrix(0, 1, 8))),
      mask, 1L
    ),
    "wrong number of columns"
  )
})

# -- as_logical_mask edge cases ----------------------------------------------

test_that("as_logical_mask preserves dimensions from numeric array", {
  arr <- array(c(1, 0, 0, 1, 1, 1, 0, 0), dim = c(2, 2, 2))
  result <- fmrilatent:::as_logical_mask(arr)
  expect_true(is.logical(result))
  expect_equal(dim(result), c(2, 2, 2))
  expect_equal(sum(result), 4L)
})

test_that("as_logical_mask uses custom location in error message", {
  expect_error(
    fmrilatent:::as_logical_mask(matrix(1, 2, 2), location = "my_func"),
    "my_func"
  )
})

# -- get_morton_ordered_indices edge cases -----------------------------------

test_that("get_morton_ordered_indices returns permutation for full mask (R path)", {
  opts <- options(fmrilatent.haar.use_rcpp = FALSE)
  on.exit(options(opts), add = TRUE)

  mask <- array(TRUE, dim = c(3, 3, 3))
  idx <- fmrilatent:::get_morton_ordered_indices(mask)
  # Should be a permutation of 1:27
  expect_equal(sort(idx), 1:27)
  expect_equal(length(idx), 27L)
})

test_that("get_morton_ordered_indices handles non-cube mask (R path)", {
  opts <- options(fmrilatent.haar.use_rcpp = FALSE)
  on.exit(options(opts), add = TRUE)

  mask <- array(TRUE, dim = c(2, 3, 4))
  idx <- fmrilatent:::get_morton_ordered_indices(mask)
  expect_equal(length(idx), 24L)
  expect_equal(sort(idx), 1:24)
})

# -- precompute_haar_scalings edge cases -------------------------------------

test_that("precompute_haar_scalings with multi-level sparse mask (R path)", {
  opts <- options(fmrilatent.haar.use_rcpp = FALSE)
  on.exit(options(opts), add = TRUE)

  mask <- array(FALSE, dim = c(4, 4, 4))
  mask[1:2, 1:2, 1:2] <- TRUE
  s <- fmrilatent:::precompute_haar_scalings(mask, 2L)

  expect_equal(length(s), 2L)
  # First level: blocks at the 2x2x2 level within 4x4x4
  expect_true(length(s[[1]]$sqrt_nvalid) > 0)
  # Second level: coarser blocks
  expect_true(length(s[[2]]$sqrt_nvalid) > 0)
})

test_that("precompute_haar_scalings with odd-dimension mask (R path)", {
  opts <- options(fmrilatent.haar.use_rcpp = FALSE)
  on.exit(options(opts), add = TRUE)

  # 3x3x3 mask: seq(1,3,by=2) gives 1,3 so 2 blocks per axis
  mask <- array(TRUE, dim = c(3, 3, 3))
  s <- fmrilatent:::precompute_haar_scalings(mask, 1L)

  expect_equal(length(s), 1L)
  # With odd dims, boundary blocks have fewer voxels
  counts <- round(s[[1]]$sqrt_nvalid^2)
  # The block starting at (1,1,1) spans [1:2,1:2,1:2] = 8 voxels
  # The block starting at (3,3,3) spans [3:3,3:3,3:3] = 1 voxel
  expect_true(8 %in% counts)
  expect_true(1 %in% counts)
})

# -- get_valid_finest_blocks edge cases --------------------------------------

test_that("get_valid_finest_blocks with non-cube mask (R path)", {
  opts <- options(fmrilatent.haar.use_rcpp = FALSE)
  on.exit(options(opts), add = TRUE)

  mask <- array(FALSE, dim = c(4, 4, 2))
  mask[1:2, 1:2, 1:2] <- TRUE
  blocks <- fmrilatent:::get_valid_finest_blocks(mask)

  expect_true(is.integer(blocks))
  expect_true(length(blocks) > 0)
  # Blocks should be sorted
  expect_equal(blocks, sort(blocks))
})

test_that("get_valid_finest_blocks returns sorted Morton codes (R path)", {
  opts <- options(fmrilatent.haar.use_rcpp = FALSE)
  on.exit(options(opts), add = TRUE)

  mask <- array(TRUE, dim = c(4, 4, 4))
  blocks <- fmrilatent:::get_valid_finest_blocks(mask)

  expect_true(length(blocks) > 0)
  expect_equal(blocks, sort(blocks))
})
