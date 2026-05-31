test_that("wavelet_active_latent roundtrip with time + space levels", {
  mask_arr <- array(TRUE, dim = c(3, 2, 2))
  mask_vol <- LogicalNeuroVol(mask_arr, NeuroSpace(c(3, 2, 2)))
  n_time <- 4L
  X <- matrix(rnorm(n_time * sum(mask_arr)), nrow = n_time)
  lv <- wavelet_active_latent(X, mask_vol, levels_space = 1L, levels_time = 1L, threshold = 0)
  reco <- predict(lv)
  expect_equal(reco, X, tolerance = 1e-8)
})

test_that("wavelet_active_latent roundtrips sparse non-contiguous masks", {
  mask_arr <- array(FALSE, dim = c(3, 3, 3))
  mask_arr[cbind(
    c(1L, 2L, 3L, 1L, 3L),
    c(1L, 2L, 3L, 3L, 1L),
    c(1L, 2L, 1L, 2L, 3L)
  )] <- TRUE
  mask_vol <- LogicalNeuroVol(mask_arr, NeuroSpace(c(3, 3, 3)))
  n_time <- 4L
  set.seed(731)
  X <- matrix(rnorm(n_time * sum(mask_arr)), nrow = n_time)

  lv <- wavelet_active_latent(
    X,
    mask_vol,
    levels_space = 1L,
    levels_time = 0L,
    threshold = 0
  )
  reco <- predict(lv)

  expect_equal(reco, X, tolerance = 1e-8)
})

test_that("wavelet_active_latent preserves time rows for a single active voxel", {
  mask_arr <- array(FALSE, dim = c(2, 2, 1))
  mask_arr[1, 1, 1] <- TRUE
  mask_vol <- LogicalNeuroVol(mask_arr, NeuroSpace(c(2, 2, 1)))
  X <- matrix(seq_len(5), nrow = 5, ncol = 1)

  lv <- wavelet_active_latent(
    X,
    mask_vol,
    levels_space = 1L,
    levels_time = 0L,
    threshold = 0
  )
  reco <- predict(lv)

  expect_equal(lv$meta$n_time, nrow(X))
  expect_equal(dim(lv$coeff), dim(X))
  expect_equal(reco, X, tolerance = 1e-8)
})

test_that("wavelet_active_latent validates degenerate inputs early", {
  empty_mask <- array(FALSE, dim = c(2, 2, 1))
  expect_error(
    wavelet_active_latent(matrix(0, nrow = 1L, ncol = 0L), empty_mask),
    class = "fmrilatent_error_empty_mask"
  )

  mask_arr <- array(TRUE, dim = c(2, 2, 1))
  expect_error(
    wavelet_active_latent(matrix(1, nrow = 2L, ncol = 4L), mask_arr, levels_space = 0L),
    class = "fmrilatent_error_invalid_count"
  )
  expect_error(
    wavelet_active_latent(matrix(1, nrow = 0L, ncol = 4L), mask_arr, levels_space = 1L),
    class = "fmrilatent_error_invalid_count"
  )
  expect_error(
    wavelet_active_latent(matrix(1, nrow = 2L, ncol = 3L), mask_arr, levels_space = 1L),
    class = "fmrilatent_error_dim"
  )
  expect_error(
    wavelet_active_latent(array(1, dim = c(3L, 2L, 1L, 2L)), mask_arr, levels_space = 1L),
    class = "fmrilatent_error_dim"
  )
  bad_x <- matrix(1, nrow = 2L, ncol = 4L)
  bad_x[1, 1] <- Inf
  expect_error(
    wavelet_active_latent(bad_x, mask_arr, levels_space = 1L),
    class = "fmrilatent_error_value"
  )
})

test_that("wavelet_active_latent partial decode respects time_idx and ROI", {
  mask_arr <- array(TRUE, dim = c(2, 2, 1))
  mask_vol <- LogicalNeuroVol(mask_arr, NeuroSpace(c(2, 2, 1)))
  n_time <- 5L
  X <- matrix(rnorm(n_time * sum(mask_arr)), nrow = n_time)
  lv <- wavelet_active_latent(X, mask_vol, levels_space = 1L, levels_time = 1L, threshold = 0)

  roi_mask <- array(FALSE, dim = c(2, 2, 1))
  roi_mask[1, 1, 1] <- TRUE
  roi_mask[2, 2, 1] <- TRUE
  t_sel <- c(2, 5)

  partial <- predict(lv, time_idx = t_sel, roi_mask = roi_mask)
  full <- predict(lv)
  idx <- which(as.logical(roi_mask))
  expect_equal(partial, full[t_sel, idx, drop = FALSE], tolerance = 1e-8)
})

test_that("wavelet_active_latent validates decoder roi_mask", {
  mask_arr <- array(TRUE, dim = c(2, 2, 1))
  mask_arr[2, 1, 1] <- FALSE
  mask_vol <- LogicalNeuroVol(mask_arr, NeuroSpace(c(2, 2, 1)))
  X <- matrix(rnorm(4L * sum(mask_arr)), nrow = 4L)
  lv <- wavelet_active_latent(X, mask_vol, levels_space = 1L, levels_time = 0L, threshold = 0)

  bad_dims <- array(TRUE, dim = c(3, 2, 1))
  expect_error(predict(lv, roi_mask = bad_dims), class = "fmrilatent_error_dim")

  outside_mask <- array(FALSE, dim = dim(mask_arr))
  outside_mask[2, 1, 1] <- TRUE
  expect_error(predict(lv, roi_mask = outside_mask), class = "fmrilatent_error_value")
})

test_that("bspline_hrbf_st factory roundtrip matches original data", {
  mask_arr <- array(TRUE, dim = c(2, 2, 1))
  mask_vol <- LogicalNeuroVol(mask_arr, NeuroSpace(c(2, 2, 1)))
  n_time <- 5L
  X <- matrix(rnorm(n_time * sum(mask_arr)), nrow = n_time)
  X <- X - matrix(colMeans(X), nrow = n_time, ncol = ncol(X), byrow = TRUE)

  lv <- latent_factory("bspline_hrbf_st", x = X, mask = mask_vol,
                       k_time = n_time, degree = 3,
                       params = list(sigma0 = 2, levels = 1L, radius_factor = 2.5, kernel_type = "gaussian", seed = 111),
                       materialize = "matrix")
  reco <- predict(lv)
  expect_equal(reco, X, tolerance = 1e-8)
})

test_that("bspline_hrbf_st partial decode (time_idx + ROI) matches subset", {
  mask_arr <- array(TRUE, dim = c(2, 2, 1))
  mask_vol <- LogicalNeuroVol(mask_arr, NeuroSpace(c(2, 2, 1)))
  n_time <- 6L
  X <- matrix(rnorm(n_time * sum(mask_arr)), nrow = n_time)
  X <- X - matrix(colMeans(X), nrow = n_time, ncol = ncol(X), byrow = TRUE)

  lv <- latent_factory("bspline_hrbf_st", x = X, mask = mask_vol,
                       k_time = n_time, degree = 3,
                       params = list(sigma0 = 2, levels = 1L, radius_factor = 2.5, kernel_type = "gaussian", seed = 222),
                       materialize = "matrix")

  t_sel <- c(2, 5)
  roi_mask <- array(FALSE, dim = c(2, 2, 1))
  roi_mask[1, 1, 1] <- TRUE
  roi_mask[2, 2, 1] <- TRUE

  partial <- predict(lv, time_idx = t_sel, roi_mask = roi_mask)
  full <- predict(lv)
  idx <- which(as.logical(roi_mask))
  expect_equal(partial, full[t_sel, idx, drop = FALSE], tolerance = 1e-8)
})

# =============================================================================
# Additional tests for wavelet_active_latent
# =============================================================================

test_that("wavelet_active_latent basic roundtrip with small mask", {
  mask_arr <- array(TRUE, dim = c(2, 2, 2))
  mask_vol <- LogicalNeuroVol(mask_arr, NeuroSpace(c(2, 2, 2)))
  n_time <- 8L
  set.seed(42)
  X <- matrix(rnorm(n_time * sum(mask_arr)), nrow = n_time)

  lv <- wavelet_active_latent(X, mask_vol, levels_space = 1L, levels_time = 0L, threshold = 0)
  reco <- predict(lv)

  expect_equal(reco, X, tolerance = 1e-8)
  expect_s3_class(lv, "ImplicitLatent")
})

test_that("wavelet_active_latent with time lifting (levels_time > 0)", {
  mask_arr <- array(TRUE, dim = c(2, 2, 1))
  mask_vol <- LogicalNeuroVol(mask_arr, NeuroSpace(c(2, 2, 1)))
  n_time <- 8L
  set.seed(123)
  X <- matrix(rnorm(n_time * sum(mask_arr)), nrow = n_time)

  lv <- wavelet_active_latent(X, mask_vol, levels_space = 1L, levels_time = 2L, threshold = 0)
  reco <- predict(lv)

  expect_equal(reco, X, tolerance = 1e-8)
})

test_that("wavelet_active_latent without time lifting (levels_time = 0)", {
  mask_arr <- array(TRUE, dim = c(2, 2, 1))
  mask_vol <- LogicalNeuroVol(mask_arr, NeuroSpace(c(2, 2, 1)))
  n_time <- 4L
  set.seed(456)
  X <- matrix(rnorm(n_time * sum(mask_arr)), nrow = n_time)

  lv <- wavelet_active_latent(X, mask_vol, levels_space = 1L, levels_time = 0L, threshold = 0)
  reco <- predict(lv)

  expect_equal(reco, X, tolerance = 1e-8)
})

test_that("wavelet_active_latent with thresholding", {
  mask_arr <- array(TRUE, dim = c(2, 2, 1))
  mask_vol <- LogicalNeuroVol(mask_arr, NeuroSpace(c(2, 2, 1)))
  n_time <- 8L
  set.seed(789)
  X <- matrix(rnorm(n_time * sum(mask_arr)), nrow = n_time)

  # With threshold, reconstruction won't be exact but should be close
  lv <- wavelet_active_latent(X, mask_vol, levels_space = 1L, levels_time = 1L, threshold = 0.1)
  reco <- predict(lv)

  # Check dimensions match
  expect_equal(dim(reco), dim(X))

  # Reconstruction should be approximate (lossy due to threshold)
  # Error should be bounded
  err <- mean((reco - X)^2)
  expect_true(err < 1)  # MSE less than 1
})

test_that("wavelet_active_latent accepts 4D array input", {
  mask_arr <- array(TRUE, dim = c(2, 2, 2))
  mask_vol <- LogicalNeuroVol(mask_arr, NeuroSpace(c(2, 2, 2)))
  n_time <- 4L
  set.seed(111)

  # Create 4D array
  X_4d <- array(rnorm(2 * 2 * 2 * n_time), dim = c(2, 2, 2, n_time))

  lv <- wavelet_active_latent(X_4d, mask_vol, levels_space = 1L, levels_time = 0L, threshold = 0)
  reco <- predict(lv)

  # Extract expected matrix from 4D array
  idx <- which(mask_arr)
  X_mat <- matrix(0, nrow = n_time, ncol = length(idx))
  for (t in seq_len(n_time)) {
    X_mat[t, ] <- X_4d[, , , t][idx]
  }

  expect_equal(reco, X_mat, tolerance = 1e-8)
})

test_that("wavelet_active_latent ROI mask in decoder", {
  mask_arr <- array(TRUE, dim = c(4, 4, 1))
  mask_vol <- LogicalNeuroVol(mask_arr, NeuroSpace(c(4, 4, 1)))
  n_time <- 8L
  set.seed(222)
  X <- matrix(rnorm(n_time * sum(mask_arr)), nrow = n_time)

  lv <- wavelet_active_latent(X, mask_vol, levels_space = 1L, levels_time = 1L, threshold = 0)

  # Create ROI that is subset of mask
  roi_mask <- array(FALSE, dim = c(4, 4, 1))
  roi_mask[1:2, 1:2, 1] <- TRUE  # Top-left quadrant

  reco_roi <- predict(lv, roi_mask = roi_mask)
  reco_full <- predict(lv)

  # ROI reconstruction should match subset of full
  idx <- which(as.logical(roi_mask))
  expect_equal(reco_roi, reco_full[, idx, drop = FALSE], tolerance = 1e-8)
})

test_that("wavelet_active_latent time_idx subsetting", {
  mask_arr <- array(TRUE, dim = c(2, 2, 1))
  mask_vol <- LogicalNeuroVol(mask_arr, NeuroSpace(c(2, 2, 1)))
  n_time <- 10L
  set.seed(333)
  X <- matrix(rnorm(n_time * sum(mask_arr)), nrow = n_time)

  lv <- wavelet_active_latent(X, mask_vol, levels_space = 1L, levels_time = 2L, threshold = 0)

  # Select subset of time points
  t_sel <- c(2, 5, 8)
  reco_subset <- predict(lv, time_idx = t_sel)
  reco_full <- predict(lv)

  expect_equal(reco_subset, reco_full[t_sel, , drop = FALSE], tolerance = 1e-8)
})

test_that("wavelet_active_latent returns ImplicitLatent object", {
  mask_arr <- array(TRUE, dim = c(2, 2, 1))
  mask_vol <- LogicalNeuroVol(mask_arr, NeuroSpace(c(2, 2, 1)))
  n_time <- 4L
  X <- matrix(rnorm(n_time * sum(mask_arr)), nrow = n_time)

  lv <- wavelet_active_latent(X, mask_vol, levels_space = 1L, levels_time = 0L)

  expect_s3_class(lv, "ImplicitLatent")
})

test_that("wavelet_active_latent meta contains correct family", {
  mask_arr <- array(TRUE, dim = c(2, 2, 1))
  mask_vol <- LogicalNeuroVol(mask_arr, NeuroSpace(c(2, 2, 1)))
  n_time <- 4L
  X <- matrix(rnorm(n_time * sum(mask_arr)), nrow = n_time)

  lv <- wavelet_active_latent(X, mask_vol, levels_space = 1L, levels_time = 0L)

  expect_equal(lv$meta$family, "wavelet_active")
})

test_that("wavelet_active_latent keeps LogicalNeuroVol geometry for array reconstruction", {
  mask_arr <- array(TRUE, dim = c(2, 2, 1))
  spacing <- c(2, 3, 4)
  origin <- c(10, 20, 30)
  mask_vol <- neuroim2::LogicalNeuroVol(
    mask_arr,
    neuroim2::NeuroSpace(dim(mask_arr), spacing = spacing, origin = origin)
  )
  X <- matrix(rnorm(4L * sum(mask_arr)), nrow = 4L)

  lv <- wavelet_active_latent(X, mask_vol, levels_space = 1L, levels_time = 0L)
  rec <- reconstruct_array(lv)

  expect_s4_class(rec, "DenseNeuroVec")
  expect_equal(neuroim2::spacing(neuroim2::space(rec)), spacing)
  expect_equal(neuroim2::origin(neuroim2::space(rec)), origin)
})
