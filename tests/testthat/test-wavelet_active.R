test_that("wavelet_active_latent roundtrip with time + space levels", {
  mask_arr <- array(TRUE, dim = c(3, 2, 2))
  mask_vol <- LogicalNeuroVol(mask_arr, NeuroSpace(c(3, 2, 2)))
  n_time <- 4L
  X <- matrix(rnorm(n_time * sum(mask_arr)), nrow = n_time)
  lv <- wavelet_active_latent(X, mask_vol, levels_space = 1L, levels_time = 1L, threshold = 0)
  reco <- predict(lv)
  expect_equal(reco, X, tolerance = 1e-8)
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
