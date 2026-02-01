test_that("slepian_temporal_latent reconstructs exactly when k = n_time", {
  mask_arr <- array(TRUE, dim = c(2, 2, 1))  # 4 vox
  mask_vol <- LogicalNeuroVol(mask_arr, NeuroSpace(c(2, 2, 1)))
  n_time <- 6L
  X <- matrix(rnorm(n_time * sum(mask_arr)), nrow = n_time)
  X <- scale(X, center = TRUE, scale = FALSE)  # make colMeans zero
  attr(X, "scaled:center") <- NULL
  attr(X, "scaled:scale") <- NULL

  lv <- slepian_temporal_latent(
    X = X,
    mask = mask_vol,
    tr = 2,
    bandwidth = 0.2,
    k = n_time,
    denoise = FALSE,
    backend = "tridiag"
  )

  recon <- as.matrix(lv)
  expect_equal(recon, X, tolerance = 1e-8)
})

test_that("slepian_spatial_latent reconstructs exactly with identity reduction", {
  mask_arr <- array(TRUE, dim = c(2, 1, 1))  # 2 voxels
  mask_vol <- LogicalNeuroVol(mask_arr, NeuroSpace(c(2, 1, 1)))
  X <- matrix(rnorm(5 * sum(mask_arr)), nrow = 5)

  map <- seq_len(sum(mask_arr))  # one cluster per voxel
  red <- make_cluster_reduction(mask_vol, map)
  lv <- slepian_spatial_latent(
    X = X,
    mask = mask_vol,
    reduction = red,
    spec = basis_slepian(k = 1),
    k_neighbors = 1
  )

  recon <- as.matrix(lv)
  expect_equal(recon, X, tolerance = 1e-8)
})

test_that("slepian_spatiotemporal_latent reconstructs exactly with full ranks", {
  mask_arr <- array(TRUE, dim = c(2, 1, 1))  # 2 voxels
  mask_vol <- LogicalNeuroVol(mask_arr, NeuroSpace(c(2, 1, 1)))
  n_time <- 4L
  X <- matrix(rnorm(n_time * sum(mask_arr)), nrow = n_time)

  map <- seq_len(sum(mask_arr))
  red <- make_cluster_reduction(mask_vol, map)

  lv <- slepian_spatiotemporal_latent(
    X = X,
    mask = mask_vol,
    tr = 2,
    bandwidth = 0.2,
    k_time = n_time,
    reduction = red,
    k_space = 1,
    k_neighbors = 1
  )

  recon <- predict(lv)
  expect_equal(recon, X, tolerance = 1e-8)
})

test_that("slepian_temporal_latent partial linear_access matches full array", {
  mask_arr <- array(TRUE, dim = c(2, 2, 1))
  mask_vol <- LogicalNeuroVol(mask_arr, NeuroSpace(c(2, 2, 1)))
  n_time <- 5L
  X <- matrix(rnorm(n_time * sum(mask_arr)), nrow = n_time)

  lv <- slepian_temporal_latent(
    X = X,
    mask = mask_vol,
    tr = 2,
    bandwidth = 0.2,
    k = n_time,
    denoise = FALSE,
    backend = "tridiag"
  )

  arr_full <- as.array(lv)
  total_len <- length(arr_full)
  idx_subset <- sample(total_len, size = 8L)

  expect_equal(linear_access(lv, idx_subset), arr_full[idx_subset], tolerance = 1e-10)
})

test_that("slepian_spatiotemporal_latent supports time_idx and ROI partial decode", {
  skip_if_not_installed("RSpectra")
  set.seed(23)
  mask_arr <- array(TRUE, dim = c(2, 2, 1))
  mask_vol <- LogicalNeuroVol(mask_arr, NeuroSpace(c(2, 2, 1)))
  n_time <- 6L
  X <- matrix(rnorm(n_time * sum(mask_arr)), nrow = n_time)

  lv <- slepian_spatiotemporal_latent(
    X = X,
    mask = mask_vol,
    tr = 2,
    bandwidth = 0.15,
    k_time = 4,
    k_space = 2,
    k_neighbors = 2
  )

  roi_mask <- array(FALSE, dim = c(2, 2, 1))
  roi_mask[1, 1, 1] <- TRUE
  roi_mask[2, 2, 1] <- TRUE

  t_sel <- c(2, 5)
  reco_partial <- predict(lv, time_idx = t_sel, roi_mask = roi_mask)

  # reconstruct via full predict then subset
  full <- predict(lv)
  mask_idx <- which(as.logical(roi_mask))
  expected <- full[t_sel, mask_idx, drop = FALSE]

  expect_equal(reco_partial, expected, tolerance = 1e-8)
})
