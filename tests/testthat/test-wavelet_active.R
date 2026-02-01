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
