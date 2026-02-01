test_that("bspline + HRBF spatiotemporal latent reconstructs exactly on tiny mask", {
  mask_arr <- array(TRUE, dim = c(2, 2, 1))
  mask_vol <- LogicalNeuroVol(mask_arr, NeuroSpace(c(2, 2, 1)))
  n_time <- 5L
  X <- matrix(rnorm(n_time * sum(mask_arr)), nrow = n_time)
  # center to make reconstruction easier
  X <- X - matrix(colMeans(X), nrow = n_time, ncol = ncol(X), byrow = TRUE)

  spec <- spec_st(
    time = spec_time_bspline(k = n_time, degree = 3, include_intercept = FALSE, orthonormalize = TRUE),
    space = spec_space_hrbf(params = list(sigma0 = 2, levels = 1L, radius_factor = 2.5, kernel_type = "gaussian", seed = 42L))
  )
  lv <- encode(X, spec, mask = mask_vol, materialize = "matrix")

  recon <- predict(lv)
  expect_equal(recon, X, tolerance = 1e-8)
})

test_that("bspline + HRBF partial decode matches subset", {
  mask_arr <- array(TRUE, dim = c(2, 2, 1))
  mask_vol <- LogicalNeuroVol(mask_arr, NeuroSpace(c(2, 2, 1)))
  n_time <- 6L
  X <- matrix(rnorm(n_time * sum(mask_arr)), nrow = n_time)
  X <- X - matrix(colMeans(X), nrow = n_time, ncol = ncol(X), byrow = TRUE)

  lv <- latent_factory("bspline_hrbf_st", x = X, mask = mask_vol,
                       k_time = n_time, degree = 3,
                       params = list(sigma0 = 2, levels = 1L, radius_factor = 2.5, kernel_type = "gaussian", seed = 123),
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
