test_that("hrbf reconstruction error declines as atom count increases", {
  mask_arr <- array(TRUE, dim = c(4, 4, 4))
  mask_vol <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(dim(mask_arr)))
  n_vox <- sum(mask_arr)
  n_time <- 6L
  X <- matrix(rnorm(n_time * n_vox), nrow = n_time)

  params_low <- list(sigma0 = 2, levels = 0L, radius_factor = 2.5, kernel_type = "gaussian", seed = 1L)
  params_high <- list(sigma0 = 2, levels = 1L, radius_factor = 2.5, kernel_type = "gaussian", seed = 1L)

  lv_low <- encode(X, spec_space_hrbf(params_low), mask = mask_vol, materialize = "matrix")
  lv_high <- encode(X, spec_space_hrbf(params_high), mask = mask_vol, materialize = "matrix")

  rec_low <- as.matrix(lv_low)
  rec_high <- as.matrix(lv_high)

  rmse <- function(a, b) sqrt(mean((a - b)^2))
  err_low <- rmse(rec_low, X)
  err_high <- rmse(rec_high, X)

  expect_lt(err_high, err_low + 1e-8)
})
