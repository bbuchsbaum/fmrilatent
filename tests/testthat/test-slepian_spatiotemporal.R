test_that("slepian spatiotemporal latent reconstructs simulated data", {
  skip_if_not_installed("RSpectra")
  set.seed(7)
  mask <- array(TRUE, dim = c(2, 2, 2))
  mask_vol <- LogicalNeuroVol(mask, NeuroSpace(c(2, 2, 2)))
  n_time <- 6L

  sim <- neuroim2::simulate_fmri(
    mask = mask_vol,
    n_time = n_time,
    spatial_fwhm = 1.0,
    ar_mean = 0,
    ar_sd = 0,
    noise_sd = 0.05,
    n_factors = 0,
    global_amp = 0,
    seed = 7,
    return_centered = FALSE
  )

  sim_arr <- as.array(sim)
  mask_idx <- which(mask)
  X <- matrix(0, nrow = n_time, ncol = length(mask_idx))
  for (t in seq_len(n_time)) {
    X[t, ] <- sim_arr[, , , t][mask_idx]
  }

  lv <- slepian_spatiotemporal_latent(
    X = X,
    mask = mask_vol,
    tr = 2,
    bandwidth = 0.1,
    k_time = 3,
    k_space = 2,
    k_neighbors = 3
  )

  reco <- predict(lv)
  # expected from stored factors
  expected <- lv$coeff$B_t %*% lv$coeff$core %*% t(lv$coeff$L_s)
  expect_equal(reco, expected, tolerance = 1e-8)
})
