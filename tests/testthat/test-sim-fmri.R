test_that("LatentNeuroVec reconstructs simulated fMRI exactly (small mask)", {
  skip_on_cran()

  # tiny mask to keep test fast/deterministic
  mask_arr <- array(TRUE, dim = c(3, 3, 3))
  mask_vol <- LogicalNeuroVol(mask_arr, NeuroSpace(c(3, 3, 3)))

  n_time <- 8L
  sim <- neuroim2::simulate_fmri(
    mask = mask_vol,
    n_time = n_time,
    spatial_fwhm = 1,   # keep small but positive
    ar_mean = 0,
    ar_sd = 0,
    noise_sd = 0.1,
    n_factors = 0,
    global_amp = 0,
    seed = 123,
    return_centered = FALSE
  )

  sim_arr <- as.array(sim)
  mask_idx <- which(mask_arr)
  n_vox <- length(mask_idx)

  # build time x vox matrix in mask order
  X <- matrix(0, nrow = n_time, ncol = n_vox)
  for (t in seq_len(n_time)) {
    vol_t <- sim_arr[, , , t]
    X[t, ] <- vol_t[mask_idx]
  }

  # exact factorization via SVD: X = U D V^T
  s <- svd(X)
  k <- length(s$d)
  B <- s$u %*% diag(s$d^0.5, nrow = k, ncol = k)
  L <- s$v %*% diag(s$d^0.5, nrow = k, ncol = k)

  spc <- neuroim2::add_dim(space(mask_vol), n_time)

  lvec <- LatentNeuroVec(
    basis = Matrix::Matrix(B, sparse = FALSE),
    loadings = Matrix::Matrix(L, sparse = FALSE),
    space = spc,
    mask = mask_vol,
    offset = numeric(0),
    label = "svd-sim"
  )

  # reconstruct all voxels/time via linear_access
  recon_vec <- linear_access(lvec, seq_len(prod(dim(lvec))))

  # original flattened in the same order (space major then time)
  original_vec <- c(sim_arr)

  expect_equal(recon_vec, original_vec, tolerance = 1e-6)
})

test_that("haar_latent roundtrip matches input on simulated data (small mask)", {
  skip_on_cran()
  set.seed(456)

  mask_arr <- array(TRUE, dim = c(2, 2, 2))
  mask_vol <- LogicalNeuroVol(mask_arr, NeuroSpace(c(2, 2, 2)))
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
    seed = 456,
    return_centered = FALSE
  )

  sim_arr <- as.array(sim)
  mask_idx <- which(mask_arr)
  X <- matrix(0, nrow = n_time, ncol = length(mask_idx))
  for (t in seq_len(n_time)) {
    X[t, ] <- sim_arr[, , , t][mask_idx]
  }

  opts <- options(fmrilatent.haar.use_rcpp = FALSE)
  hl <- haar_latent(X, mask_arr, levels = 1, z_seed = 99L, threshold = list(type = "none", value = 0))
  options(opts)

  recon <- as.matrix(hl)
  expect_equal(recon, X, tolerance = 1e-6)
})
