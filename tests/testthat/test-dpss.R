test_that("dpss_time_basis tridiag returns an orthonormal basis", {
  n <- 64
  tr <- 2
  bw <- 0.05
  k <- 5
  B_tri   <- dpss_time_basis(n, tr, bw, k = k, backend = "tridiag")
  I_est <- crossprod(B_tri)
  expect_equal(I_est, diag(k), tolerance = 1e-6)
})

test_that("dpss_time_basis rejects dense backend", {
  expect_error(
    dpss_time_basis(64, tr = 2, bandwidth = 0.05, k = 5, backend = "dense"),
    class = "fmrilatent_error_unsupported_dpss_backend"
  )
})

test_that("slepian_temporal_latent projects correctly on toy data", {
  n_time <- 50
  n_vox <- 10
  set.seed(123)
  X <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  mask_arr <- array(TRUE, dim = c(2, 5, 1))
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(2, 5, 1)))
  lat <- slepian_temporal_latent(X, mask, tr = 1, bandwidth = 0.1, k = 4, backend = "tridiag")
  recon_centered <- as.matrix(lat@basis %*% t(lat@loadings))
  X_centered <- sweep(X, 2L, colMeans(X), "-")
  err_proj <- norm(X_centered - recon_centered, type = "F")
  err_zero <- norm(X_centered, type = "F")
  expect_lt(err_proj, err_zero)        # projection reduces energy
  # projection formula matches direct projection
  B <- lat@basis
  proj_direct <- as.matrix(B %*% crossprod(B, X_centered))
  expect_equal(recon_centered, proj_direct, tolerance = 1e-6)
})
