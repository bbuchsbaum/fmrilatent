test_that("dpss_time_basis dense vs tridiag agree on small N", {
  skip_if_not_installed("Matrix")
  n <- 64
  tr <- 2
  bw <- 0.05
  k <- 5
  B_dense <- dpss_time_basis(n, tr, bw, k = k, backend = "dense")
  B_tri   <- dpss_time_basis(n, tr, bw, k = k, backend = "tridiag")
  # Compare subspace via principal angle (should be ~0)
  proj <- svd(crossprod(B_dense, B_tri))$d
  expect_true(all(abs(proj) > 1 - 1e-6))
  # Orthonormality
  I_est <- crossprod(B_tri)
  expect_equal(I_est, diag(k), tolerance = 1e-6)
})

test_that("slepian_temporal_latent projects correctly on toy data", {
  n_time <- 50
  n_vox <- 10
  set.seed(123)
  X <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  mask_arr <- array(TRUE, dim = c(2, 5, 1))
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(2, 5, 1)))
  lat <- slepian_temporal_latent(X, mask, tr = 1, bandwidth = 0.1, k = 4, backend = "tridiag")
  recon <- as.matrix(lat@basis %*% t(lat@loadings))
  err_proj <- norm(X - recon, type = "F")
  err_zero <- norm(X, type = "F")
  expect_lt(err_proj, err_zero)        # projection reduces energy
  # projection formula matches direct projection
  B <- lat@basis
  proj_direct <- as.matrix(B %*% crossprod(B, X))
  expect_equal(recon, proj_direct, tolerance = 1e-6)
})
