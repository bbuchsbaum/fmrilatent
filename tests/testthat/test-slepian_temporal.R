# Tests for slepian_temporal.R
# Testing dpss_time_basis (already tested) and slepian_temporal_latent (new tests)

# =============================================================================
# slepian_temporal_latent tests
# =============================================================================

test_that("slepian_temporal_latent basic construction with small mask", {
  mask_arr <- array(TRUE, dim = c(2, 2, 2))
  mask_vol <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(2, 2, 2)))
  n_time <- 20L
  n_vox <- sum(mask_arr)

  set.seed(42)
  X <- matrix(rnorm(n_time * n_vox), nrow = n_time, ncol = n_vox)

  # Create slepian temporal latent
  lv <- slepian_temporal_latent(X, mask_vol, tr = 2.0, bandwidth = 0.08, k = 4L)

  expect_s4_class(lv, "LatentNeuroVec")
  expect_true(validObject(lv))
})

test_that("slepian_temporal_latent returns LatentNeuroVec", {
  mask_arr <- array(TRUE, dim = c(3, 3, 1))
  mask_vol <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(3, 3, 1)))
  n_time <- 16L
  n_vox <- sum(mask_arr)

  set.seed(123)
  X <- matrix(rnorm(n_time * n_vox), nrow = n_time, ncol = n_vox)

  lv <- slepian_temporal_latent(X, mask_vol, tr = 2.0, bandwidth = 0.1, k = 3L)

  expect_s4_class(lv, "LatentNeuroVec")
})

test_that("slepian_temporal_latent meta contains correct family", {
  mask_arr <- array(TRUE, dim = c(2, 2, 1))
  mask_vol <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(2, 2, 1)))
  n_time <- 24L
  n_vox <- sum(mask_arr)

  set.seed(456)
  X <- matrix(rnorm(n_time * n_vox), nrow = n_time, ncol = n_vox)

  lv <- slepian_temporal_latent(X, mask_vol, tr = 2.0, bandwidth = 0.05, k = 4L)

  expect_equal(lv@meta$family, "slepian_temporal")
  expect_equal(lv@meta$tr, 2.0)
  expect_equal(lv@meta$bandwidth, 0.05)
})

test_that("slepian_temporal_latent has correct dimensions", {
  mask_arr <- array(TRUE, dim = c(3, 3, 2))
  mask_vol <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(3, 3, 2)))
  n_time <- 32L
  n_vox <- sum(mask_arr)
  k <- 5L

  set.seed(789)
  X <- matrix(rnorm(n_time * n_vox), nrow = n_time, ncol = n_vox)

  lv <- slepian_temporal_latent(X, mask_vol, tr = 2.0, bandwidth = 0.08, k = k)

  # Check basis dimensions
  B <- basis(lv)
  expect_equal(nrow(B), n_time)
  expect_equal(ncol(B), k)

  # Check loadings dimensions
  L <- loadings(lv)
  expect_equal(nrow(L), n_vox)
  expect_equal(ncol(L), k)
})

test_that("slepian_temporal_latent with automatic k selection", {
  mask_arr <- array(TRUE, dim = c(2, 2, 1))
  mask_vol <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(2, 2, 1)))
  n_time <- 64L
  n_vox <- sum(mask_arr)
  tr <- 2.0
  bandwidth <- 0.05

  set.seed(111)
  X <- matrix(rnorm(n_time * n_vox), nrow = n_time, ncol = n_vox)

  # k = NULL should auto-select based on Shannon number
  lv <- slepian_temporal_latent(X, mask_vol, tr = tr, bandwidth = bandwidth, k = NULL)

  expect_s4_class(lv, "LatentNeuroVec")

  # Check that k was auto-selected
  B <- basis(lv)
  W <- bandwidth * tr
  NW <- n_time * W
  expected_k <- floor(2 * NW) - 1
  expect_equal(ncol(B), max(1L, min(as.integer(expected_k), n_time)))
})

test_that("slepian_temporal_latent includes offset", {
  mask_arr <- array(TRUE, dim = c(2, 2, 1))
  mask_vol <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(2, 2, 1)))
  n_time <- 20L
  n_vox <- sum(mask_arr)

  set.seed(222)
  X <- matrix(rnorm(n_time * n_vox), nrow = n_time, ncol = n_vox)

  lv <- slepian_temporal_latent(X, mask_vol, tr = 2.0, bandwidth = 0.1, k = 3L)

  # Offset should be column means
  off <- offset(lv)
  expect_length(off, n_vox)
  expect_equal(off, colMeans(X), tolerance = 1e-10)
})

test_that("slepian_temporal_latent with different backend options", {
  mask_arr <- array(TRUE, dim = c(2, 2, 1))
  mask_vol <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(2, 2, 1)))
  n_time <- 16L
  n_vox <- sum(mask_arr)

  set.seed(333)
  X <- matrix(rnorm(n_time * n_vox), nrow = n_time, ncol = n_vox)

  # Test tridiag backend
  lv_tri <- slepian_temporal_latent(X, mask_vol, tr = 2.0, bandwidth = 0.08,
                                     k = 3L, backend = "tridiag")
  expect_s4_class(lv_tri, "LatentNeuroVec")

  # Test dense backend
  lv_dense <- slepian_temporal_latent(X, mask_vol, tr = 2.0, bandwidth = 0.08,
                                       k = 3L, backend = "dense")
  expect_s4_class(lv_dense, "LatentNeuroVec")

  # Results should be very similar (up to numerical precision and sign flips)
  B_tri <- as.matrix(basis(lv_tri))
  B_dense <- as.matrix(basis(lv_dense))

  # Check orthonormality of both
  G_tri <- crossprod(B_tri)
  G_dense <- crossprod(B_dense)
  expect_equal(G_tri, diag(3), tolerance = 1e-6)
  expect_equal(G_dense, diag(3), tolerance = 1e-6)
})

test_that("slepian_temporal_latent with custom label", {
  mask_arr <- array(TRUE, dim = c(2, 2, 1))
  mask_vol <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(2, 2, 1)))
  n_time <- 20L
  n_vox <- sum(mask_arr)

  set.seed(444)
  X <- matrix(rnorm(n_time * n_vox), nrow = n_time, ncol = n_vox)

  lv <- slepian_temporal_latent(X, mask_vol, tr = 2.0, bandwidth = 0.1,
                                 k = 3L, label = "my_slepian")

  expect_equal(lv@label, "my_slepian")
})

test_that("slepian_temporal_latent validates input", {
  mask_arr <- array(TRUE, dim = c(2, 2, 1))
  mask_vol <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(2, 2, 1)))

  # Missing n_time
  expect_error(
    slepian_temporal_latent(matrix(numeric(0), 0, 4), mask_vol, tr = 2.0, bandwidth = 0.1),
    "must have time in rows"
  )
})

test_that("slepian_temporal_latent basis is orthonormal", {
  mask_arr <- array(TRUE, dim = c(3, 3, 1))
  mask_vol <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(3, 3, 1)))
  n_time <- 48L
  n_vox <- sum(mask_arr)
  k <- 6L

  set.seed(555)
  X <- matrix(rnorm(n_time * n_vox), nrow = n_time, ncol = n_vox)

  lv <- slepian_temporal_latent(X, mask_vol, tr = 2.0, bandwidth = 0.06, k = k)

  B <- as.matrix(basis(lv))
  G <- crossprod(B)

  # Should be approximately identity
  expect_equal(G, diag(k), tolerance = 1e-6)
})

test_that("slepian_temporal_latent stores NW in meta", {
  mask_arr <- array(TRUE, dim = c(2, 2, 1))
  mask_vol <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(2, 2, 1)))
  n_time <- 30L
  n_vox <- sum(mask_arr)
  tr <- 2.0
  bandwidth <- 0.1

  set.seed(666)
  X <- matrix(rnorm(n_time * n_vox), nrow = n_time, ncol = n_vox)

  lv <- slepian_temporal_latent(X, mask_vol, tr = tr, bandwidth = bandwidth, k = 4L)

  W <- bandwidth * tr
  expected_NW <- n_time * W

  expect_equal(lv@meta$nw, expected_NW, tolerance = 1e-10)
})

test_that("slepian_temporal_latent with very small bandwidth", {
  mask_arr <- array(TRUE, dim = c(2, 2, 1))
  mask_vol <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(2, 2, 1)))
  n_time <- 40L
  n_vox <- sum(mask_arr)

  set.seed(777)
  X <- matrix(rnorm(n_time * n_vox), nrow = n_time, ncol = n_vox)

  # Very small bandwidth -> small NW -> small k
  lv <- slepian_temporal_latent(X, mask_vol, tr = 2.0, bandwidth = 0.01, k = NULL)

  expect_s4_class(lv, "LatentNeuroVec")

  # k should be at least 1
  B <- basis(lv)
  expect_true(ncol(B) >= 1)
})

test_that("slepian_temporal_latent reconstruction quality", {
  mask_arr <- array(TRUE, dim = c(2, 2, 1))
  mask_vol <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(2, 2, 1)))
  n_time <- 32L
  n_vox <- sum(mask_arr)

  set.seed(888)
  X <- matrix(rnorm(n_time * n_vox), nrow = n_time, ncol = n_vox)
  X_centered <- X - matrix(colMeans(X), nrow = n_time, ncol = n_vox, byrow = TRUE)

  # Use high k for good reconstruction
  lv <- slepian_temporal_latent(X, mask_vol, tr = 2.0, bandwidth = 0.15, k = 15L)

  # Reconstruct
  B <- as.matrix(basis(lv))
  L <- as.matrix(loadings(lv))
  X_recon <- B %*% t(L)

  # Add back offset
  X_recon <- X_recon + matrix(offset(lv), nrow = n_time, ncol = n_vox, byrow = TRUE)

  # Check reconstruction error
  err <- mean((X - X_recon)^2)

  # With k=15 out of 32, reconstruction captures ~47% of variance for random data
  # Error should be less than original data variance
  orig_var <- mean(X^2)
  expect_true(err < orig_var)
})

# =============================================================================
# Additional dpss_time_basis edge case tests
# =============================================================================

test_that("dpss_time_basis validates missing n_time", {
  expect_error(dpss_time_basis(tr = 2, bandwidth = 0.1), "n_time must be positive")
})

test_that("dpss_time_basis validates n_time < 1", {
  expect_error(dpss_time_basis(n_time = 0, tr = 2, bandwidth = 0.1), "n_time must be positive")
  expect_error(dpss_time_basis(n_time = -5, tr = 2, bandwidth = 0.1), "n_time must be positive")
})

test_that("dpss_time_basis validates missing tr", {
  expect_error(dpss_time_basis(n_time = 10, bandwidth = 0.1), "tr must be positive")
})

test_that("dpss_time_basis validates tr <= 0", {
  expect_error(dpss_time_basis(n_time = 10, tr = 0, bandwidth = 0.1), "tr must be positive")
  expect_error(dpss_time_basis(n_time = 10, tr = -1, bandwidth = 0.1), "tr must be positive")
})

test_that("dpss_time_basis validates missing bandwidth", {
  expect_error(dpss_time_basis(n_time = 10, tr = 2), "bandwidth must be positive")
})

test_that("dpss_time_basis validates bandwidth <= 0", {
  expect_error(dpss_time_basis(n_time = 10, tr = 2, bandwidth = 0), "bandwidth must be positive")
  expect_error(dpss_time_basis(n_time = 10, tr = 2, bandwidth = -0.05), "bandwidth must be positive")
})

test_that("dpss_time_basis with n_time = 1 returns a single-row matrix", {
  B <- dpss_time_basis(n_time = 1, tr = 2, bandwidth = 0.1, k = 1)
  expect_true(is.matrix(B))
  expect_equal(nrow(B), 1L)
  expect_equal(ncol(B), 1L)
})

test_that("dpss_time_basis clamps k to [1, n_time]", {
  # k = 0 should be clamped to 1
  B <- dpss_time_basis(n_time = 10, tr = 2, bandwidth = 0.01, k = 0)
  expect_equal(ncol(B), 1L)

  # k larger than n_time should be clamped
  B2 <- dpss_time_basis(n_time = 5, tr = 2, bandwidth = 0.1, k = 100)
  expect_equal(ncol(B2), 5L)
})

test_that("dpss_time_basis auto-selects k from NW when k is NULL", {
  n_time <- 50
  tr <- 2.0
  bandwidth <- 0.05
  W <- bandwidth * tr
  NW <- n_time * W
  expected_k <- max(1L, min(as.integer(floor(2 * NW) - 1), n_time))

  B <- dpss_time_basis(n_time = n_time, tr = tr, bandwidth = bandwidth, k = NULL)
  expect_equal(ncol(B), expected_k)
})

test_that("dpss_time_basis dense backend returns orthonormal matrix", {
  B <- dpss_time_basis(n_time = 20, tr = 2, bandwidth = 0.08, k = 3, backend = "dense")
  expect_true(is.matrix(B))
  expect_equal(nrow(B), 20L)
  expect_equal(ncol(B), 3L)
  G <- crossprod(B)
  expect_equal(G, diag(3), tolerance = 1e-6)
})

test_that("dpss_time_basis tridiag backend returns orthonormal matrix", {
  B <- dpss_time_basis(n_time = 20, tr = 2, bandwidth = 0.08, k = 3, backend = "tridiag")
  expect_true(is.matrix(B))
  G <- crossprod(B)
  expect_equal(G, diag(3), tolerance = 1e-6)
})

test_that("dpss_time_basis rejects invalid backend", {
  expect_error(
    dpss_time_basis(n_time = 10, tr = 2, bandwidth = 0.1, backend = "invalid"),
    "should be one of"
  )
})

test_that("dpss_time_basis with very small NW yields k = 1", {
  # NW = n_time * bandwidth * tr = 10 * 0.001 * 1 = 0.01

  # floor(2 * 0.01) - 1 = floor(0.02) - 1 = 0 - 1 = -1 => clamped to 1
  B <- dpss_time_basis(n_time = 10, tr = 1, bandwidth = 0.001, k = NULL)
  expect_equal(ncol(B), 1L)
})
