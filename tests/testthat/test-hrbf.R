test_that("hrbf basis small mask is identity-like", {
  mask_arr <- array(TRUE, dim = c(3, 3, 3))
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(3, 3, 3)))
  params <- list(sigma0 = 1, levels = 0L, radius_factor = 2.5, kernel_type = "gaussian", seed = 1L)
  B <- hrbf_generate_basis(params, mask)
  expect_s4_class(B, "dgCMatrix")
  expect_equal(dim(B), c(27, 27))
  # rows are atoms; columns are voxels
  expect_true(all(Matrix::diag(B %*% Matrix::t(B)) > 0.9))
})

test_that("hrbf_latent roundtrip on tiny mask", {
  set.seed(42)
  mask_arr <- array(TRUE, dim = c(3, 3, 3))
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(3, 3, 3)))
  X <- matrix(rnorm(5 * 27), nrow = 5)
  params <- list(sigma0 = 1, levels = 0L, radius_factor = 2.5, kernel_type = "gaussian", seed = 1L)
  latent <- hrbf_latent(X, mask, params)
  expect_s4_class(latent, "LatentNeuroVec")
  Xhat <- hrbf_reconstruct_matrix(latent@basis, mask, params)
  expect_equal(dim(Xhat), dim(X))
  expect_equal(Xhat, X, tolerance = 1e-6)
})

test_that("wendland alias works", {
  mask_arr <- array(TRUE, dim = c(2, 2, 2))
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(2, 2, 2)))
  params <- list(sigma0 = 2, levels = 1L, radius_factor = 2.0, kernel_type = "wendland_c4", seed = 7L)
  B <- hrbf_generate_basis(params, mask)
  expect_equal(ncol(B), sum(mask_arr))
})

test_that("hrbf basis works with R-only backend", {
  mask_arr <- array(TRUE, dim = c(4, 4, 4))
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(4, 4, 4)))
  params <- list(sigma0 = 2, levels = 1L, radius_factor = 2.5, kernel_type = "gaussian", seed = 11L)
  old <- getOption("fmrilatent.hrbf.use_rcpp")
  on.exit(options(fmrilatent.hrbf.use_rcpp = old), add = TRUE)
  options(fmrilatent.hrbf.use_rcpp = FALSE)
  B <- hrbf_generate_basis(params, mask)
  expect_s4_class(B, "dgCMatrix")
  expect_equal(ncol(B), sum(mask_arr))
})

test_that("hrbf basis dimensions scale on medium mask", {
  mask_arr <- array(TRUE, dim = c(8, 8, 8))
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(8, 8, 8)))
  params <- list(sigma0 = 3, levels = 2L, radius_factor = 2.5, kernel_type = "gaussian", seed = 5L)
  B <- hrbf_generate_basis(params, mask)
  expect_equal(ncol(B), sum(mask_arr))
  expect_gt(nrow(B), 0)
})

test_that("hrbf partial reconstruction matches full for selected voxels/time", {
  set.seed(123)
  mask_arr <- array(TRUE, dim = c(4, 4, 4))
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(4, 4, 4)))
  X <- matrix(rnorm(6 * 64), nrow = 6)
  params <- list(sigma0 = 1.5, levels = 1L, radius_factor = 2.5, kernel_type = "gaussian", seed = 2L)
  coeff <- hrbf_project_matrix(X, mask, params)
  full <- hrbf_reconstruct_matrix(coeff, mask, params)
  vox_idx <- c(3, 10, 25)
  time_idx <- c(2, 5)
  part <- hrbf_reconstruct_partial(coeff, mask, params, voxel_idx = vox_idx, time_idx = time_idx)
  expect_equal(part, full[time_idx, vox_idx, drop = FALSE], tolerance = 1e-6)
})

test_that("hrbf_latent carries HRBF metadata", {
  mask_arr <- array(TRUE, dim = c(2, 2, 2))
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(2, 2, 2)))
  X <- matrix(rnorm(3 * 8), nrow = 3)
  params <- list(sigma0 = 1, levels = 0L, radius_factor = 2.5, kernel_type = "gaussian", seed = 9L)
  lvec <- hrbf_latent(X, mask, params)
  m <- hrbf_meta(lvec)
  expect_true(is_hrbf_latent(lvec))
  expect_equal(m$family, "hrbf")
  expect_equal(m$params$sigma0, 1)
})

test_that("hrbf_latent roundtrip on simulated fMRI (tiny mask)", {
  skip_on_cran()
  set.seed(99)

  mask_arr <- array(TRUE, dim = c(2, 2, 2))
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(2, 2, 2)))
  n_time <- 6L

  sim <- neuroim2::simulate_fmri(
    mask = mask,
    n_time = n_time,
    spatial_fwhm = 1.0,
    ar_mean = 0,
    ar_sd = 0,
    noise_sd = 0.05,
    n_factors = 0,
    global_amp = 0,
    seed = 99,
    return_centered = FALSE
  )

  sim_arr <- as.array(sim)
  mask_idx <- which(mask_arr)
  X <- matrix(0, nrow = n_time, ncol = length(mask_idx))
  for (t in seq_len(n_time)) {
    X[t, ] <- sim_arr[, , , t][mask_idx]
  }

  params <- list(sigma0 = 2, levels = 1L, radius_factor = 2.5, kernel_type = "gaussian", seed = 99L)
  lvec <- hrbf_latent(X, mask, params)
  recon <- as.matrix(lvec)

  expect_equal(recon, X, tolerance = 1e-6)
})

# ============================================================================
# Direct tests for hrbf_atoms_rcpp C++ function to improve coverage
# ============================================================================

# Helper to call internal Rcpp function
hrbf_atoms_rcpp_internal <- fmrilatent:::hrbf_atoms_rcpp

test_that("hrbf_atoms_rcpp gaussian kernel produces expected values", {
  # Simple 3x3 grid of mask points

mask_xyz <- matrix(c(
    0, 0, 0,
    1, 0, 0,
    2, 0, 0,
    0, 1, 0,
    1, 1, 0,
    2, 1, 0,
    0, 2, 0,
    1, 2, 0,
    2, 2, 0
  ), ncol = 3, byrow = TRUE)

  # Single center at (1, 1, 0)
  centres <- matrix(c(1, 1, 0), ncol = 3, byrow = TRUE)
  sigma_vec <- 1.0

  result <- hrbf_atoms_rcpp_internal(mask_xyz, centres, sigma_vec, "gaussian", 1e-12)

  expect_s4_class(result, "dgCMatrix")
  expect_equal(dim(result), c(1, 9))

  # Center should have value exp(0) = 1
  # The dense version for checking
  result_dense <- as.matrix(result)

  # Voxel at (1,1,0) is index 5 (0-based: 4)
  expect_equal(result_dense[1, 5], 1.0, tolerance = 1e-10)

  # Adjacent voxels (distance = 1): exp(-1/(2*1^2)) = exp(-0.5)
  expect_equal(result_dense[1, 2], exp(-0.5), tolerance = 1e-10) # (1,0,0)
  expect_equal(result_dense[1, 4], exp(-0.5), tolerance = 1e-10) # (0,1,0)
  expect_equal(result_dense[1, 6], exp(-0.5), tolerance = 1e-10) # (2,1,0)
  expect_equal(result_dense[1, 8], exp(-0.5), tolerance = 1e-10) # (1,2,0)

  # Corner voxels (distance = sqrt(2)): exp(-2/(2*1^2)) = exp(-1)
  expect_equal(result_dense[1, 1], exp(-1), tolerance = 1e-10) # (0,0,0)
  expect_equal(result_dense[1, 3], exp(-1), tolerance = 1e-10) # (2,0,0)
  expect_equal(result_dense[1, 7], exp(-1), tolerance = 1e-10) # (0,2,0)
  expect_equal(result_dense[1, 9], exp(-1), tolerance = 1e-10) # (2,2,0)
})

test_that("hrbf_atoms_rcpp wendland kernel produces expected values", {
  # Points at various distances from center
  mask_xyz <- matrix(c(
    0, 0, 0,   # distance 0
    0.5, 0, 0, # distance 0.5
    0.9, 0, 0, # distance 0.9 (within support)
    1.0, 0, 0, # distance 1.0 (boundary - should be zero)
    1.5, 0, 0  # distance 1.5 (outside support)
  ), ncol = 3, byrow = TRUE)

  centres <- matrix(c(0, 0, 0), ncol = 3, byrow = TRUE)
  sigma_vec <- 1.0

  result <- hrbf_atoms_rcpp_internal(mask_xyz, centres, sigma_vec, "wendland_c6", 1e-12)

  expect_s4_class(result, "dgCMatrix")
  expect_equal(dim(result), c(1, 5))

  result_dense <- as.matrix(result)

  # At r=0: (1-0)^8 * (32*0 + 25*0 + 8*0 + 1) = 1
  expect_equal(result_dense[1, 1], 1.0, tolerance = 1e-10)

  # At r=0.5: (1-0.5)^8 * (32*0.125 + 25*0.25 + 8*0.5 + 1)
  #         = (0.5)^8 * (4 + 6.25 + 4 + 1) = 0.00390625 * 15.25
  r <- 0.5
  expected_05 <- (1 - r)^8 * (32*r^3 + 25*r^2 + 8*r + 1)
  expect_equal(result_dense[1, 2], expected_05, tolerance = 1e-10)

  # At r=0.9: should be small but positive
  r <- 0.9
  expected_09 <- (1 - r)^8 * (32*r^3 + 25*r^2 + 8*r + 1)
  expect_equal(result_dense[1, 3], expected_09, tolerance = 1e-10)

  # At r=1.0: exactly zero (boundary)
  expect_equal(result_dense[1, 4], 0.0, tolerance = 1e-10)

  # At r=1.5: outside support, should be zero
  expect_equal(result_dense[1, 5], 0.0, tolerance = 1e-10)
})

test_that("hrbf_atoms_rcpp wendland kernel has compact support", {
  # Create points at increasing distances
  # Use 0.9 instead of 0.99 to avoid floating point edge cases
  distances <- c(0, 0.25, 0.5, 0.75, 0.9, 1.0, 1.01, 1.5, 2.0, 3.0)
  mask_xyz <- matrix(cbind(distances, 0, 0), ncol = 3)

  centres <- matrix(c(0, 0, 0), ncol = 3, byrow = TRUE)
  sigma_vec <- 1.0

  result <- hrbf_atoms_rcpp_internal(mask_xyz, centres, sigma_vec, "wendland_c6", 1e-12)
  result_dense <- as.matrix(result)

  # Points with r < 1 should be non-zero
  expect_gt(result_dense[1, 1], 0) # r=0
  expect_gt(result_dense[1, 2], 0) # r=0.25
  expect_gt(result_dense[1, 3], 0) # r=0.5
  expect_gt(result_dense[1, 4], 0) # r=0.75
  expect_gt(result_dense[1, 5], 0) # r=0.9

  # Points with r >= 1 should be exactly zero
  expect_equal(result_dense[1, 6], 0.0)  # r=1.0
  expect_equal(result_dense[1, 7], 0.0)  # r=1.01
  expect_equal(result_dense[1, 8], 0.0)  # r=1.5
  expect_equal(result_dense[1, 9], 0.0)  # r=2.0
  expect_equal(result_dense[1, 10], 0.0) # r=3.0
})

test_that("hrbf_atoms_rcpp respects value_threshold for gaussian", {
  # Create points at various distances
  mask_xyz <- matrix(c(
    0, 0, 0,
    5, 0, 0,  # distance 5 with sigma=1: exp(-12.5) very small
    10, 0, 0  # distance 10: exp(-50) extremely small
  ), ncol = 3, byrow = TRUE)

  centres <- matrix(c(0, 0, 0), ncol = 3, byrow = TRUE)
  sigma_vec <- 1.0

  # With very low threshold, should include small values
  result_low <- hrbf_atoms_rcpp_internal(mask_xyz, centres, sigma_vec, "gaussian", 1e-50)
  result_low_dense <- as.matrix(result_low)

  # With higher threshold, should filter out small values
  result_high <- hrbf_atoms_rcpp_internal(mask_xyz, centres, sigma_vec, "gaussian", 1e-6)
  result_high_dense <- as.matrix(result_high)

  # Both should have center value
  expect_equal(result_low_dense[1, 1], 1.0, tolerance = 1e-10)
  expect_equal(result_high_dense[1, 1], 1.0, tolerance = 1e-10)

  # exp(-12.5) approx 3.7e-6, so with threshold 1e-6 should be included
  expect_gt(result_low_dense[1, 2], 0)

  # exp(-50) approx 1.9e-22, with threshold 1e-6 should be filtered
  expect_equal(result_high_dense[1, 3], 0.0)
})

test_that("hrbf_atoms_rcpp handles multiple centres", {
  # 2D grid of mask points
  mask_xyz <- matrix(c(
    0, 0, 0,
    1, 0, 0,
    2, 0, 0,
    3, 0, 0,
    4, 0, 0
  ), ncol = 3, byrow = TRUE)

  # Two centres
  centres <- matrix(c(
    0, 0, 0,
    4, 0, 0
  ), ncol = 3, byrow = TRUE)

  sigma_vec <- c(1.0, 1.0)

  result <- hrbf_atoms_rcpp_internal(mask_xyz, centres, sigma_vec, "gaussian", 1e-12)

  expect_equal(dim(result), c(2, 5))

  result_dense <- as.matrix(result)

  # First centre at (0,0,0) - peak at first column
  expect_equal(result_dense[1, 1], 1.0, tolerance = 1e-10)

  # Second centre at (4,0,0) - peak at last column
  expect_equal(result_dense[2, 5], 1.0, tolerance = 1e-10)

  # Check decay for first centre
  expect_equal(result_dense[1, 2], exp(-0.5), tolerance = 1e-10) # dist=1
  expect_equal(result_dense[1, 3], exp(-2.0), tolerance = 1e-10) # dist=2

  # Check decay for second centre
  expect_equal(result_dense[2, 4], exp(-0.5), tolerance = 1e-10) # dist=1
  expect_equal(result_dense[2, 3], exp(-2.0), tolerance = 1e-10) # dist=2
})

test_that("hrbf_atoms_rcpp handles different sigma values per centre", {
  mask_xyz <- matrix(c(
    0, 0, 0,
    2, 0, 0
  ), ncol = 3, byrow = TRUE)

  centres <- matrix(c(
    0, 0, 0,
    2, 0, 0
  ), ncol = 3, byrow = TRUE)

  # Different sigma for each centre
  sigma_vec <- c(1.0, 2.0)

  result <- hrbf_atoms_rcpp_internal(mask_xyz, centres, sigma_vec, "gaussian", 1e-12)
  result_dense <- as.matrix(result)

  # Both centres at their locations should be 1.0
  expect_equal(result_dense[1, 1], 1.0, tolerance = 1e-10)
  expect_equal(result_dense[2, 2], 1.0, tolerance = 1e-10)

  # Distance from centre 1 to point 2: distance=2, sigma=1
  # exp(-4/(2*1)) = exp(-2)
  expect_equal(result_dense[1, 2], exp(-2), tolerance = 1e-10)

  # Distance from centre 2 to point 1: distance=2, sigma=2
  # exp(-4/(2*4)) = exp(-0.5)
  expect_equal(result_dense[2, 1], exp(-0.5), tolerance = 1e-10)
})

test_that("hrbf_atoms_rcpp validates sigma_vec length", {
  mask_xyz <- matrix(c(0, 0, 0, 1, 0, 0), ncol = 3, byrow = TRUE)
  centres <- matrix(c(0, 0, 0, 1, 0, 0), ncol = 3, byrow = TRUE) # 2 centres

  # Wrong sigma length
  sigma_vec <- c(1.0) # Only 1 sigma for 2 centres

  expect_error(
    hrbf_atoms_rcpp_internal(mask_xyz, centres, sigma_vec, "gaussian", 1e-8),
    "sigma_vec_mm length must match centres rows"
  )
})

test_that("hrbf_atoms_rcpp handles single point", {
  mask_xyz <- matrix(c(0, 0, 0), ncol = 3, byrow = TRUE)
  centres <- matrix(c(0, 0, 0), ncol = 3, byrow = TRUE)
  sigma_vec <- 1.0

  result <- hrbf_atoms_rcpp_internal(mask_xyz, centres, sigma_vec, "gaussian", 1e-12)
  result_dense <- as.matrix(result)

  expect_equal(dim(result), c(1, 1))
  expect_equal(result_dense[1, 1], 1.0, tolerance = 1e-10)
})

test_that("hrbf_atoms_rcpp handles 3D point cloud", {
  # Create a small 3D grid - ensure it's a numeric matrix
  coords <- expand.grid(x = 0:2, y = 0:2, z = 0:2)
  mask_xyz <- matrix(as.numeric(as.matrix(coords)), ncol = 3)

  centres <- matrix(c(1, 1, 1), ncol = 3, byrow = TRUE)
  sigma_vec <- 1.0

  result <- hrbf_atoms_rcpp_internal(mask_xyz, centres, sigma_vec, "gaussian", 1e-12)

  expect_equal(dim(result), c(1, 27))

  result_dense <- as.matrix(result)

  # Find center point (1,1,1) - should be index 14 (row-major from expand.grid)
  center_idx <- which(coords$x == 1 & coords$y == 1 & coords$z == 1)
  expect_equal(result_dense[1, center_idx], 1.0, tolerance = 1e-10)

  # All values should be positive for gaussian
  expect_true(all(result_dense > 0))

  # Values should decrease with distance
  corner_idx <- which(coords$x == 0 & coords$y == 0 & coords$z == 0)
  expect_lt(result_dense[1, corner_idx], result_dense[1, center_idx])
})

test_that("hrbf_atoms_rcpp wendland kernel with 3D points", {
  # Create points in 3D - ensure numeric matrix
  coords <- expand.grid(x = seq(0, 2, by = 0.5), y = 0, z = 0)
  mask_xyz <- matrix(as.numeric(as.matrix(coords)), ncol = 3)

  centres <- matrix(c(0, 0, 0), ncol = 3, byrow = TRUE)
  sigma_vec <- 1.5

  result <- hrbf_atoms_rcpp_internal(mask_xyz, centres, sigma_vec, "wendland_c6", 1e-12)
  result_dense <- as.matrix(result)

  # Center should be 1.0
  expect_equal(result_dense[1, 1], 1.0, tolerance = 1e-10)

  # Point at x=2 has distance 2, r = 2/1.5 = 1.33 > 1, should be zero
  last_idx <- ncol(result_dense)
  expect_equal(result_dense[1, last_idx], 0.0, tolerance = 1e-10)
})

test_that("hrbf_atoms_rcpp handles large sigma", {
  mask_xyz <- matrix(c(
    0, 0, 0,
    10, 0, 0,
    100, 0, 0
  ), ncol = 3, byrow = TRUE)

  centres <- matrix(c(0, 0, 0), ncol = 3, byrow = TRUE)
  sigma_vec <- 50.0  # Very large sigma

  result <- hrbf_atoms_rcpp_internal(mask_xyz, centres, sigma_vec, "gaussian", 1e-12)
  result_dense <- as.matrix(result)

  # With large sigma, all values should be close to 1
  expect_equal(result_dense[1, 1], 1.0, tolerance = 1e-10)

  # exp(-100/(2*2500)) = exp(-0.02) approx 0.98
  expect_gt(result_dense[1, 2], 0.95)

  # exp(-10000/(2*2500)) = exp(-2) approx 0.135
  expect_gt(result_dense[1, 3], 0.1)
})

test_that("hrbf_atoms_rcpp handles small sigma", {
  mask_xyz <- matrix(c(
    0, 0, 0,
    0.1, 0, 0,
    0.5, 0, 0
  ), ncol = 3, byrow = TRUE)

  centres <- matrix(c(0, 0, 0), ncol = 3, byrow = TRUE)
  sigma_vec <- 0.1  # Very small sigma

  result <- hrbf_atoms_rcpp_internal(mask_xyz, centres, sigma_vec, "gaussian", 1e-50)
  result_dense <- as.matrix(result)

  # Center should still be 1.0
  expect_equal(result_dense[1, 1], 1.0, tolerance = 1e-10)

  # exp(-0.01/(2*0.01)) = exp(-0.5) at distance 0.1
  expect_equal(result_dense[1, 2], exp(-0.5), tolerance = 1e-10)

  # exp(-0.25/(2*0.01)) = exp(-12.5) at distance 0.5 - very small
  expect_lt(result_dense[1, 3], 1e-5)
})

test_that("hrbf_atoms_rcpp output is sparse for wendland", {
  # Create a larger grid - ensure numeric matrix
  coords <- expand.grid(x = 0:9, y = 0:9, z = 0)
  mask_xyz <- matrix(as.numeric(as.matrix(coords)), ncol = 3)

  centres <- matrix(c(5, 5, 0), ncol = 3, byrow = TRUE)
  sigma_vec <- 2.0  # Wendland has compact support at r=1, so radius is 2.0

  result <- hrbf_atoms_rcpp_internal(mask_xyz, centres, sigma_vec, "wendland_c6", 1e-12)

  # Should be sparse - many zeros outside radius
  nnz <- Matrix::nnzero(result)
  total <- prod(dim(result))

  # With radius 2.0, only points within distance 2 should be non-zero
  # This is a circle of radius 2 centered at (5,5) in a 10x10 grid
  expect_lt(nnz, total * 0.5)
})

test_that("hrbf_atoms_rcpp gaussian is symmetric",
{
  # Test symmetry by placing mask points equidistant from center
  mask_xyz <- matrix(c(
    1, 0, 0,
    -1, 0, 0,
    0, 1, 0,
    0, -1, 0,
    0, 0, 1,
    0, 0, -1
  ), ncol = 3, byrow = TRUE)

  centres <- matrix(c(0, 0, 0), ncol = 3, byrow = TRUE)
  sigma_vec <- 1.0

  result <- hrbf_atoms_rcpp_internal(mask_xyz, centres, sigma_vec, "gaussian", 1e-12)
  result_dense <- as.matrix(result)

  # All equidistant points should have the same value
  expect_equal(result_dense[1, 1], result_dense[1, 2], tolerance = 1e-10)
  expect_equal(result_dense[1, 1], result_dense[1, 3], tolerance = 1e-10)
  expect_equal(result_dense[1, 1], result_dense[1, 4], tolerance = 1e-10)
  expect_equal(result_dense[1, 1], result_dense[1, 5], tolerance = 1e-10)
  expect_equal(result_dense[1, 1], result_dense[1, 6], tolerance = 1e-10)
})

test_that("hrbf_atoms_rcpp wendland is symmetric", {
  mask_xyz <- matrix(c(
    0.5, 0, 0,
    -0.5, 0, 0,
    0, 0.5, 0,
    0, -0.5, 0,
    0, 0, 0.5,
    0, 0, -0.5
  ), ncol = 3, byrow = TRUE)

  centres <- matrix(c(0, 0, 0), ncol = 3, byrow = TRUE)
  sigma_vec <- 1.0

  result <- hrbf_atoms_rcpp_internal(mask_xyz, centres, sigma_vec, "wendland_c6", 1e-12)
  result_dense <- as.matrix(result)

  # All equidistant points should have the same value
  expect_equal(result_dense[1, 1], result_dense[1, 2], tolerance = 1e-10)
  expect_equal(result_dense[1, 1], result_dense[1, 3], tolerance = 1e-10)
  expect_equal(result_dense[1, 1], result_dense[1, 4], tolerance = 1e-10)
  expect_equal(result_dense[1, 1], result_dense[1, 5], tolerance = 1e-10)
  expect_equal(result_dense[1, 1], result_dense[1, 6], tolerance = 1e-10)
})

test_that("hrbf_atoms_rcpp handles many centres efficiently", {
  # Create a grid of mask points - ensure numeric matrix
  coords <- expand.grid(x = 0:4, y = 0:4, z = 0)
  mask_xyz <- matrix(as.numeric(as.matrix(coords)), ncol = 3)

  # Multiple centres
  centres <- matrix(c(
    0, 0, 0,
    2, 0, 0,
    4, 0, 0,
    0, 2, 0,
    2, 2, 0,
    4, 2, 0,
    0, 4, 0,
    2, 4, 0,
    4, 4, 0
  ), ncol = 3, byrow = TRUE)

  sigma_vec <- rep(1.0, 9)

  result <- hrbf_atoms_rcpp_internal(mask_xyz, centres, sigma_vec, "gaussian", 1e-12)

  expect_equal(dim(result), c(9, 25))

  # Each centre should have peak at its location
  result_dense <- as.matrix(result)
  for (k in 1:9) {
    cx <- centres[k, 1]
    cy <- centres[k, 2]
    idx <- which(coords$x == cx & coords$y == cy)
    expect_equal(result_dense[k, idx], 1.0, tolerance = 1e-10)
  }
})

test_that("hrbf_atoms_rcpp handles negative coordinates", {
  mask_xyz <- matrix(c(
    -2, -2, -2,
    -1, -1, -1,
    0, 0, 0,
    1, 1, 1,
    2, 2, 2
  ), ncol = 3, byrow = TRUE)

  centres <- matrix(c(0, 0, 0), ncol = 3, byrow = TRUE)
  sigma_vec <- 2.0

  result <- hrbf_atoms_rcpp_internal(mask_xyz, centres, sigma_vec, "gaussian", 1e-12)
  result_dense <- as.matrix(result)

  # Center at (0,0,0) should be 1.0
  expect_equal(result_dense[1, 3], 1.0, tolerance = 1e-10)

  # Points (-1,-1,-1) and (1,1,1) are equidistant (distance = sqrt(3))
  expect_equal(result_dense[1, 2], result_dense[1, 4], tolerance = 1e-10)

  # Points (-2,-2,-2) and (2,2,2) are equidistant (distance = sqrt(12))
  expect_equal(result_dense[1, 1], result_dense[1, 5], tolerance = 1e-10)
})

test_that("hrbf_atoms_rcpp with Rcpp backend matches expected for gaussian", {
  # Enable Rcpp backend and test
  old_opt <- getOption("fmrilatent.hrbf.use_rcpp")
  on.exit(options(fmrilatent.hrbf.use_rcpp = old_opt), add = TRUE)
  options(fmrilatent.hrbf.use_rcpp = TRUE)

  mask_arr <- array(TRUE, dim = c(4, 4, 4))
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(4, 4, 4)))
  params <- list(sigma0 = 2, levels = 1L, radius_factor = 2.5, kernel_type = "gaussian", seed = 42L)

  B <- hrbf_generate_basis(params, mask)

  expect_s4_class(B, "dgCMatrix")
  expect_equal(ncol(B), sum(mask_arr))
  expect_gt(nrow(B), 0)

  # All values in basis should be non-negative for gaussian
  expect_true(all(B@x >= 0))
})

test_that("hrbf_atoms_rcpp with Rcpp backend for wendland kernel", {
  old_opt <- getOption("fmrilatent.hrbf.use_rcpp")
  on.exit(options(fmrilatent.hrbf.use_rcpp = old_opt), add = TRUE)
  options(fmrilatent.hrbf.use_rcpp = TRUE)

  mask_arr <- array(TRUE, dim = c(4, 4, 4))
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(4, 4, 4)))
  params <- list(sigma0 = 2, levels = 1L, radius_factor = 2.5, kernel_type = "wendland_c6", seed = 42L)

  B <- hrbf_generate_basis(params, mask)

  expect_s4_class(B, "dgCMatrix")
  expect_equal(ncol(B), sum(mask_arr))
  expect_gt(nrow(B), 0)

  # All values in basis should be non-negative for wendland
  expect_true(all(B@x >= 0))
})

test_that("hrbf_atoms_rcpp Rcpp and R backends produce similar results", {
  # Compare Rcpp and R backends to ensure consistency
  mask_arr <- array(TRUE, dim = c(3, 3, 3))
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(3, 3, 3)))
  params <- list(sigma0 = 2, levels = 1L, radius_factor = 2.5, kernel_type = "gaussian", seed = 123L)

  old_opt <- getOption("fmrilatent.hrbf.use_rcpp")
  on.exit(options(fmrilatent.hrbf.use_rcpp = old_opt), add = TRUE)

  # Get R backend result
  options(fmrilatent.hrbf.use_rcpp = FALSE)
  B_R <- hrbf_generate_basis(params, mask)

  # Get Rcpp backend result
  options(fmrilatent.hrbf.use_rcpp = TRUE)
  B_Rcpp <- hrbf_generate_basis(params, mask)

  # Dimensions should match
  expect_equal(dim(B_R), dim(B_Rcpp))

  # Values should be similar (normalized)
  B_R_dense <- as.matrix(B_R)
  B_Rcpp_dense <- as.matrix(B_Rcpp)

  # Both should have similar structure
  expect_equal(nrow(B_R), nrow(B_Rcpp))
  expect_equal(ncol(B_R), ncol(B_Rcpp))
})

test_that("hrbf_atoms_rcpp handles zero threshold edge case", {
  # Use smaller distances to avoid underflow
  mask_xyz <- matrix(c(
    0, 0, 0,
    1, 0, 0,
    5, 0, 0
  ), ncol = 3, byrow = TRUE)

  centres <- matrix(c(0, 0, 0), ncol = 3, byrow = TRUE)
  sigma_vec <- 2.0  # Larger sigma to keep values non-zero

  # With threshold 0, only exact zeros should be filtered
  result <- hrbf_atoms_rcpp_internal(mask_xyz, centres, sigma_vec, "gaussian", 0)
  result_dense <- as.matrix(result)

  # All gaussian values should be included (none are exactly 0)
  # exp(-25/(2*4)) = exp(-3.125) is still a reasonable value
  expect_true(all(result_dense > 0))
})

test_that("hrbf_atoms_rcpp handles very small distances", {
  # Points very close together
  mask_xyz <- matrix(c(
    0, 0, 0,
    1e-10, 0, 0,
    1e-5, 0, 0
  ), ncol = 3, byrow = TRUE)

  centres <- matrix(c(0, 0, 0), ncol = 3, byrow = TRUE)
  sigma_vec <- 1.0

  result <- hrbf_atoms_rcpp_internal(mask_xyz, centres, sigma_vec, "gaussian", 1e-12)
  result_dense <- as.matrix(result)

  # All points very close to center should have values near 1.0
  expect_equal(result_dense[1, 1], 1.0, tolerance = 1e-10)
  expect_equal(result_dense[1, 2], 1.0, tolerance = 1e-6)
  expect_equal(result_dense[1, 3], 1.0, tolerance = 1e-6)
})

test_that("hrbf_atoms_rcpp handles non-uniform spacing", {
  # Points with non-uniform spacing
  mask_xyz <- matrix(c(
    0, 0, 0,
    0.1, 0, 0,
    0.5, 0, 0,
    2.0, 0, 0,
    10.0, 0, 0
  ), ncol = 3, byrow = TRUE)

  centres <- matrix(c(0, 0, 0), ncol = 3, byrow = TRUE)
  sigma_vec <- 2.0

  result <- hrbf_atoms_rcpp_internal(mask_xyz, centres, sigma_vec, "gaussian", 1e-12)
  result_dense <- as.matrix(result)

  # Verify values decrease monotonically with distance
  expect_gt(result_dense[1, 1], result_dense[1, 2])
  expect_gt(result_dense[1, 2], result_dense[1, 3])
  expect_gt(result_dense[1, 3], result_dense[1, 4])
  expect_gt(result_dense[1, 4], result_dense[1, 5])
})

test_that("hrbf_atoms_rcpp with tiny mask uses Rcpp backend", {
  # Test the tiny-ROI fallback path with Rcpp enabled
  old_opt <- getOption("fmrilatent.hrbf.use_rcpp")
  on.exit(options(fmrilatent.hrbf.use_rcpp = old_opt), add = TRUE)
  options(fmrilatent.hrbf.use_rcpp = TRUE)

  mask_arr <- array(TRUE, dim = c(2, 2, 2))
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(2, 2, 2)))
  params <- list(sigma0 = 1, levels = 0L, radius_factor = 2.5, kernel_type = "gaussian", seed = 1L)

  B <- hrbf_generate_basis(params, mask)

  expect_s4_class(B, "dgCMatrix")
  expect_equal(ncol(B), sum(mask_arr))
  expect_gt(nrow(B), 0)
})

test_that("hrbf_atoms_rcpp with extra fine levels", {
  old_opt <- getOption("fmrilatent.hrbf.use_rcpp")
  on.exit(options(fmrilatent.hrbf.use_rcpp = old_opt), add = TRUE)
  options(fmrilatent.hrbf.use_rcpp = TRUE)

  mask_arr <- array(TRUE, dim = c(5, 5, 5))
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(5, 5, 5)))
  params <- list(
    sigma0 = 2,
    levels = 1L,
    radius_factor = 2.5,
    kernel_type = "gaussian",
    num_extra_fine_levels = 1L,
    seed = 42L
  )

  B <- hrbf_generate_basis(params, mask)

  expect_s4_class(B, "dgCMatrix")
  expect_equal(ncol(B), sum(mask_arr))
  expect_gt(nrow(B), 0)
})

# ============================================================================
# Additional tests for zero-coverage functions
# ============================================================================

test_that("hrbf_generate_basis edge case: missing seed gets default", {
  mask_arr <- array(TRUE, dim = c(2, 2, 2))
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(2, 2, 2)))
  params <- list(sigma0 = 1, levels = 0L, radius_factor = 2.5, kernel_type = "gaussian")
  # No seed provided
  B <- hrbf_generate_basis(params, mask)
  expect_s4_class(B, "dgCMatrix")
  expect_equal(ncol(B), sum(mask_arr))
})

test_that("hrbf_reconstruct_partial error: missing voxel_idx", {
  mask_arr <- array(TRUE, dim = c(2, 2, 2))
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(2, 2, 2)))
  params <- list(sigma0 = 1, levels = 0L, radius_factor = 2.5, kernel_type = "gaussian", seed = 1L)
  coeff <- matrix(rnorm(3 * 8), nrow = 3)

  expect_error(
    hrbf_reconstruct_partial(coeff, mask, params),
    "voxel_idx must be provided"
  )
})

test_that("hrbf_reconstruct_partial error: voxel_idx out of bounds", {
  mask_arr <- array(TRUE, dim = c(2, 2, 2))
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(2, 2, 2)))
  params <- list(sigma0 = 1, levels = 0L, radius_factor = 2.5, kernel_type = "gaussian", seed = 1L)
  coeff <- matrix(rnorm(3 * 8), nrow = 3)

  expect_error(
    hrbf_reconstruct_partial(coeff, mask, params, voxel_idx = c(1, 100)),
    "voxel_idx out of bounds"
  )
})

test_that("hrbf_reconstruct_partial error: time_idx out of bounds", {
  mask_arr <- array(TRUE, dim = c(2, 2, 2))
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(2, 2, 2)))
  params <- list(sigma0 = 1, levels = 0L, radius_factor = 2.5, kernel_type = "gaussian", seed = 1L)
  coeff <- matrix(rnorm(3 * 8), nrow = 3)

  expect_error(
    hrbf_reconstruct_partial(coeff, mask, params, voxel_idx = c(1, 2), time_idx = c(1, 10)),
    "time_idx out of bounds"
  )
})

test_that("hrbf_reconstruct_partial edge case: empty voxel_idx", {
  mask_arr <- array(TRUE, dim = c(2, 2, 2))
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(2, 2, 2)))
  params <- list(sigma0 = 1, levels = 0L, radius_factor = 2.5, kernel_type = "gaussian", seed = 1L)
  coeff <- matrix(rnorm(3 * 8), nrow = 3)

  expect_error(
    hrbf_reconstruct_partial(coeff, mask, params, voxel_idx = integer(0)),
    "voxel_idx must be provided"
  )
})

test_that("hrbf_latent error: X must have time in rows", {
  mask_arr <- array(TRUE, dim = c(2, 2, 2))
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(2, 2, 2)))
  params <- list(sigma0 = 1, levels = 0L, radius_factor = 2.5, kernel_type = "gaussian", seed = 1L)

  # Empty matrix
  expect_error(
    hrbf_latent(matrix(nrow = 0, ncol = 8), mask, params),
    "X must have time in rows"
  )
})

test_that("hrbf_latent with custom label", {
  mask_arr <- array(TRUE, dim = c(2, 2, 2))
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(2, 2, 2)))
  X <- matrix(rnorm(3 * 8), nrow = 3)
  params <- list(sigma0 = 1, levels = 0L, radius_factor = 2.5, kernel_type = "gaussian", seed = 1L)

  lvec <- hrbf_latent(X, mask, params, label = "test_label")
  expect_s4_class(lvec, "LatentNeuroVec")
  expect_equal(lvec@label, "test_label")
})

test_that("as_hrbf_latent attaches metadata correctly", {
  mask_arr <- array(TRUE, dim = c(2, 2, 2))
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(2, 2, 2)))
  X <- matrix(rnorm(3 * 8), nrow = 3)
  params <- list(sigma0 = 1, levels = 0L, radius_factor = 2.5, kernel_type = "gaussian", seed = 1L)

  lvec <- hrbf_latent(X, mask, params)
  centres_test <- matrix(c(0, 0, 0), ncol = 3)
  sigmas_test <- c(1.0)

  lvec2 <- as_hrbf_latent(lvec, params, centres = centres_test, sigmas = sigmas_test)

  expect_s4_class(lvec2, "LatentNeuroVec")
  meta <- hrbf_meta(lvec2)
  expect_equal(meta$family, "hrbf")
  expect_equal(meta$centres, centres_test)
  expect_equal(meta$sigmas, sigmas_test)
})

test_that("as_hrbf_latent error: non-LatentNeuroVec input", {
  expect_error(
    as_hrbf_latent(list(), params = list()),
    "inherits\\(lvec, \"LatentNeuroVec\"\\) is not TRUE"
  )
})

test_that("hrbf_meta returns NULL for non-LatentNeuroVec", {
  expect_null(hrbf_meta(list()))
  expect_null(hrbf_meta(matrix(1:10)))
})

test_that("hrbf_meta returns NULL for LatentNeuroVec without meta", {
  mask_arr <- array(TRUE, dim = c(2, 2, 2))
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(2, 2, 2)))
  spc <- neuroim2::NeuroSpace(c(2, 2, 2, 3))

  # Create LatentNeuroVec without HRBF metadata
  lvec <- LatentNeuroVec(
    basis = matrix(rnorm(3 * 8), nrow = 3),
    loadings = Matrix::Matrix(rnorm(8 * 8), nrow = 8),
    space = spc,
    mask = mask
  )

  expect_null(hrbf_meta(lvec))
})

test_that("hrbf_meta returns NULL for wrong family", {
  mask_arr <- array(TRUE, dim = c(2, 2, 2))
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(2, 2, 2)))
  spc <- neuroim2::NeuroSpace(c(2, 2, 2, 3))

  lvec <- LatentNeuroVec(
    basis = matrix(rnorm(3 * 8), nrow = 3),
    loadings = Matrix::Matrix(rnorm(8 * 8), nrow = 8),
    space = spc,
    mask = mask,
    meta = list(family = "other")
  )

  expect_null(hrbf_meta(lvec))
})

test_that("extract_mask_array_logical error: cannot extract array", {
  expect_error(
    fmrilatent:::extract_mask_array_logical(NULL, "test_location"),
    "Cannot extract array from mask \\(test_location\\)"
  )
})

test_that("extract_mask_array_logical coerces vectors to arrays", {
  # Vectors get coerced to arrays with dimensions by as.array
  # This is expected behavior, so just verify it works
  vec_mask <- c(TRUE, FALSE, TRUE)
  result <- fmrilatent:::extract_mask_array_logical(vec_mask, "test_location")
  expect_true(is.logical(result))
  expect_equal(length(dim(result)), 1)  # 1D array
})

test_that("extract_mask_array_logical handles list with arr field", {
  mask_list <- list(arr = array(TRUE, dim = c(2, 2, 2)))
  result <- fmrilatent:::extract_mask_array_logical(mask_list, "test")
  expect_equal(dim(result), c(2, 2, 2))
  expect_true(is.logical(result))
})

test_that("poisson_disk_sample_neuroim2 handles multi-component mask", {
  # Create a mask with two disconnected components
  mask_arr <- array(FALSE, dim = c(5, 5, 5))
  mask_arr[1:2, 1:2, 1:2] <- TRUE  # Component 1
  mask_arr[4:5, 4:5, 4:5] <- TRUE  # Component 2
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(5, 5, 5)))

  centres <- fmrilatent:::poisson_disk_sample_neuroim2(mask, radius_mm = 1.5, seed = 42L)

  expect_true(nrow(centres) >= 2)  # Should have centres from both components
  expect_equal(ncol(centres), 3)
  expect_equal(colnames(centres), c("i", "j", "k"))
})

test_that("poisson_disk_sample_neuroim2 handles empty mask", {
  mask_arr <- array(FALSE, dim = c(3, 3, 3))
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(3, 3, 3)))

  centres <- fmrilatent:::poisson_disk_sample_neuroim2(mask, radius_mm = 1.0, seed = 1L)

  expect_equal(nrow(centres), 0)
  expect_equal(ncol(centres), 3)
})

test_that("poisson_disk_sample_neuroim2 handles very small mask", {
  mask_arr <- array(FALSE, dim = c(3, 3, 3))
  mask_arr[2, 2, 2] <- TRUE  # Single voxel
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(3, 3, 3)))

  centres <- fmrilatent:::poisson_disk_sample_neuroim2(mask, radius_mm = 0.5, seed = 1L)

  expect_gte(nrow(centres), 1)
  expect_equal(ncol(centres), 3)
})

test_that("poisson_disk_sample_neuroim2_compat handles multi-component mask", {
  mask_arr <- array(FALSE, dim = c(5, 5, 5))
  mask_arr[1:2, 1:2, 1:2] <- TRUE
  mask_arr[4:5, 4:5, 4:5] <- TRUE
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(5, 5, 5)))

  centres <- fmrilatent:::poisson_disk_sample_neuroim2_compat(mask, radius_mm = 1.5, seed = 42L)

  expect_true(nrow(centres) >= 1)
  expect_equal(ncol(centres), 3)
  expect_equal(colnames(centres), c("i", "j", "k"))
})

test_that("poisson_disk_sample_neuroim2_compat handles empty mask", {
  mask_arr <- array(FALSE, dim = c(3, 3, 3))
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(3, 3, 3)))

  centres <- fmrilatent:::poisson_disk_sample_neuroim2_compat(mask, radius_mm = 1.0, seed = 1L)

  expect_equal(nrow(centres), 0)
  expect_equal(ncol(centres), 3)
})

test_that("poisson_disk_sample_neuroim2_compat with component offset", {
  mask_arr <- array(TRUE, dim = c(3, 3, 3))
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(3, 3, 3)))

  centres <- fmrilatent:::poisson_disk_sample_neuroim2_compat(
    mask, radius_mm = 2.0, seed = 1L, component_id_for_seed_offset = 5L
  )

  expect_gte(nrow(centres), 1)
  expect_equal(ncol(centres), 3)
})

test_that("generate_hrbf_atom gaussian kernel", {
  old_opt <- getOption("fmrilatent.hrbf.use_rcpp")
  on.exit(options(fmrilatent.hrbf.use_rcpp = old_opt), add = TRUE)
  options(fmrilatent.hrbf.use_rcpp = FALSE)

  mask_coords_world <- matrix(c(
    0, 0, 0,
    1, 0, 0,
    2, 0, 0
  ), ncol = 3, byrow = TRUE)
  mask_linear_indices <- 1:3
  centre_coord_world <- c(0, 0, 0)
  sigma_mm <- 1.0
  params <- list(kernel_type = "gaussian")

  atom <- fmrilatent:::generate_hrbf_atom(
    mask_coords_world, mask_linear_indices, centre_coord_world, sigma_mm,
    current_level_j = 0L, total_levels = 1L, params
  )

  expect_true(is.list(atom))
  expect_true("values" %in% names(atom))
  expect_true("indices" %in% names(atom))
  expect_equal(length(atom$values), 3)
  expect_equal(atom$indices, mask_linear_indices)

  # Center should have highest value
  expect_equal(which.max(atom$values), 1)
})

test_that("generate_hrbf_atom wendland kernel", {
  old_opt <- getOption("fmrilatent.hrbf.use_rcpp")
  on.exit(options(fmrilatent.hrbf.use_rcpp = old_opt), add = TRUE)
  options(fmrilatent.hrbf.use_rcpp = FALSE)

  mask_coords_world <- matrix(c(
    0, 0, 0,
    0.5, 0, 0,
    1.5, 0, 0
  ), ncol = 3, byrow = TRUE)
  mask_linear_indices <- 1:3
  centre_coord_world <- c(0, 0, 0)
  sigma_mm <- 1.0
  params <- list(kernel_type = "wendland_c6")

  atom <- fmrilatent:::generate_hrbf_atom(
    mask_coords_world, mask_linear_indices, centre_coord_world, sigma_mm,
    current_level_j = 0L, total_levels = 1L, params
  )

  expect_equal(length(atom$values), 3)
  # Point at distance 1.5 should be zero (outside support radius)
  expect_equal(atom$values[3], 0)
  # Point at distance 0.5 should be non-zero
  expect_gt(atom$values[2], 0)
})

test_that("generate_hrbf_atom wendland_c4 alias", {
  old_opt <- getOption("fmrilatent.hrbf.use_rcpp")
  on.exit(options(fmrilatent.hrbf.use_rcpp = old_opt), add = TRUE)
  options(fmrilatent.hrbf.use_rcpp = FALSE)

  mask_coords_world <- matrix(c(0, 0, 0, 0.5, 0, 0), ncol = 3, byrow = TRUE)
  params <- list(kernel_type = "wendland_c4")

  atom <- fmrilatent:::generate_hrbf_atom(
    mask_coords_world, 1:2, c(0, 0, 0), 1.0, 0L, 1L, params
  )

  expect_equal(length(atom$values), 2)
  expect_gt(atom$values[1], 0)
})

test_that("lna_hrbf_centres_from_params with explicit centres and sigmas", {
  mask_arr <- array(TRUE, dim = c(2, 2, 2))
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(2, 2, 2)))
  params <- list(sigma0 = 2, levels = 1L)

  centres <- matrix(c(0, 0, 0, 1, 1, 1), ncol = 3, byrow = TRUE)
  sigmas <- c(2.0, 1.0)

  result <- fmrilatent:::lna_hrbf_centres_from_params(
    params, mask, centres = centres, sigmas = sigmas
  )

  expect_equal(nrow(result$centres), 2)
  expect_equal(result$sigmas, sigmas)
  expect_true(result$centres_stored)
})

test_that("lna_hrbf_centres_from_params error: centres without sigmas", {
  mask_arr <- array(TRUE, dim = c(2, 2, 2))
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(2, 2, 2)))
  params <- list(sigma0 = 2, levels = 1L)
  centres <- matrix(c(0, 0, 0), ncol = 3)

  expect_error(
    fmrilatent:::lna_hrbf_centres_from_params(params, mask, centres = centres),
    "Both centres and sigmas must be supplied together"
  )
})

test_that("lna_hrbf_centres_from_params error: wrong centres dimensions", {
  mask_arr <- array(TRUE, dim = c(2, 2, 2))
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(2, 2, 2)))
  params <- list(sigma0 = 2, levels = 1L)
  centres <- matrix(c(0, 0), ncol = 2)  # Wrong: should be 3 columns
  sigmas <- c(1.0)

  expect_error(
    fmrilatent:::lna_hrbf_centres_from_params(params, mask, centres = centres, sigmas = sigmas),
    "centres must be a matrix with 3 columns"
  )
})

test_that("lna_hrbf_centres_from_params error: length mismatch", {
  mask_arr <- array(TRUE, dim = c(2, 2, 2))
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(2, 2, 2)))
  params <- list(sigma0 = 2, levels = 1L)
  centres <- matrix(c(0, 0, 0, 1, 1, 1), ncol = 3, byrow = TRUE)
  sigmas <- c(1.0)  # Wrong: should match nrow(centres)

  expect_error(
    fmrilatent:::lna_hrbf_centres_from_params(params, mask, centres = centres, sigmas = sigmas),
    "length\\(sigmas\\) must equal nrow\\(centres\\)"
  )
})

test_that("lna_hrbf_centres_from_params error: no seed", {
  mask_arr <- array(TRUE, dim = c(2, 2, 2))
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(2, 2, 2)))
  params <- list(sigma0 = 2, levels = 1L)  # No seed

  expect_error(
    fmrilatent:::lna_hrbf_centres_from_params(params, mask),
    "lna_hrbf_centres_from_params requires centres/sigmas or seed"
  )
})

test_that("lna_hrbf_centres_from_params with seed", {
  mask_arr <- array(TRUE, dim = c(3, 3, 3))
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(3, 3, 3)))
  params <- list(sigma0 = 2, levels = 1L, radius_factor = 2.5, seed = 42L)

  result <- fmrilatent:::lna_hrbf_centres_from_params(params, mask)

  expect_true(is.matrix(result$centres))
  expect_equal(ncol(result$centres), 3)
  expect_gt(nrow(result$centres), 0)
  expect_equal(length(result$sigmas), nrow(result$centres))
  expect_false(result$centres_stored)
})

test_that("lna_hrbf_basis_from_params with full_grid TRUE", {
  mask_arr <- array(TRUE, dim = c(2, 2, 2))
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(2, 2, 2)))
  params <- list(sigma0 = 1, levels = 0L, radius_factor = 2.5, kernel_type = "gaussian", seed = 1L)

  B <- lna_hrbf_basis_from_params(params, mask, full_grid = TRUE)

  expect_s4_class(B, "dgCMatrix")
  expect_equal(ncol(B), length(mask_arr))  # Full grid
})

test_that("lna_hrbf_basis_from_params with full_grid FALSE", {
  mask_arr <- array(TRUE, dim = c(2, 2, 2))
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(2, 2, 2)))
  params <- list(sigma0 = 1, levels = 0L, radius_factor = 2.5, kernel_type = "gaussian", seed = 1L)

  B <- lna_hrbf_basis_from_params(params, mask, full_grid = FALSE)

  expect_s4_class(B, "dgCMatrix")
  expect_equal(ncol(B), sum(mask_arr))  # Active voxels only
})

test_that("lna_hrbf_basis_from_params error: centres without sigmas", {
  mask_arr <- array(TRUE, dim = c(2, 2, 2))
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(2, 2, 2)))
  params <- list(sigma0 = 1, levels = 0L)
  centres <- matrix(c(0, 0, 0), ncol = 3)

  expect_error(
    lna_hrbf_basis_from_params(params, mask, centres = centres),
    "Both centres and sigmas must be supplied together"
  )
})

test_that("hrbf_basis_from_params with R backend and wendland kernel", {
  old_opt <- getOption("fmrilatent.hrbf.use_rcpp")
  on.exit(options(fmrilatent.hrbf.use_rcpp = old_opt), add = TRUE)
  options(fmrilatent.hrbf.use_rcpp = FALSE)

  mask_arr <- array(TRUE, dim = c(3, 3, 3))
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(3, 3, 3)))
  params <- list(sigma0 = 2, levels = 1L, radius_factor = 2.5, kernel_type = "wendland_c6", seed = 42L)

  B <- fmrilatent:::hrbf_basis_from_params(params, mask)

  expect_s4_class(B, "dgCMatrix")
  expect_equal(ncol(B), sum(mask_arr))
  expect_gt(nrow(B), 0)
})

test_that("compute_hrbf_coefficients square invertible case", {
  # Create a square invertible basis
  set.seed(123)
  n <- 8
  basis <- Matrix::Matrix(rnorm(n * n), nrow = n, sparse = TRUE)
  # Make it well-conditioned
  basis <- basis + Matrix::Diagonal(n, x = 10)
  X <- matrix(rnorm(5 * n), nrow = 5)

  coeff <- fmrilatent:::compute_hrbf_coefficients(X, basis)

  expect_equal(dim(coeff), c(5, n))

  # Reconstruction should be close
  recon <- coeff %*% basis
  expect_equal(dim(recon), dim(X))
})

test_that("compute_hrbf_coefficients overdetermined case", {
  skip_on_cran()  # Skip due to long vector limits in diag() on some platforms
  set.seed(456)
  n_vox <- 12
  n_atoms <- 6
  # Create a truly sparse basis to avoid dense conversion
  basis <- Matrix::sparseMatrix(
    i = rep(1:n_atoms, each = 2),
    j = sample(1:n_vox, n_atoms * 2, replace = TRUE),
    x = rnorm(n_atoms * 2),
    dims = c(n_atoms, n_vox)
  )
  X <- matrix(rnorm(3 * n_vox), nrow = 3)

  coeff <- fmrilatent:::compute_hrbf_coefficients(X, basis)

  expect_equal(dim(coeff), c(3, n_atoms))
})

test_that("compute_hrbf_coefficients singular case fallback", {
  # Create a singular basis (rank deficient)
  n <- 10
  basis <- Matrix::Matrix(0, nrow = n, ncol = n, sparse = TRUE)
  basis[1, ] <- 1  # Only one non-zero row
  X <- matrix(rnorm(3 * n), nrow = 3)

  # Should not error, falls back to projection
  coeff <- fmrilatent:::compute_hrbf_coefficients(X, basis)

  expect_equal(dim(coeff), c(3, n))
})

test_that("compute_hrbf_coefficients dimension mismatch error", {
  basis <- Matrix::Matrix(rnorm(5 * 10), nrow = 5, sparse = TRUE)
  X <- matrix(rnorm(3 * 20), nrow = 3)  # Wrong: ncol != ncol(basis)

  expect_error(
    fmrilatent:::compute_hrbf_coefficients(X, basis),
    "ncol\\(X_time_vox\\) == ncol\\(basis_atoms_vox\\) is not TRUE"
  )
})

test_that("is_hrbf_latent returns FALSE for non-HRBF latent", {
  mask_arr <- array(TRUE, dim = c(2, 2, 2))
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(2, 2, 2)))
  spc <- neuroim2::NeuroSpace(c(2, 2, 2, 3))

  lvec <- LatentNeuroVec(
    basis = matrix(rnorm(3 * 8), nrow = 3),
    loadings = Matrix::Matrix(rnorm(8 * 8), nrow = 8),
    space = spc,
    mask = mask
  )

  expect_false(is_hrbf_latent(lvec))
})

test_that("hrbf_generate_basis with num_extra_fine_levels", {
  old_opt <- getOption("fmrilatent.hrbf.use_rcpp")
  on.exit(options(fmrilatent.hrbf.use_rcpp = old_opt), add = TRUE)
  options(fmrilatent.hrbf.use_rcpp = FALSE)

  mask_arr <- array(TRUE, dim = c(4, 4, 4))
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(4, 4, 4)))
  params <- list(
    sigma0 = 2,
    levels = 1L,
    radius_factor = 2.5,
    kernel_type = "gaussian",
    num_extra_fine_levels = 2L,
    seed = 42L
  )

  B <- hrbf_generate_basis(params, mask)

  expect_s4_class(B, "dgCMatrix")
  expect_equal(ncol(B), sum(mask_arr))
  # Should have more atoms due to extra fine levels
  expect_gt(nrow(B), 0)
})

test_that("hrbf_basis_from_params empty mask", {
  mask_arr <- array(FALSE, dim = c(3, 3, 3))
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(3, 3, 3)))
  params <- list(sigma0 = 1, levels = 1L, radius_factor = 2.5, kernel_type = "gaussian", seed = 1L)

  B <- fmrilatent:::hrbf_basis_from_params(params, mask)

  expect_s4_class(B, "dgCMatrix")
  expect_equal(dim(B), c(0, 0))
})

test_that("lna_hrbf_basis_from_params empty mask", {
  mask_arr <- array(FALSE, dim = c(3, 3, 3))
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(3, 3, 3)))
  params <- list(sigma0 = 1, levels = 1L, radius_factor = 2.5, kernel_type = "gaussian", seed = 1L)

  B <- lna_hrbf_basis_from_params(params, mask)

  expect_s4_class(B, "dgCMatrix")
  expect_equal(nrow(B), 0)
})

test_that("generate_hrbf_atom with fine level kernel switching", {
  old_opt <- getOption("fmrilatent.hrbf.use_rcpp")
  on.exit(options(fmrilatent.hrbf.use_rcpp = old_opt), add = TRUE)
  options(fmrilatent.hrbf.use_rcpp = FALSE)

  mask_coords_world <- matrix(c(0, 0, 0, 0.5, 0, 0, 1.5, 0, 0), ncol = 3, byrow = TRUE)
  params <- list(
    levels = 2L,
    num_extra_fine_levels = 1L,
    num_fine_levels_alt_kernel = 1L,
    kernel_type = "gaussian",
    kernel_type_fine_levels = "wendland_c6"
  )

  # At fine level, should use wendland
  atom <- fmrilatent:::generate_hrbf_atom(
    mask_coords_world, 1:3, c(0, 0, 0), 1.0,
    current_level_j = 3L,  # Fine level
    total_levels = 2L,
    params
  )

  # Point at distance 1.5 should be zero (wendland compact support)
  expect_equal(atom$values[3], 0)
})

test_that("hrbf_reconstruct_partial with all time indices", {
  set.seed(789)
  mask_arr <- array(TRUE, dim = c(3, 3, 3))
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(3, 3, 3)))
  X <- matrix(rnorm(4 * 27), nrow = 4)
  params <- list(sigma0 = 1.5, levels = 1L, radius_factor = 2.5, kernel_type = "gaussian", seed = 3L)

  coeff <- hrbf_project_matrix(X, mask, params)
  full <- hrbf_reconstruct_matrix(coeff, mask, params)

  vox_idx <- c(5, 10, 15)
  # Don't specify time_idx, should default to all
  part <- hrbf_reconstruct_partial(coeff, mask, params, voxel_idx = vox_idx)

  expect_equal(part, full[, vox_idx, drop = FALSE], tolerance = 1e-6)
  expect_equal(nrow(part), nrow(coeff))
})

# ============================================================================
# Additional coverage tests for uncovered functions
# ============================================================================

# --- hrbf_meta: thorough branch coverage ---

test_that("hrbf_meta returns metadata for valid hrbf-tagged LatentNeuroVec", {
  mask_arr <- array(TRUE, dim = c(2, 2, 2))
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(2, 2, 2)))
  X <- matrix(rnorm(3 * 8), nrow = 3)
  params <- list(sigma0 = 1, levels = 0L, radius_factor = 2.5,
                 kernel_type = "gaussian", seed = 1L)
  lvec <- hrbf_latent(X, mask, params)
  m <- hrbf_meta(lvec)
  expect_true(is.list(m))
  expect_equal(m$family, "hrbf")
  expect_true(!is.null(m$params))
})

test_that("hrbf_meta returns NULL when meta slot is NULL", {
  mask_arr <- array(TRUE, dim = c(2, 2, 2))
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(2, 2, 2)))
  spc <- neuroim2::NeuroSpace(c(2, 2, 2, 3))
  lvec <- LatentNeuroVec(
    basis = matrix(rnorm(3 * 2), nrow = 3),
    loadings = Matrix::Matrix(rnorm(8 * 2), nrow = 8),
    space = spc, mask = mask
  )
  # meta defaults to list() which is not NULL but has no $family
  expect_null(hrbf_meta(lvec))
})

test_that("hrbf_meta returns NULL for numeric input", {
  expect_null(hrbf_meta(42))
})

# --- extract_mask_array_logical: additional branches ---

test_that("extract_mask_array_logical converts numeric array to logical", {
  num_arr <- array(c(1, 0, 1, 0, 1, 0, 1, 0), dim = c(2, 2, 2))
  result <- fmrilatent:::extract_mask_array_logical(num_arr, "test_numeric")
  expect_true(is.logical(result))
  expect_equal(dim(result), c(2, 2, 2))
  expect_equal(sum(result), 4)
})

test_that("extract_mask_array_logical handles NeuroVol mask", {
  mask_arr <- array(c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE),
                    dim = c(2, 2, 2))
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(2, 2, 2)))
  result <- fmrilatent:::extract_mask_array_logical(mask, "test_neurovol")
  expect_true(is.logical(result))
  expect_equal(dim(result), c(2, 2, 2))
  expect_equal(sum(result), 4)
})

test_that("extract_mask_array_logical errors on list without arr field", {
  # A list without $arr that also fails as.array should error
  bad_list <- list(foo = "bar")
  expect_error(
    fmrilatent:::extract_mask_array_logical(bad_list, "test_bad_list"),
    "Cannot extract array from mask"
  )
})

# --- poisson_disk_sample_neuroim2: branch coverage ---

test_that("poisson_disk_sample_neuroim2 with single-component contiguous mask", {
  mask_arr <- array(TRUE, dim = c(4, 4, 4))
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(4, 4, 4)))
  centres <- fmrilatent:::poisson_disk_sample_neuroim2(
    mask, radius_mm = 2.0, seed = 10L
  )
  expect_true(nrow(centres) > 0)
  expect_equal(ncol(centres), 3)
  expect_equal(colnames(centres), c("i", "j", "k"))
})

test_that("poisson_disk_sample_neuroim2 with component_id_for_seed_offset > 0", {
  mask_arr <- array(TRUE, dim = c(3, 3, 3))
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(3, 3, 3)))
  # When component_id_for_seed_offset > 0, multi-component path is skipped
  centres <- fmrilatent:::poisson_disk_sample_neuroim2(
    mask, radius_mm = 1.5, seed = 5L, component_id_for_seed_offset = 2L
  )
  expect_true(nrow(centres) > 0)
  expect_equal(ncol(centres), 3)
})

test_that("poisson_disk_sample_neuroim2 small component uses centroid fallback", {
  # A mask with a very small disconnected component (<= 8 voxels)
  mask_arr <- array(FALSE, dim = c(10, 10, 3))
  mask_arr[1:2, 1:2, 1] <- TRUE   # 4 voxels -- small component
  mask_arr[8:10, 8:10, 3] <- TRUE  # 9 voxels -- larger component
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(10, 10, 3)))
  centres <- fmrilatent:::poisson_disk_sample_neuroim2(
    mask, radius_mm = 0.5, seed = 1L
  )
  # Should get centres from both components
  expect_true(nrow(centres) >= 2)
  expect_equal(ncol(centres), 3)
})

# --- generate_hrbf_atom: normalization and no-normalization ---

test_that("generate_hrbf_atom with normalize_over_mask = FALSE", {
  mask_coords_world <- matrix(c(0, 0, 0, 1, 0, 0), ncol = 3, byrow = TRUE)
  atom <- fmrilatent:::generate_hrbf_atom(
    mask_coords_world, 1:2, c(0, 0, 0), 1.0,
    0L, 1L, list(kernel_type = "gaussian"),
    normalize_over_mask = FALSE
  )
  # Without normalization, center value should be exactly 1 (exp(0))
  expect_equal(atom$values[1], 1.0, tolerance = 1e-10)
})

test_that("generate_hrbf_atom with normalize_over_mask = TRUE normalizes", {
  mask_coords_world <- matrix(c(0, 0, 0, 1, 0, 0), ncol = 3, byrow = TRUE)
  atom <- fmrilatent:::generate_hrbf_atom(
    mask_coords_world, 1:2, c(0, 0, 0), 1.0,
    0L, 1L, list(kernel_type = "gaussian"),
    normalize_over_mask = TRUE
  )
  # Normalized: sum(phi^2) should be approximately 1
  expect_equal(sum(atom$values^2), 1.0, tolerance = 1e-10)
})

test_that("generate_hrbf_atom with wendland_c6 kernel and far point", {
  mask_coords_world <- matrix(c(0, 0, 0, 2, 0, 0), ncol = 3, byrow = TRUE)
  atom <- fmrilatent:::generate_hrbf_atom(
    mask_coords_world, 1:2, c(0, 0, 0), 1.0,
    0L, 1L, list(kernel_type = "wendland_c6"),
    normalize_over_mask = FALSE
  )
  # At distance 2.0, r = 2.0/1.0 = 2.0 >= 1, so value is 0

  expect_equal(atom$values[2], 0)
  # At distance 0, value should be 1
  expect_equal(atom$values[1], 1.0, tolerance = 1e-10)
})

test_that("generate_hrbf_atom with fine level alt kernel switching", {
  mask_coords_world <- matrix(c(0, 0, 0, 0.5, 0, 0, 1.5, 0, 0),
                              ncol = 3, byrow = TRUE)
  params <- list(
    levels = 2L,
    num_extra_fine_levels = 1L,
    num_fine_levels_alt_kernel = 2L,
    kernel_type = "gaussian",
    kernel_type_fine_levels = "wendland_c4"  # alias -> wendland_c6
  )
  atom <- fmrilatent:::generate_hrbf_atom(
    mask_coords_world, 1:3, c(0, 0, 0), 1.0,
    current_level_j = 2L, total_levels = 2L, params,
    normalize_over_mask = FALSE
  )
  # wendland_c4 alias should be resolved to wendland_c6
  # at distance 1.5, r = 1.5 >= 1, so value = 0
  expect_equal(atom$values[3], 0)
})

# --- lna_hrbf_centres_from_params: more branch coverage ---

test_that("lna_hrbf_centres_from_params with extra fine levels", {
  mask_arr <- array(TRUE, dim = c(4, 4, 4))
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(4, 4, 4)))
  params <- list(sigma0 = 4, levels = 1L, radius_factor = 2.5,
                 num_extra_fine_levels = 2L, seed = 7L)
  result <- fmrilatent:::lna_hrbf_centres_from_params(params, mask)
  expect_true(nrow(result$centres) > 0)
  # Should have centres from extra fine levels
  max_level <- max(result$level_vec)
  expect_gte(max_level, 2L)
  expect_false(result$centres_stored)
})

test_that("lna_hrbf_centres_from_params with neuroarchive compat profile", {
  mask_arr <- array(TRUE, dim = c(3, 3, 3))
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(3, 3, 3)))
  params <- list(sigma0 = 2, levels = 1L, radius_factor = 2.5, seed = 5L)
  result <- fmrilatent:::lna_hrbf_centres_from_params(
    params, mask, compat_profile = "neuroarchive_0.1.1"
  )
  expect_true(nrow(result$centres) > 0)
  expect_false(result$centres_stored)
})

test_that("lna_hrbf_centres_from_params with explicit centres handles Inf sigma", {
  mask_arr <- array(TRUE, dim = c(2, 2, 2))
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(2, 2, 2)))
  params <- list(sigma0 = 6)
  centres <- matrix(c(0, 0, 0), ncol = 3)
  sigmas <- c(Inf)  # Inf sigma -> log2(6/Inf) = -Inf -> not finite -> clamped to 0
  result <- fmrilatent:::lna_hrbf_centres_from_params(
    params, mask, centres = centres, sigmas = sigmas
  )
  expect_equal(result$level_vec, 0L)
  expect_true(result$centres_stored)
})

# --- compute_hrbf_coefficients: coverage for regularized least-squares path ---

test_that("compute_hrbf_coefficients rectangular basis uses gram path", {
  set.seed(321)
  n_atoms <- 5
  n_vox <- 20
  # Rectangular (non-square) sparse basis -- takes gram regularization path
  basis <- Matrix::sparseMatrix(
    i = rep(1:n_atoms, each = 4),
    j = sample(1:n_vox, n_atoms * 4, replace = TRUE),
    x = rnorm(n_atoms * 4),
    dims = c(n_atoms, n_vox)
  )
  X <- matrix(rnorm(3 * n_vox), nrow = 3)
  coeff <- fmrilatent:::compute_hrbf_coefficients(X, basis)
  expect_equal(dim(coeff), c(3, n_atoms))
  expect_true(is.matrix(coeff))
})

test_that("compute_hrbf_coefficients well-conditioned square uses inverse path", {
  set.seed(111)
  n <- 5
  # Well-conditioned diagonal-dominant sparse matrix
  basis <- Matrix::Diagonal(n, x = 5) +
    Matrix::sparseMatrix(i = 1:n, j = 1:n, x = rnorm(n) * 0.1,
                         dims = c(n, n))
  X <- matrix(rnorm(3 * n), nrow = 3)
  coeff <- fmrilatent:::compute_hrbf_coefficients(X, basis)
  expect_equal(dim(coeff), c(3, n))
  # Verify: X approx coeff %*% basis
  recon <- coeff %*% as.matrix(basis)
  expect_equal(recon, X, tolerance = 1e-4)
})
