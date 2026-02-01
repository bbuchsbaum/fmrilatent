library(testthat)

set.seed(1)

make_mask_cube <- function(nx = 2, ny = 2, nz = 2, fill = TRUE) {
  if (fill) {
    array(TRUE, dim = c(nx, ny, nz))
  } else {
    m <- array(FALSE, dim = c(nx, ny, nz))
    m[1, 1, 1] <- TRUE
    m
  }
}

test_that("precompute_haar_scalings counts voxels", {
  mask_full <- make_mask_cube(2, 2, 2, TRUE)
  s_full <- fmrilatent:::precompute_haar_scalings(mask_full, 1)
  expect_equal(length(s_full), 1L)
  expect_equal(s_full[[1]]$sqrt_nvalid, sqrt(8))

  mask_single <- make_mask_cube(2, 2, 2, FALSE)
  s_single <- fmrilatent:::precompute_haar_scalings(mask_single, 1)
  expect_equal(sum(round(s_single[[1]]$sqrt_nvalid^2)), 1)
})

test_that("single-level Haar roundtrip works (R path)", {
  mask <- make_mask_cube(2, 2, 2, TRUE)
  X <- matrix(rnorm(16), nrow = 2)
  opts <- options(fmrilatent.hwt.use_rcpp = FALSE)
  coeff <- fmrilatent:::perform_haar_lift_analysis(X, mask, levels = 1, z_order_seed = 42L)
  reco  <- fmrilatent:::perform_haar_lift_synthesis(coeff, mask, levels = 1, z_order_seed = 42L)
  options(opts)
  expect_equal(reco, X, tolerance = 1e-7)
})

test_that("multi-level Haar roundtrip works (R path)", {
  mask <- array(TRUE, dim = c(3, 3, 3))
  X <- matrix(rnorm(5 * sum(mask)), nrow = 5)
  opts <- options(fmrilatent.hwt.use_rcpp = FALSE)
  coeff <- fmrilatent:::perform_haar_lift_analysis(X, mask, levels = 2, z_order_seed = 42L)
  reco  <- fmrilatent:::perform_haar_lift_synthesis(coeff, mask, levels = 2, z_order_seed = 42L)
  options(opts)
  expect_equal(reco, X, tolerance = 1e-6)
})

test_that("ROI subsetting returns matching subset", {
  mask <- array(FALSE, dim = c(3, 3, 3))
  mask[1:2, 1:2, 1:2] <- TRUE
  mask[3, 3, 3] <- TRUE
  X <- matrix(rnorm(4 * sum(mask)), nrow = 4)

  fw <- haar_wavelet_forward(X, mask, levels = 2, z_seed = 42L)
  full <- haar_wavelet_inverse(fw, mask, levels = 2, z_seed = 42L)

  roi <- array(FALSE, dim = c(3, 3, 3))
  roi[1:2, 1:2, 1:2] <- TRUE
  roi_reco <- haar_wavelet_inverse(fw, mask, levels = 2, z_seed = 42L, roi_mask = roi)

  global_idx <- which(as.logical(mask))
  roi_global <- which(as.logical(roi))
  col_keep <- which(global_idx %in% roi_global)

  expect_equal(roi_reco, full[, col_keep, drop = FALSE], tolerance = 1e-6)
})

test_that("levels_keep zeros selected detail levels", {
  mask <- array(TRUE, dim = c(2, 2, 2))
  X <- matrix(rnorm(12 * sum(mask)), nrow = 12)
  fw <- haar_wavelet_forward(X, mask, levels = 1, z_seed = 42L)

  full <- haar_wavelet_inverse(fw, mask)
  root_only <- haar_wavelet_inverse(fw, mask, levels_keep = integer(0))

  expect_equal(full, X, tolerance = 1e-6)
  expect_gt(sum(abs(full - root_only)), 1e-8)
})

test_that("ImplicitLatent predict matches HaarLatent predict", {
  mask <- array(TRUE, dim = c(2, 2, 2))
  X <- matrix(rnorm(6 * sum(mask)), nrow = 6)
  hl <- haar_latent(X, mask, levels = 1, z_seed = 42L)
  expect_true(is_implicit_latent(hl))
  expect_true(is_haar_latent(hl))
  full1 <- predict(hl)
  full2 <- predict.ImplicitLatent(hl)
  expect_equal(full1, full2, tolerance = 1e-6)
})
