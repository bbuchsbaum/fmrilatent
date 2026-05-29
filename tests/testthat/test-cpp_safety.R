test_that("active_pencil_wavelet validates dims length", {
  expect_error(
    active_pencil_wavelet(
      data_voxels = c(1),
      coords = matrix(c(1L, 1L, 1L), nrow = 1L),
      dims = as.integer(c(1L, 1L)),
      levels = 1L,
      forward = TRUE
    ),
    "dims must have length 3"
  )
})

test_that("Haar Rcpp entry points validate levels", {
  mask <- array(TRUE, dim = c(2, 2, 2))
  scalings <- precompute_haar_scalings_rcpp(mask, 1L)

  expect_error(
    precompute_haar_scalings_rcpp(mask, -1L),
    "levels must be non-negative"
  )
  expect_error(
    forward_lift_rcpp(
      data_morton = as.numeric(seq_len(8)),
      mask_flat_morton = rep(TRUE, 8),
      mask_dims = as.integer(c(2L, 2L, 2L)),
      levels = 2L,
      scalings = scalings
    ),
    "levels cannot exceed length"
  )
  expect_error(
    inverse_lift_rcpp(
      root_coeff = 1,
      detail_vecs = list(rep(0, 8)),
      mask_flat_morton = rep(TRUE, 8),
      mask_dims = as.integer(c(2L, 2L, 2L)),
      levels = 2L,
      scalings = scalings
    ),
    "levels cannot exceed length"
  )
})

test_that("get_valid_finest_blocks_rcpp guards integer Morton code width", {
  mask <- array(TRUE, dim = c(2049, 1, 1))
  expect_error(
    get_valid_finest_blocks_rcpp(mask),
    "too large for 32-bit Morton codes"
  )
})
