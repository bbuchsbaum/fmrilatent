test_that("active_pencil_wavelet validates dims length", {
  expect_error(
    active_pencil_wavelet(
      data_voxels = c(1),
      coords = matrix(c(1L, 1L, 1L), nrow = 1L),
      dims = as.integer(c(1L, 1L)),
      levels = 1L,
      forward = TRUE
    ),
    "dims must have length 3",
    class = "fmrilatent_error_cpp_boundary"
  )
})

test_that("active_pencil_wavelet rejects invalid levels and duplicate coordinates", {
  coords <- matrix(c(
    1L, 1L, 1L,
    1L, 1L, 1L
  ), ncol = 3L, byrow = TRUE)

  expect_error(
    active_pencil_wavelet(
      data_voxels = c(1, 2),
      coords = coords,
      dims = as.integer(c(2L, 2L, 2L)),
      levels = 0L,
      forward = TRUE
    ),
    "levels must be positive",
    class = "fmrilatent_error_cpp_boundary"
  )
  expect_error(
    active_pencil_wavelet(
      data_voxels = c(1, 2),
      coords = coords,
      dims = as.integer(c(2L, 2L, 2L)),
      levels = 1L,
      forward = TRUE
    ),
    "coords must be unique",
    class = "fmrilatent_error_cpp_boundary"
  )
})

test_that("Haar Rcpp entry points validate levels", {
  mask <- array(TRUE, dim = c(2, 2, 2))
  scalings <- precompute_haar_scalings_rcpp(mask, 1L)

  expect_error(
    precompute_haar_scalings_rcpp(mask, -1L),
    "levels must be non-negative",
    class = "fmrilatent_error_cpp_boundary"
  )
  expect_error(
    forward_lift_rcpp(
      data_morton = as.numeric(seq_len(8)),
      mask_flat_morton = rep(TRUE, 8),
      mask_dims = as.integer(c(2L, 2L, 2L)),
      levels = 2L,
      scalings = scalings
    ),
    "levels cannot exceed length",
    class = "fmrilatent_error_cpp_boundary"
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
    "levels cannot exceed length",
    class = "fmrilatent_error_cpp_boundary"
  )
})

test_that("Haar Rcpp entry points validate malformed direct-call lengths", {
  mask_flat <- rep(TRUE, 8)
  mask_dims <- as.integer(c(2L, 2L, 2L))
  scalings <- precompute_haar_scalings_rcpp(array(TRUE, dim = c(2, 2, 2)), 1L)

  expect_error(
    forward_lift_rcpp(
      data_morton = as.numeric(seq_len(7)),
      mask_flat_morton = mask_flat,
      mask_dims = mask_dims,
      levels = 1L,
      scalings = scalings
    ),
    "data_morton length must match active mask count",
    class = "fmrilatent_error_cpp_boundary"
  )
  expect_error(
    forward_lift_rcpp(
      data_morton = as.numeric(seq_len(8)),
      mask_flat_morton = rep(TRUE, 7),
      mask_dims = mask_dims,
      levels = 1L,
      scalings = scalings
    ),
    "mask_flat_morton length must match mask_dims product",
    class = "fmrilatent_error_cpp_boundary"
  )
  expect_error(
    inverse_lift_rcpp(
      root_coeff = c(1, 2),
      detail_vecs = list(rep(0, 8)),
      mask_flat_morton = mask_flat,
      mask_dims = mask_dims,
      levels = 1L,
      scalings = scalings
    ),
    "root_coeff length must match coarsest active count",
    class = "fmrilatent_error_cpp_boundary"
  )
})

test_that("get_valid_finest_blocks_rcpp guards integer Morton code width", {
  mask <- array(TRUE, dim = c(2049, 1, 1))
  expect_error(
    get_valid_finest_blocks_rcpp(mask),
    "too large for Morton codes \\(max 10 bits / 1024 per axis\\)",
    class = "fmrilatent_error_morton_overflow"
  )
})

test_that("HRBF and DPSS Rcpp wrappers translate validation errors", {
  expect_error(
    hrbf_atoms_rcpp(
      mask_xyz_world = matrix(0, nrow = 1L, ncol = 3L),
      centres_xyz_world = matrix(0, nrow = 1L, ncol = 3L),
      sigma_vec_mm = -1,
      kernel_type = "gaussian"
    ),
    "sigma values must be positive",
    class = "fmrilatent_error_cpp_boundary"
  )
  expect_error(
    hrbf_atoms_rcpp(
      mask_xyz_world = matrix(0, nrow = 1L, ncol = 3L),
      centres_xyz_world = matrix(0, nrow = 1L, ncol = 3L),
      sigma_vec_mm = 1,
      kernel_type = "not-a-kernel"
    ),
    "kernel_type must be one of",
    class = "fmrilatent_error_cpp_boundary"
  )
  expect_error(
    generate_dpss_basis_rcpp(0L, 1, 1L),
    "n must be positive",
    class = "fmrilatent_error_cpp_boundary"
  )
  expect_error(
    generate_dpss_tridiag_rcpp(0L, 1, 1L),
    "n must be positive",
    class = "fmrilatent_error_cpp_boundary"
  )
  expect_error(
    generate_dpss_basis_rcpp(4097L, 1, 1L),
    "dense-allocation safety limit",
    class = "fmrilatent_error_cpp_boundary"
  )
  expect_error(
    generate_dpss_tridiag_rcpp(4097L, 1, 1L),
    "dense-allocation safety limit",
    class = "fmrilatent_error_cpp_boundary"
  )
})
