test_that("basis_mat materializes DCT handle and caches in registry", {
  n_time <- 5L
  k <- 3L
  bh <- dct_basis_handle(n_time = n_time, k = k)

  mat_handle <- basis_mat(bh)
  mat_exp <- build_dct_basis(n_time = n_time, k = k)

  expect_equal(mat_handle, mat_exp)
  expect_true(fmrilatent:::`.latent_has_matrix`(bh@id, type = "basis"))
})

test_that("bspline handle builds correct dimensions", {
  n_time <- 10L
  k <- 6L
  bh <- bspline_basis_handle(n_time = n_time, k = k, degree = 3L)

  mat <- basis_mat(bh)
  expect_equal(dim(mat), c(n_time, k))
  expect_true(all(is.finite(mat)))
})

test_that("slepian temporal handle matches DPSS basis", {
  n_time <- 12L
  tr <- 2
  bandwidth <- 0.08
  bh <- slepian_temporal_handle(n_time = n_time, tr = tr, bandwidth = bandwidth, backend = "tridiag")
  mat_handle <- as.matrix(basis_mat(bh))
  mat_exp <- dpss_time_basis(n_time, tr = tr, bandwidth = bandwidth, k = ncol(mat_handle), backend = "tridiag")
  expect_equal(mat_handle, mat_exp, tolerance = 1e-8)
})

test_that("LatentNeuroVec with handle-backed basis/loadings reconstructs correctly", {
  # small 2x2x1 space, 3 time points, 2 components
  mask_arr <- array(TRUE, dim = c(2, 2, 1))
  mask_vol <- LogicalNeuroVol(mask_arr, NeuroSpace(c(2, 2, 1)))
  space <- NeuroSpace(c(2, 2, 1, 3))

  b_handle <- dct_basis_handle(n_time = 3L, k = 2L)

  load_mat <- Matrix::Matrix(
    rbind(
      c(1.0, 0.0),
      c(0.5, 0.2),
      c(0.1, 0.8),
      c(0.0, 1.0)
    ),
    sparse = FALSE
  )

  l_handle <- new("LoadingsHandle",
    id = "explicit-loadings-test",
    dim = as.integer(dim(load_mat)),
    kind = "explicit",
    spec = list(matrix = load_mat),
    label = "explicit-loadings"
  )

  lvec <- LatentNeuroVec(
    basis = b_handle,
    loadings = l_handle,
    space = space,
    mask = mask_vol,
    offset = numeric(0),
    label = "test-handle"
  )

  expected_mat <- as.matrix(basis_mat(lvec) %*% t(loadings_mat(lvec)))
  expected_vec <- c(t(expected_mat)) # time-major (t fast) to match 4D linear order

  lin_idx <- seq_len(prod(dim(lvec)))
  observed_vec <- linear_access(lvec, lin_idx)

  expect_equal(observed_vec, expected_vec, tolerance = 1e-10)
})

test_that("linear_access returns correct subset without full reconstruction", {
  set.seed(2025)
  mask_arr <- array(TRUE, dim = c(3, 3, 2)) # 18 vox
  mask_vol <- LogicalNeuroVol(mask_arr, NeuroSpace(c(3, 3, 2)))
  n_time <- 5L
  k <- 3L

  basis <- Matrix::Matrix(matrix(rnorm(n_time * k), n_time, k), sparse = FALSE)
  loadings <- Matrix::Matrix(matrix(rnorm(sum(mask_arr) * k), sum(mask_arr), k), sparse = FALSE)

  spc <- neuroim2::NeuroSpace(c(3, 3, 2, n_time))
  lvec <- LatentNeuroVec(
    basis = basis,
    loadings = loadings,
    space = spc,
    mask = mask_vol,
    offset = numeric(0),
    label = "subset-test"
  )

  full_arr <- as.array(lvec)
  total_len <- length(full_arr)

  idx_subset <- sample(total_len, size = 10L)
  expect_equal(linear_access(lvec, idx_subset), full_arr[idx_subset], tolerance = 1e-10)

  # also exercise matricized_access with (time, spatial) pairs
  spatial_total <- prod(dim(lvec)[1:3])
  t_idx <- sample(n_time, size = 4, replace = TRUE)
  s_idx <- sample(spatial_total, size = 4, replace = TRUE)
  pair_idx <- cbind(t_idx, s_idx)
  mat_vals <- matricized_access(lvec, pair_idx)
  # map to full_arr positions
  lin_idx <- (s_idx - 1L) + (t_idx - 1L) * spatial_total + 1L
  expect_equal(mat_vals, full_arr[lin_idx], tolerance = 1e-10)
})
