library(testthat)

test_that("diffusion wavelet lift builds dense loadings", {
  skip_if_not_installed("rgsp")
  set.seed(42)
  mask <- array(TRUE, dim = c(3, 3, 1))
  map <- seq_len(sum(mask))
  red <- make_cluster_reduction(mask, map)
  spec <- basis_diffusion_wavelet(target_rank = 3L, oversample = 2L, threshold = 0, max_scales = 1L)
  loadings <- lift(red, spec, k_neighbors = 2L)
  expect_s4_class(loadings, "dgeMatrix")
  expect_equal(nrow(loadings), length(map))
  expect_gt(ncol(loadings), 0L)
})

test_that("diffusion wavelet cluster loadings are near-orthonormal", {
  skip_if_not_installed("rgsp")
  set.seed(24)
  mask <- array(TRUE, dim = c(3, 3, 1))
  map <- seq_len(sum(mask))
  red <- make_cluster_reduction(mask, map)
  spec <- basis_diffusion_wavelet(target_rank = 3L, oversample = 2L, threshold = 0, max_scales = 1L)
  T_sparse <- fmrilatent:::build_cluster_graph(red, k_neighbors = 2L)
  cluster_loadings <- fmrilatent:::diffusion_wavelet_loadings(T_sparse, spec)
  gram <- crossprod(cluster_loadings)
  expect_lt(max(abs(gram - diag(ncol(cluster_loadings)))), 1e-5)
})

test_that("diffusion_wavelet_latent constructs LatentNeuroVec", {
  skip_if_not_installed("rgsp")
  set.seed(7)
  mask <- array(TRUE, dim = c(3, 3, 1))
  X <- matrix(rnorm(4 * sum(mask)), nrow = 4)
  lv <- diffusion_wavelet_latent(X, mask, spec = basis_diffusion_wavelet(target_rank = 2L, oversample = 1L, max_scales = 1L))
  expect_s4_class(lv, "LatentNeuroVec")
  expect_equal(nrow(loadings(lv)), sum(mask))
  expect_equal(nrow(basis(lv)), nrow(X))
})

test_that("diffusion wavelet loadings handle materializes and reconstructs", {
  skip_if_not_installed("rgsp")
  set.seed(123)
  mask <- array(TRUE, dim = c(2, 2, 1))
  map <- seq_len(sum(mask))
  red <- make_cluster_reduction(mask, map)
  spec <- basis_diffusion_wavelet(target_rank = 2L, oversample = 1L, threshold = 0, max_scales = 1L)

  l_handle <- diffusion_wavelet_loadings_handle(red, spec, k_neighbors = 2L, id = "dw-handle-test")
  expect_s4_class(l_handle, "LoadingsHandle")
  expect_true(fmrilatent:::`.latent_has_matrix`(l_handle@id, type = "loadings"))

  load_mat <- loadings_mat(l_handle)
  n_time <- 3L
  basis_mat <- Matrix::Matrix(matrix(rnorm(n_time * ncol(load_mat)), n_time, ncol(load_mat)), sparse = FALSE)
  spc <- neuroim2::NeuroSpace(c(dim(mask), n_time))
  mask_vol <- LogicalNeuroVol(mask, neuroim2::NeuroSpace(dim(mask)))

  lvec <- LatentNeuroVec(
    basis = basis_mat,
    loadings = l_handle,
    space = spc,
    mask = mask_vol,
    offset = numeric(0),
    label = "dw-handle"
  )

  expected <- as.matrix(basis_mat %*% t(load_mat))
  expect_equal(as.matrix(lvec), expected)
})

test_that("diffusion wavelet loadings handle defaults to diffusion basis spec", {
  skip_if_not_installed("rgsp")
  mask <- array(TRUE, dim = c(2, 2, 1))
  map <- seq_len(sum(mask))
  red <- make_cluster_reduction(mask, map)

  handle <- diffusion_wavelet_loadings_handle(red, k_neighbors = 2L)

  expect_s4_class(handle, "LoadingsHandle")
  expect_s3_class(handle@spec$basis_spec, "spec_diffusion_wavelet")
  expect_identical(handle@spec$basis_spec, basis_diffusion_wavelet())
})

test_that("diffusion wavelet specs are deterministic across RNG state", {
  skip_if_not_installed("rgsp")
  mask <- array(TRUE, dim = c(3, 3, 1))
  map <- seq_len(sum(mask))
  red <- make_cluster_reduction(mask, map)
  spec <- basis_diffusion_wavelet(
    target_rank = 2L,
    oversample = 1L,
    threshold = 0,
    max_scales = 1L,
    seed = 99L
  )

  set.seed(1)
  loadings_1 <- lift(red, spec, k_neighbors = 2L)
  set.seed(999)
  loadings_2 <- lift(red, spec, k_neighbors = 2L)

  expect_equal(as.matrix(loadings_1), as.matrix(loadings_2), tolerance = 1e-12)
})

test_that("diffusion wavelet handles derive deterministic ids", {
  skip_if_not_installed("rgsp")
  mask <- array(TRUE, dim = c(3, 3, 1))
  map <- seq_len(sum(mask))
  red <- make_cluster_reduction(mask, map)
  spec <- basis_diffusion_wavelet(
    target_rank = 2L,
    oversample = 1L,
    threshold = 0,
    max_scales = 1L,
    seed = 99L
  )

  handle_1 <- diffusion_wavelet_loadings_handle(red, spec, k_neighbors = 2L)
  handle_2 <- diffusion_wavelet_loadings_handle(red, spec, k_neighbors = 2L)

  expect_equal(handle_1@id, handle_2@id)
})

test_that("diffusion wavelet specs reject non-positive and non-finite scalars", {
  expect_error(
    basis_diffusion_wavelet(max_scales = 0L),
    class = "fmrilatent_error_invalid_count"
  )
  expect_error(
    basis_diffusion_wavelet(target_rank = Inf),
    class = "fmrilatent_error_invalid_count"
  )
  expect_error(
    basis_diffusion_wavelet(threshold = NA_real_),
    class = "fmrilatent_error_invalid_scalar"
  )
})

test_that("diffusion wavelet spec supports sparsify_eps with threshold compatibility", {
  spec <- basis_diffusion_wavelet(sparsify_eps = 0, max_scales = 1L)
  expect_equal(spec$sparsify_eps, 0)
  expect_equal(spec$threshold, 0)

  spec_old <- basis_diffusion_wavelet(threshold = 0, max_scales = 1L)
  expect_equal(spec_old$sparsify_eps, 0)
  expect_equal(spec_old$threshold, 0)

  expect_error(
    basis_diffusion_wavelet(sparsify_eps = 0, threshold = 1e-5),
    class = "fmrilatent_error_invalid_argument"
  )
})

test_that("diffusion wavelet spec preserves legacy positional threshold and max_scales", {
  spec <- basis_diffusion_wavelet(10L, 2L, 0, 1L)
  expect_equal(spec$target_rank, 10L)
  expect_equal(spec$oversample, 2L)
  expect_equal(spec$threshold, 0)
  expect_equal(spec$sparsify_eps, 0)
  expect_equal(spec$max_scales, 1L)
})

test_that("diffusion wavelet loadings validate manually assembled specs", {
  T_sparse <- Matrix::Diagonal(3)
  bad_spec <- structure(
    list(target_rank = 1L, oversample = 0L, threshold = 0, max_scales = 0L),
    class = "spec_diffusion_wavelet"
  )
  expect_error(
    fmrilatent:::diffusion_wavelet_loadings(T_sparse, bad_spec),
    class = "fmrilatent_error_invalid_count"
  )
})

test_that("diffusion step uses Galerkin-compressed operator", {
  T_mat <- matrix(c(0.8, 0.1, 0.2, 0.7), nrow = 2L)
  T_op <- function(X) T_mat %*% X

  set.seed(123)
  step <- fmrilatent:::randomized_diffusion_step(
    T_op,
    n = 2L,
    target_rank = 1L,
    oversample = 0L,
    threshold = 0
  )

  expect_equal(
    as.matrix(step$T_compressed),
    as.matrix(crossprod(step$Q, T_op(step$Q))),
    tolerance = 1e-12
  )
})

test_that("diffusion wavelet loadings accept function operators with explicit n", {
  T_mat <- Matrix::Matrix(matrix(c(0.8, 0.1, 0.2, 0.7), nrow = 2L), sparse = TRUE)
  spec <- basis_diffusion_wavelet(
    target_rank = 1L,
    oversample = 0L,
    threshold = 0,
    max_scales = 1L,
    seed = 42L
  )
  T_fun <- function(X) T_mat %*% X

  from_matrix <- fmrilatent:::diffusion_wavelet_loadings(T_mat, spec)
  from_function <- fmrilatent:::diffusion_wavelet_loadings(T_fun, spec, n = nrow(T_mat))

  expect_equal(as.matrix(from_function), as.matrix(from_matrix), tolerance = 1e-12)
})

test_that("diffusion wavelet seed NULL leaves RNG un-restored", {
  T_mat <- Matrix::Diagonal(3)
  spec <- basis_diffusion_wavelet(
    target_rank = 1L,
    oversample = 0L,
    threshold = 0,
    max_scales = 1L,
    seed = NULL
  )

  set.seed(99)
  before <- .Random.seed
  fmrilatent:::diffusion_wavelet_loadings(T_mat, spec)
  expect_false(identical(.Random.seed, before))
})
