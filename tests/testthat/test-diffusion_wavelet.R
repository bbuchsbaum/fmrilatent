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
