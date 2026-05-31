library(testthat)

test_that("heat wavelet specs validate numeric controls at construction", {
  expect_error(
    basis_heat_wavelet(scales = c(1, 0)),
    "positive finite"
  )
  expect_error(
    basis_heat_wavelet(scales = c(1, NA_real_)),
    "positive finite"
  )
  expect_error(
    basis_heat_wavelet(order = 0L),
    "positive integer"
  )
  expect_error(
    basis_heat_wavelet(threshold = -1e-6),
    "non-negative finite"
  )

  spec <- basis_heat_wavelet(scales = c(1, 2), order = 3L, threshold = 0)
  expect_equal(spec$scales, c(1, 2))
  expect_equal(spec$order, 3L)
  expect_equal(spec$threshold, 0)
})

test_that("heat wavelet lift builds sparse loadings", {
  skip_if_not_installed("rgsp")
  mask <- array(TRUE, dim = c(2, 2, 2))
  map <- seq_len(sum(mask))
  red <- make_cluster_reduction(mask, map)
  spec <- basis_heat_wavelet(scales = c(1, 2), order = 10, threshold = 0)
  loadings <- lift(red, spec, k_neighbors = 3L)
  expect_s4_class(loadings, "dgCMatrix")
  expect_equal(nrow(loadings), length(map))
  expect_equal(ncol(loadings), length(map) * length(spec$scales))
})

test_that("heat wavelet lift clamps k_neighbors for two-voxel clusters", {
  skip_if_not_installed("rgsp")
  mask <- array(TRUE, dim = c(2, 1, 1))
  map <- c(1L, 1L)
  red <- make_cluster_reduction(mask, map)
  spec <- basis_heat_wavelet(scales = c(1, 2), order = 10, threshold = 0)

  loadings <- lift(red, spec, k_neighbors = 6L)

  expect_s4_class(loadings, "dgCMatrix")
  expect_equal(dim(loadings), c(2L, 4L))
})

test_that("heat wavelet lift guards dense identity extraction by cluster size", {
  skip_if_not_installed("rgsp")
  mask <- array(TRUE, dim = c(3, 1, 1))
  map <- c(1L, 1L, 1L)
  red <- make_cluster_reduction(mask, map)
  spec <- basis_heat_wavelet(scales = 1, order = 8, threshold = 0)
  old <- options(fmrilatent.heat_wavelet.max_dense_cluster = 2L)
  on.exit(options(old), add = TRUE)

  expect_error(
    lift(red, spec, k_neighbors = 2L),
    class = "fmrilatent_error_dense_cluster_limit"
  )
})

test_that("heat_wavelet_latent defaults to a single global cluster", {
  skip_if_not_installed("rgsp")
  mask <- array(TRUE, dim = c(2, 1, 1))
  X <- matrix(rnorm(4 * sum(mask)), nrow = 4)
  spec <- basis_heat_wavelet(scales = 1, order = 10, threshold = 0)

  lv <- heat_wavelet_latent(X, mask, spec = spec, k_neighbors = 6L)
  global_red <- make_cluster_reduction(mask, rep(1L, sum(mask)))
  singleton_red <- make_cluster_reduction(mask, seq_len(sum(mask)))
  global_loadings <- lift(global_red, spec, k_neighbors = 6L)
  singleton_loadings <- lift(singleton_red, spec, k_neighbors = 6L)

  expect_equal(as.matrix(loadings(lv)), as.matrix(global_loadings), tolerance = 1e-8)
  expect_false(isTRUE(all.equal(
    as.matrix(loadings(lv)),
    as.matrix(singleton_loadings),
    tolerance = 1e-8
  )))
})

test_that("heat_wavelet_latent constructs LatentNeuroVec", {
  skip_if_not_installed("rgsp")
  mask <- array(TRUE, dim = c(2, 2, 2))
  X <- matrix(rnorm(5 * sum(mask)), nrow = 5)
  lv <- heat_wavelet_latent(X, mask, spec = basis_heat_wavelet(scales = 1, order = 8))
  expect_s4_class(lv, "LatentNeuroVec")
  expect_equal(nrow(loadings(lv)), sum(mask))
  expect_equal(nrow(basis(lv)), nrow(X))
})
