library(testthat)

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

test_that("heat_wavelet_latent constructs LatentNeuroVec", {
  skip_if_not_installed("rgsp")
  mask <- array(TRUE, dim = c(2, 2, 2))
  X <- matrix(rnorm(5 * sum(mask)), nrow = 5)
  lv <- heat_wavelet_latent(X, mask, spec = basis_heat_wavelet(scales = 1, order = 8))
  expect_s4_class(lv, "LatentNeuroVec")
  expect_equal(nrow(loadings(lv)), sum(mask))
  expect_equal(nrow(basis(lv)), nrow(X))
})

