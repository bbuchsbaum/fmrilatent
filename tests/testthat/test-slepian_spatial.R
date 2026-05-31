test_that("slepian spatial lift builds sparse loadings", {
  skip_if_not_installed("RSpectra")
  mask <- array(TRUE, dim = c(2, 2, 2))
  map <- seq_len(sum(mask))
  red <- make_cluster_reduction(mask, map)
  spec <- basis_slepian(k = 2)
  L <- lift(red, spec, k_neighbors = 3L)
  expect_s4_class(L, "dgCMatrix")
  expect_equal(nrow(L), length(map))
  expect_gt(ncol(L), 0)
})

test_that("slepian_spatial_latent constructs LatentNeuroVec", {
  skip_if_not_installed("RSpectra")
  mask <- array(TRUE, dim = c(2, 2, 2))
  X <- matrix(rnorm(5 * sum(mask)), nrow = 5)
  lv <- slepian_spatial_latent(X, mask, spec = basis_slepian(k = 1))
  expect_s4_class(lv, "LatentNeuroVec")
  expect_equal(nrow(loadings(lv)), sum(mask))
  expect_equal(nrow(basis(lv)), nrow(X))
})

test_that("slepian spatial rejects non-positive component counts", {
  expect_error(basis_slepian(k = 0L), class = "fmrilatent_error_invalid_count")

  mask <- array(TRUE, dim = c(2, 2, 1))
  red <- make_cluster_reduction(mask, seq_len(sum(mask)))
  bad_spec <- structure(list(k = 0L, type = "laplacian"), class = "spec_slepian")
  expect_error(
    lift(red, bad_spec, k_neighbors = 2L),
    class = "fmrilatent_error_invalid_count"
  )
})

test_that("slepian spatial uses base eigen when rank equals cluster size", {
  skip_if_not_installed("rgsp")
  mask <- array(TRUE, dim = c(2, 2, 1))
  red <- make_cluster_reduction(mask, rep(1L, sum(mask)))
  L <- lift(red, basis_slepian(k = sum(mask)), k_neighbors = 2L)
  expect_equal(dim(L), c(sum(mask), sum(mask)))
  expect_true(all(is.finite(as.matrix(L))))
})
