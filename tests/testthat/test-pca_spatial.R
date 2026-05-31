library(testthat)

# =============================================================================
# Tests for cluster-local PCA spatial encoder (lift + encode spec)
# =============================================================================

as.matrix.pca_bad_as_matrix <- function(x, ...) {
  stop("base matrix data was unnecessarily coerced", call. = FALSE)
}

test_that("lift(ClusterReduction, spec_pca) returns block-sparse loadings with singular values", {
  set.seed(1)
  dims <- c(2, 2, 1)
  mask_arr <- array(TRUE, dim = dims)
  mask_vol <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(dims))
  n_time <- 6L
  n_vox <- sum(mask_arr)

  X <- matrix(rnorm(n_time * n_vox), nrow = n_time)
  map <- as.integer(c(1L, 1L, 2L, 2L))
  red <- make_cluster_reduction(mask_vol, map)

  L <- lift(red, basis_pca(k = 2L), data = X, center = TRUE, backend = "svd")

  expect_true(inherits(L, "Matrix"))
  expect_equal(dim(L), c(n_vox, 4L))

  d <- attr(L, "fmrilatent.singular_values")
  expect_true(is.numeric(d))
  expect_length(d, ncol(L))
})

test_that("lift(ClusterReduction, spec_pca) does not coerce base matrix data", {
  set.seed(12)
  dims <- c(2, 2, 1)
  mask_arr <- array(TRUE, dim = dims)
  mask_vol <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(dims))
  X <- structure(
    matrix(rnorm(6L * sum(mask_arr)), nrow = 6L),
    class = c("pca_bad_as_matrix", "matrix")
  )
  red <- make_cluster_reduction(mask_vol, as.integer(c(1L, 1L, 2L, 2L)))

  L <- lift(red, basis_pca(k = 2L), data = X, center = TRUE, backend = "svd")

  expect_true(inherits(L, "Matrix"))
  expect_equal(nrow(L), sum(mask_arr))
})

test_that("lift(ClusterReduction, spec_pca) warns that whitening is caller-owned", {
  set.seed(11)
  dims <- c(2, 2, 1)
  mask_arr <- array(TRUE, dim = dims)
  mask_vol <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(dims))
  X <- matrix(rnorm(6L * sum(mask_arr)), nrow = 6L)
  red <- make_cluster_reduction(mask_vol, as.integer(c(1L, 1L, 2L, 2L)))

  expect_warning(
    lift(red, basis_pca(k = 2L, whiten = TRUE), data = X, backend = "svd"),
    class = "fmrilatent_warning_pca_whiten_ignored"
  )
})

test_that("encode spec_space_pca reconstructs exactly at full per-cluster rank (centered)", {
  set.seed(2)
  dims <- c(2, 2, 1)
  mask_arr <- array(TRUE, dim = dims)
  mask_vol <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(dims))
  n_time <- 5L
  n_vox <- sum(mask_arr)

  X <- matrix(rnorm(n_time * n_vox), nrow = n_time)
  X <- sweep(X, 2, seq_len(n_vox), "+") # add distinct voxel means

  map <- as.integer(c(1L, 1L, 2L, 2L))
  red <- make_cluster_reduction(mask_vol, map)

  lv <- encode(
    X,
    spec_space_pca(k = 2L, center = TRUE, whiten = FALSE, backend = "svd"),
    mask = mask_vol,
    reduction = red,
    materialize = "matrix"
  )

  rec <- as.matrix(lv)
  expect_equal(dim(rec), dim(X))
  expect_lt(max(abs(rec - X)), 1e-8)
  expect_length(offset(lv), n_vox)
})

test_that("encode spec_space_pca whitening preserves reconstruction", {
  set.seed(3)
  dims <- c(2, 2, 1)
  mask_arr <- array(TRUE, dim = dims)
  mask_vol <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(dims))
  n_time <- 6L
  n_vox <- sum(mask_arr)

  X <- matrix(rnorm(n_time * n_vox), nrow = n_time)
  X <- sweep(X, 2, seq_len(n_vox), "+")

  map <- as.integer(c(1L, 1L, 2L, 2L))
  red <- make_cluster_reduction(mask_vol, map)

  lv <- encode(
    X,
    spec_space_pca(k = 2L, center = TRUE, whiten = TRUE, backend = "svd"),
    mask = mask_vol,
    reduction = red,
    materialize = "matrix"
  )

  rec <- as.matrix(lv)
  expect_lt(max(abs(rec - X)), 1e-8)
})

test_that("cluster PCA rejects non-positive component counts", {
  expect_error(basis_pca(k = 0L), class = "fmrilatent_error_invalid_count")

  dims <- c(2, 2, 1)
  mask_arr <- array(TRUE, dim = dims)
  mask_vol <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(dims))
  X <- matrix(rnorm(5L * sum(mask_arr)), nrow = 5L)
  red <- make_cluster_reduction(mask_vol, as.integer(c(1L, 1L, 2L, 2L)))
  bad_spec <- structure(list(k = 0L, whiten = FALSE), class = "spec_pca")

  expect_error(
    lift(red, bad_spec, data = X, backend = "svd"),
    class = "fmrilatent_error_invalid_count"
  )
})
