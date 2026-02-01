# Tests for searchlight_utils.R
# These functions compute Gram matrices and apply user functions over neighborhoods
# in latent space. They work with Matrix objects for basis and loadings.

# Helper to create Matrix objects for testing
make_matrix <- function(nrow, ncol, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  Matrix::Matrix(matrix(rnorm(nrow * ncol), nrow = nrow, ncol = ncol), sparse = FALSE)
}

# ------------------------------------------------------------------------------
# Tests for compute_local_gram()
# ------------------------------------------------------------------------------

test_that("compute_local_gram returns list of Gram matrices with simplify=FALSE", {
  set.seed(42)
  n_vox <- 20
  k <- 4
  L <- make_matrix(n_vox, k, seed = 42)

  neighborhoods <- list(
    1:5,
    6:10,
    11:15
  )

  result <- compute_local_gram(L, neighborhoods, simplify = FALSE)

  # Should return a list

  expect_type(result, "list")
  expect_length(result, 3)

  # Each element should be a k x k matrix
  for (i in seq_along(result)) {
    expect_equal(nrow(result[[i]]), k)
    expect_equal(ncol(result[[i]]), k)
  }

  # Verify Gram matrix is computed correctly for first neighborhood
  L_subset <- L[neighborhoods[[1]], , drop = FALSE]
  expected_gram <- Matrix::crossprod(L_subset)
  expect_equal(as.matrix(result[[1]]), as.matrix(expected_gram), tolerance = 1e-10)
})

test_that("compute_local_gram with simplify=TRUE returns array when all neighborhoods same size", {
  set.seed(123)
  n_vox <- 15
  k <- 3
  L <- make_matrix(n_vox, k, seed = 123)

  # All neighborhoods have same size (5 voxels each)
  neighborhoods <- list(
    1:5,
    6:10,
    11:15
  )

  result <- compute_local_gram(L, neighborhoods, simplify = TRUE)

  # Should return an array when simplify=TRUE and all same size
  expect_true(is.array(result))
  expect_equal(dim(result), c(k, k, 3))

  # Verify the array slices match expected Gram matrices
  for (i in seq_along(neighborhoods)) {
    L_sub <- L[neighborhoods[[i]], , drop = FALSE]
    expected <- Matrix::crossprod(L_sub)
    expect_equal(result[, , i], as.matrix(expected), tolerance = 1e-10)
  }
})

test_that("compute_local_gram simplify returns list when neighborhoods differ in size", {
  set.seed(456)
  n_vox <- 15
  k <- 3
  L <- make_matrix(n_vox, k, seed = 456)

  # Neighborhoods with different sizes
  neighborhoods <- list(
    1:5,      # 5 voxels
    6:8,      # 3 voxels
    9:15      # 7 voxels
  )

  # Even with simplify=TRUE, should return list because sizes differ
  result <- compute_local_gram(L, neighborhoods, simplify = TRUE)

  # The simplify condition checks if all Gram matrices have same dimensions.
  # Since Gram is always k x k regardless of neighborhood size, it should still
  # return an array (the Gram matrices are always k x k)
  expect_true(is.array(result))
  expect_equal(dim(result), c(k, k, 3))
})

test_that("compute_local_gram handles single voxel neighborhoods", {
  set.seed(789)
  n_vox <- 10
  k <- 2
  L <- make_matrix(n_vox, k, seed = 789)

  neighborhoods <- list(
    1L,           # Single voxel
    c(2L, 3L)     # Two voxels
  )

  result <- compute_local_gram(L, neighborhoods, simplify = FALSE)

  expect_length(result, 2)
  # Both should still produce k x k Gram matrices
  expect_equal(dim(result[[1]]), c(k, k))
  expect_equal(dim(result[[2]]), c(k, k))

  # Verify single voxel case: Gram is outer product of loading vector with itself
  L_single <- L[1, , drop = FALSE]
  expected_single <- Matrix::crossprod(L_single)
  expect_equal(as.matrix(result[[1]]), as.matrix(expected_single), tolerance = 1e-10)
})

test_that("compute_local_gram handles empty neighborhoods list", {
  set.seed(101)
  n_vox <- 10
  k <- 3
  L <- make_matrix(n_vox, k, seed = 101)

  neighborhoods <- list()

  result <- compute_local_gram(L, neighborhoods, simplify = FALSE)

  expect_type(result, "list")
  expect_length(result, 0)
})

test_that("compute_local_gram with simplify=FALSE returns list even when all same size", {
  set.seed(202)
  n_vox <- 9
  k <- 2
  L <- make_matrix(n_vox, k, seed = 202)

  neighborhoods <- list(
    1:3,
    4:6,
    7:9
  )

  result <- compute_local_gram(L, neighborhoods, simplify = FALSE)

  # Should always return list when simplify=FALSE
  expect_type(result, "list")
  expect_length(result, 3)
})

test_that("compute_local_gram produces symmetric Gram matrices", {
  set.seed(303)
  n_vox <- 12
  k <- 4
  L <- make_matrix(n_vox, k, seed = 303)

  neighborhoods <- list(1:4, 5:8, 9:12)

  result <- compute_local_gram(L, neighborhoods, simplify = FALSE)

  for (i in seq_along(result)) {
    M <- as.matrix(result[[i]])
    # Gram matrices should be symmetric
    expect_equal(M, t(M), tolerance = 1e-10)
    # Gram matrices should be positive semi-definite (all eigenvalues >= 0)
    eigenvalues <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
    expect_true(all(eigenvalues >= -1e-10))  # Allow small numerical error
  }
})

test_that("compute_local_gram handles overlapping neighborhoods", {
  set.seed(404)
  n_vox <- 10
  k <- 3
  L <- make_matrix(n_vox, k, seed = 404)

  # Overlapping neighborhoods
  neighborhoods <- list(
    1:5,
    3:7,   # Overlaps with first
    6:10   # Overlaps with second
  )

  result <- compute_local_gram(L, neighborhoods, simplify = FALSE)

  expect_length(result, 3)

  # Verify each is computed correctly despite overlap
  for (i in seq_along(neighborhoods)) {
    L_sub <- L[neighborhoods[[i]], , drop = FALSE]
    expected <- Matrix::crossprod(L_sub)
    expect_equal(as.matrix(result[[i]]), as.matrix(expected), tolerance = 1e-10)
  }
})

test_that("compute_local_gram works with k=1 (single component)", {
  set.seed(505)
  n_vox <- 10
  k <- 1
  L <- make_matrix(n_vox, k, seed = 505)

  neighborhoods <- list(1:3, 4:6, 7:10)

  result <- compute_local_gram(L, neighborhoods, simplify = TRUE)

  # Should be 1 x 1 x 3 array
  expect_equal(dim(result), c(1, 1, 3))

  # Each Gram matrix is just the sum of squared loadings
  for (i in seq_along(neighborhoods)) {
    L_sub <- L[neighborhoods[[i]], , drop = FALSE]
    expected <- sum(L_sub^2)
    expect_equal(result[1, 1, i], expected, tolerance = 1e-10)
  }
})

# ------------------------------------------------------------------------------
# Tests for latent_searchlight()
# ------------------------------------------------------------------------------

test_that("latent_searchlight applies function over neighborhoods", {
  set.seed(606)
  n_time <- 10
  n_vox <- 20
  k <- 4

  B <- make_matrix(n_time, k, seed = 606)
  L <- make_matrix(n_vox, k, seed = 607)

  neighborhoods <- list(
    1:5,
    6:10,
    11:15,
    16:20
  )

  # Function that captures all arguments and returns info about them
  test_fun <- function(B, L_V, M_V, idx, ...) {
    list(
      B_rows = nrow(B),
      B_cols = ncol(B),
      L_V_rows = nrow(L_V),
      L_V_cols = ncol(L_V),
      M_V_dim = dim(M_V),
      idx_len = length(idx),
      idx_vals = idx
    )
  }

  result <- latent_searchlight(B, L, neighborhoods, test_fun)

  expect_type(result, "list")
  expect_length(result, 4)

  # Check all results have correct structure
  for (i in seq_along(result)) {
    expect_equal(result[[i]]$B_rows, n_time)
    expect_equal(result[[i]]$B_cols, k)
    expect_equal(result[[i]]$L_V_rows, 5)  # Each neighborhood has 5 voxels
    expect_equal(result[[i]]$L_V_cols, k)
    expect_equal(result[[i]]$M_V_dim, c(k, k))
    expect_equal(result[[i]]$idx_len, 5)
    expect_equal(result[[i]]$idx_vals, neighborhoods[[i]])
  }
})

test_that("latent_searchlight passes extra arguments to function", {
  set.seed(707)
  n_time <- 5
  n_vox <- 10
  k <- 2

  B <- make_matrix(n_time, k, seed = 707)
  L <- make_matrix(n_vox, k, seed = 708)

  neighborhoods <- list(1:5, 6:10)

  # Function that uses extra arguments
  test_fun <- function(B, L_V, M_V, idx, multiplier, offset) {
    sum(L_V) * multiplier + offset
  }

  result <- latent_searchlight(B, L, neighborhoods, test_fun,
                               multiplier = 2, offset = 10)

  expect_length(result, 2)
  expect_type(result[[1]], "double")
  expect_type(result[[2]], "double")

  # Verify calculation for each neighborhood
  for (i in seq_along(neighborhoods)) {
    L_sub <- L[neighborhoods[[i]], , drop = FALSE]
    expected <- sum(L_sub) * 2 + 10
    expect_equal(result[[i]], expected, tolerance = 1e-10)
  }
})

test_that("latent_searchlight handles empty neighborhoods list", {
  set.seed(808)
  n_time <- 5
  n_vox <- 10
  k <- 2

  B <- make_matrix(n_time, k, seed = 808)
  L <- make_matrix(n_vox, k, seed = 809)

  neighborhoods <- list()

  test_fun <- function(B, L_V, M_V, idx, ...) {
    nrow(L_V)
  }

  result <- latent_searchlight(B, L, neighborhoods, test_fun)

  expect_type(result, "list")
  expect_length(result, 0)
})

test_that("latent_searchlight computes M_V correctly as Gram matrix", {
  set.seed(909)
  n_time <- 8
  n_vox <- 15
  k <- 3

  B <- make_matrix(n_time, k, seed = 909)
  L <- make_matrix(n_vox, k, seed = 910)

  neighborhoods <- list(1:5, 6:10, 11:15)

  # Function that returns the M_V matrix for verification
  test_fun <- function(B, L_V, M_V, idx, ...) {
    as.matrix(M_V)
  }

  result <- latent_searchlight(B, L, neighborhoods, test_fun)

  # Verify M_V matches crossprod(L_V) for each neighborhood
  for (i in seq_along(neighborhoods)) {
    L_sub <- L[neighborhoods[[i]], , drop = FALSE]
    expected_M <- as.matrix(Matrix::crossprod(L_sub))
    expect_equal(result[[i]], expected_M, tolerance = 1e-10)
  }
})

test_that("latent_searchlight provides correct B matrix to all neighborhoods", {
  set.seed(1010)
  n_time <- 6
  n_vox <- 12
  k <- 2

  B <- make_matrix(n_time, k, seed = 1010)
  L <- make_matrix(n_vox, k, seed = 1011)

  neighborhoods <- list(1:4, 5:8, 9:12)

  # Function that returns B matrix
  test_fun <- function(B, L_V, M_V, idx, ...) {
    as.matrix(B)
  }

  result <- latent_searchlight(B, L, neighborhoods, test_fun)

  # All results should have identical B matrices
  expected_B <- as.matrix(B)
  for (i in seq_along(result)) {
    expect_equal(result[[i]], expected_B, tolerance = 1e-10)
  }
})

test_that("latent_searchlight handles single voxel neighborhoods", {
  set.seed(1111)
  n_time <- 5
  n_vox <- 5
  k <- 3

  B <- make_matrix(n_time, k, seed = 1111)
  L <- make_matrix(n_vox, k, seed = 1112)

  neighborhoods <- list(1L, 2L, 3L, 4L, 5L)

  # Function that returns L_V dimensions
  test_fun <- function(B, L_V, M_V, idx, ...) {
    list(L_V_dim = dim(L_V), M_V_dim = dim(M_V))
  }

  result <- latent_searchlight(B, L, neighborhoods, test_fun)

  expect_length(result, 5)
  for (i in seq_along(result)) {
    # L_V should be 1 x k (single voxel)
    expect_equal(result[[i]]$L_V_dim, c(1, k))
    # M_V should still be k x k
    expect_equal(result[[i]]$M_V_dim, c(k, k))
  }
})

test_that("latent_searchlight function can access all components for reconstruction", {
  set.seed(1212)
  n_time <- 4
  n_vox <- 8
  k <- 2

  B <- make_matrix(n_time, k, seed = 1212)
  L <- make_matrix(n_vox, k, seed = 1213)

  neighborhoods <- list(1:4, 5:8)

  # Function that reconstructs voxel time series using B and L_V
  test_fun <- function(B, L_V, M_V, idx, ...) {
    # Reconstruct: X_V = B %*% t(L_V) gives time x voxels matrix
    as.matrix(B %*% t(L_V))
  }

  result <- latent_searchlight(B, L, neighborhoods, test_fun)

  expect_length(result, 2)

  # Verify reconstruction matches expected
  for (i in seq_along(neighborhoods)) {
    L_sub <- L[neighborhoods[[i]], , drop = FALSE]
    expected_recon <- as.matrix(B %*% t(L_sub))
    expect_equal(result[[i]], expected_recon, tolerance = 1e-10)
  }
})

test_that("latent_searchlight handles neighborhoods with repeated indices", {
  set.seed(1313)
  n_time <- 5
  n_vox <- 10
  k <- 2

  B <- make_matrix(n_time, k, seed = 1313)
  L <- make_matrix(n_vox, k, seed = 1314)

  # Neighborhoods with repeated indices (unusual but should work)
  neighborhoods <- list(
    c(1L, 1L, 2L),   # voxel 1 appears twice
    c(3L, 4L, 5L)
  )

  test_fun <- function(B, L_V, M_V, idx, ...) {
    list(L_V_rows = nrow(L_V), idx = idx)
  }

  result <- latent_searchlight(B, L, neighborhoods, test_fun)

  expect_length(result, 2)
  # L_V should have 3 rows (including duplicate)
  expect_equal(result[[1]]$L_V_rows, 3)
  expect_equal(result[[1]]$idx, c(1L, 1L, 2L))
})

test_that("latent_searchlight works with k=1 (single component)", {
  set.seed(1414)
  n_time <- 6
  n_vox <- 9
  k <- 1

  B <- make_matrix(n_time, k, seed = 1414)
  L <- make_matrix(n_vox, k, seed = 1415)

  neighborhoods <- list(1:3, 4:6, 7:9)

  test_fun <- function(B, L_V, M_V, idx, ...) {
    list(
      B_dim = dim(B),
      L_V_dim = dim(L_V),
      M_V_dim = dim(M_V),
      M_V_val = as.numeric(M_V)
    )
  }

  result <- latent_searchlight(B, L, neighborhoods, test_fun)

  for (i in seq_along(result)) {
    expect_equal(result[[i]]$B_dim, c(n_time, 1))
    expect_equal(result[[i]]$L_V_dim, c(3, 1))
    expect_equal(result[[i]]$M_V_dim, c(1, 1))

    # M_V should be sum of squared loadings for k=1
    L_sub <- L[neighborhoods[[i]], , drop = FALSE]
    expect_equal(result[[i]]$M_V_val, sum(L_sub^2), tolerance = 1e-10)
  }
})

test_that("latent_searchlight function can return various types", {
  set.seed(1515)
  n_time <- 4
  n_vox <- 8
  k <- 2

  B <- make_matrix(n_time, k, seed = 1515)
  L <- make_matrix(n_vox, k, seed = 1516)

  neighborhoods <- list(1:4, 5:8)

  # Function returning different types
  test_fun <- function(B, L_V, M_V, idx, type) {
    switch(type,
      "scalar" = sum(L_V),
      "vector" = rowMeans(as.matrix(L_V)),  # Convert Matrix to matrix for rowMeans
      "matrix" = as.matrix(M_V),
      "list" = list(a = 1, b = 2),
      "character" = paste(idx, collapse = "-")
    )
  }

  # Test scalar return
  result_scalar <- latent_searchlight(B, L, neighborhoods, test_fun, type = "scalar")
  expect_type(result_scalar[[1]], "double")
  expect_length(result_scalar[[1]], 1)

  # Test vector return
  result_vector <- latent_searchlight(B, L, neighborhoods, test_fun, type = "vector")
  expect_length(result_vector[[1]], 4)

  # Test matrix return
  result_matrix <- latent_searchlight(B, L, neighborhoods, test_fun, type = "matrix")
  expect_true(is.matrix(result_matrix[[1]]))

  # Test list return
  result_list <- latent_searchlight(B, L, neighborhoods, test_fun, type = "list")
  expect_type(result_list[[1]], "list")

  # Test character return
  result_char <- latent_searchlight(B, L, neighborhoods, test_fun, type = "character")
  expect_equal(result_char[[1]], "1-2-3-4")
  expect_equal(result_char[[2]], "5-6-7-8")
})

test_that("latent_searchlight with large neighborhoods works correctly", {
  set.seed(1616)
  n_time <- 10
  n_vox <- 100
  k <- 5

  B <- make_matrix(n_time, k, seed = 1616)
  L <- make_matrix(n_vox, k, seed = 1617)

  # Large neighborhoods
  neighborhoods <- list(
    1:50,
    51:100
  )

  test_fun <- function(B, L_V, M_V, idx, ...) {
    list(
      n_vox = nrow(L_V),
      trace_M = sum(Matrix::diag(M_V))  # Use Matrix::diag for Matrix objects
    )
  }

  result <- latent_searchlight(B, L, neighborhoods, test_fun)

  expect_equal(result[[1]]$n_vox, 50)
  expect_equal(result[[2]]$n_vox, 50)

  # Verify trace of Gram matrix
  for (i in 1:2) {
    L_sub <- L[neighborhoods[[i]], , drop = FALSE]
    expected_trace <- sum(Matrix::diag(Matrix::crossprod(L_sub)))
    expect_equal(result[[i]]$trace_M, expected_trace, tolerance = 1e-10)
  }
})

# ------------------------------------------------------------------------------
# Edge case tests
# ------------------------------------------------------------------------------

test_that("compute_local_gram handles single neighborhood", {
  set.seed(1717)
  n_vox <- 10
  k <- 3
  L <- make_matrix(n_vox, k, seed = 1717)

  neighborhoods <- list(1:5)

  result_list <- compute_local_gram(L, neighborhoods, simplify = FALSE)
  result_array <- compute_local_gram(L, neighborhoods, simplify = TRUE)

  expect_length(result_list, 1)
  expect_equal(dim(result_array), c(k, k, 1))
  expect_equal(as.matrix(result_list[[1]]), result_array[, , 1], tolerance = 1e-10)
})

test_that("latent_searchlight handles single neighborhood", {
  set.seed(1818)
  n_time <- 5
  n_vox <- 10
  k <- 2

  B <- make_matrix(n_time, k, seed = 1818)
  L <- make_matrix(n_vox, k, seed = 1819)

  neighborhoods <- list(1:5)

  test_fun <- function(B, L_V, M_V, idx, ...) {
    sum(L_V)
  }

  result <- latent_searchlight(B, L, neighborhoods, test_fun)

  expect_length(result, 1)
  expect_equal(result[[1]], sum(L[1:5, ]), tolerance = 1e-10)
})

test_that("functions work with non-contiguous neighborhood indices", {
  set.seed(1919)
  n_vox <- 20
  k <- 3
  L <- make_matrix(n_vox, k, seed = 1919)

  # Non-contiguous indices
  neighborhoods <- list(
    c(1L, 3L, 5L, 7L, 9L),
    c(2L, 4L, 6L, 8L, 10L)
  )

  result <- compute_local_gram(L, neighborhoods, simplify = FALSE)

  expect_length(result, 2)

  # Verify computation
  for (i in seq_along(neighborhoods)) {
    L_sub <- L[neighborhoods[[i]], , drop = FALSE]
    expected <- Matrix::crossprod(L_sub)
    expect_equal(as.matrix(result[[i]]), as.matrix(expected), tolerance = 1e-10)
  }
})
