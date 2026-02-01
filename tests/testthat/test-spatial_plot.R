# Tests for spatial_plot.R

test_that("plot_spatial_atom works with matrix loadings", {
  set.seed(42)
  mask_arr <- array(TRUE, dim = c(4, 4, 2))
  n_vox <- sum(mask_arr)
  k <- 3

  loadings <- matrix(rnorm(n_vox * k), nrow = n_vox, ncol = k)

  result <- plot_spatial_atom(loadings, mask_arr, idx = 1)

  # Returns the 3D array invisibly
  expect_true(is.array(result))
  expect_equal(dim(result), c(4, 4, 2))
})

test_that("plot_spatial_atom respects idx parameter", {
  set.seed(123)
  mask_arr <- array(TRUE, dim = c(3, 3, 2))
  n_vox <- sum(mask_arr)
  k <- 5

  loadings <- matrix(rnorm(n_vox * k), nrow = n_vox, ncol = k)

  # Different idx values should produce different outputs
  result1 <- plot_spatial_atom(loadings, mask_arr, idx = 1)
  result2 <- plot_spatial_atom(loadings, mask_arr, idx = 2)

  expect_false(identical(result1, result2))
})

test_that("plot_spatial_atom validates idx range", {
  mask_arr <- array(TRUE, dim = c(2, 2, 2))
  n_vox <- sum(mask_arr)
  k <- 3

  loadings <- matrix(rnorm(n_vox * k), nrow = n_vox, ncol = k)

  expect_error(plot_spatial_atom(loadings, mask_arr, idx = 0), "out of range")
  expect_error(plot_spatial_atom(loadings, mask_arr, idx = 4), "out of range")
  expect_error(plot_spatial_atom(loadings, mask_arr, idx = -1), "out of range")
})

test_that("plot_spatial_atom rejects invalid loadings", {
  mask_arr <- array(TRUE, dim = c(2, 2, 2))

  expect_error(plot_spatial_atom(NULL, mask_arr, idx = 1), "must be a matrix")
  expect_error(plot_spatial_atom(c(1, 2, 3), mask_arr, idx = 1), "must be a matrix")
})

test_that("plot_spatial_atom rejects invalid mask", {
  loadings <- matrix(rnorm(8 * 2), nrow = 8, ncol = 2)

  # Use a function which cannot be converted to an array
  expect_error(plot_spatial_atom(loadings, function(x) x, idx = 1), "must be array-like")
})

test_that("plot_spatial_atom works with LogicalNeuroVol mask", {
  skip_if_not_installed("neuroim2")

  set.seed(456)
  mask_arr <- array(TRUE, dim = c(3, 3, 3))
  mask <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(3, 3, 3)))

  n_vox <- sum(mask_arr)
  k <- 2
  loadings <- matrix(rnorm(n_vox * k), nrow = n_vox, ncol = k)

  result <- plot_spatial_atom(loadings, mask, idx = 1)

  expect_true(is.array(result))
  expect_equal(dim(result), c(3, 3, 3))
})

test_that("plot_spatial_atom uses custom title", {
  set.seed(789)
  mask_arr <- array(TRUE, dim = c(2, 2, 2))
  n_vox <- sum(mask_arr)
  k <- 2
  loadings <- matrix(rnorm(n_vox * k), nrow = n_vox, ncol = k)

  # Should not error with custom main
  result <- plot_spatial_atom(loadings, mask_arr, idx = 1, main = "Custom Title")

  expect_true(is.array(result))
})

test_that("plot_spatial_atom correctly maps voxels", {
  # Create a mask with specific pattern
  mask_arr <- array(FALSE, dim = c(3, 3, 1))
  mask_arr[2, 2, 1] <- TRUE
  mask_arr[1, 1, 1] <- TRUE
  mask_arr[3, 3, 1] <- TRUE

  n_vox <- sum(mask_arr)  # 3 voxels
  k <- 1

  # Loadings with known values
  loadings <- matrix(c(10, 20, 30), nrow = 3, ncol = 1)

  result <- plot_spatial_atom(loadings, mask_arr, idx = 1)

  # Check that values are placed correctly
  # The order follows which(mask_arr) which goes column-major
  idx_order <- which(mask_arr)

  # Values should match at mask positions
  expect_equal(result[mask_arr], loadings[, 1])

  # Non-mask positions should be zero
  expect_equal(result[!mask_arr], rep(0, sum(!mask_arr)))
})

test_that("plot_spatial_atom extracts correct component column", {
  mask_arr <- array(TRUE, dim = c(2, 2, 1))  # 4 voxels, all active
  n_vox <- 4
  k <- 3

  # Different values in each component column
  loadings <- matrix(c(1, 2, 3, 4,       # component 1
                       10, 20, 30, 40,   # component 2
                       100, 200, 300, 400), nrow = n_vox, ncol = k)

  # Test idx = 1
  result1 <- plot_spatial_atom(loadings, mask_arr, idx = 1L)
  expect_equal(as.vector(result1), c(1, 2, 3, 4))

  # Test idx = 2
  result2 <- plot_spatial_atom(loadings, mask_arr, idx = 2L)
  expect_equal(as.vector(result2), c(10, 20, 30, 40))

  # Test idx = 3
  result3 <- plot_spatial_atom(loadings, mask_arr, idx = 3L)
  expect_equal(as.vector(result3), c(100, 200, 300, 400))
})

test_that("plot_spatial_atom works when idx equals ncol (boundary case)", {
  mask_arr <- array(TRUE, dim = c(2, 2, 1))
  loadings <- matrix(1:12, nrow = 4, ncol = 3)

  # idx = 3 should work (boundary of ncol)
  result <- plot_spatial_atom(loadings, mask_arr, idx = 3L)
  expect_true(is.array(result))
  expect_equal(as.vector(result), 9:12)
})

test_that("plot_spatial_atom works with single-voxel mask", {
  mask_arr <- array(FALSE, dim = c(2, 2, 2))
  mask_arr[1, 1, 1] <- TRUE

  loadings <- matrix(c(5.5, 10.2), nrow = 1, ncol = 2)

  result <- plot_spatial_atom(loadings, mask_arr, idx = 1L)

  expect_equal(dim(result), c(2, 2, 2))
  expect_equal(result[1, 1, 1], 5.5)
  expect_equal(sum(result), 5.5)  # Only one non-zero value
})

test_that("plot_spatial_atom handles negative loading values", {
  mask_arr <- array(TRUE, dim = c(2, 2, 1))
  loadings <- matrix(c(-1.5, 2.0, -3.5, 4.0), nrow = 4, ncol = 1)

  result <- plot_spatial_atom(loadings, mask_arr, idx = 1L)

  expect_equal(as.vector(result), c(-1.5, 2.0, -3.5, 4.0))
})

test_that("plot_spatial_atom handles sparse Matrix loadings", {
  skip_if_not_installed("Matrix")

  mask_arr <- array(TRUE, dim = c(2, 2, 1))
  n_vox <- 4
  k <- 2

  # Create a sparse Matrix
  loadings <- Matrix::Matrix(
    c(1, 0, 0, 2,
      0, 3, 4, 0),
    nrow = n_vox, ncol = k, sparse = TRUE
  )

  result <- plot_spatial_atom(loadings, mask_arr, idx = 1L)
  expect_true(is.array(result))
  expect_equal(dim(result), dim(mask_arr))
  expect_equal(as.vector(result), c(1, 0, 0, 2))
})

test_that("plot_spatial_atom handles 3D mask with multiple z slices", {
  # Create a 3x3x3 mask
  mask_arr <- array(FALSE, dim = c(3, 3, 3))
  # Set specific voxels to TRUE
  mask_arr[1, 1, 1] <- TRUE
  mask_arr[2, 2, 2] <- TRUE
  mask_arr[3, 3, 3] <- TRUE
  n_vox <- sum(mask_arr)  # 3 voxels

  loadings <- matrix(c(1.0, 2.0, 3.0), nrow = n_vox, ncol = 1)

  result <- plot_spatial_atom(loadings, mask_arr, idx = 1L)

  expect_equal(dim(result), c(3, 3, 3))
  expect_equal(result[1, 1, 1], 1.0)
  expect_equal(result[2, 2, 2], 2.0)
  expect_equal(result[3, 3, 3], 3.0)

  # Check that non-mask voxels are 0
  expect_equal(result[1, 2, 1], 0)
  expect_equal(result[2, 1, 3], 0)
})

test_that("plot_spatial_atom does not error with ggplot2 available", {
  skip_if_not_installed("ggplot2")

  mask_arr <- array(TRUE, dim = c(3, 3, 2))
  loadings <- matrix(rnorm(18 * 2), nrow = 18, ncol = 2)

  # The function returns the array invisibly regardless of ggplot availability
  expect_no_error({
    result <- plot_spatial_atom(loadings, mask_arr, idx = 1L, main = "Test Plot")
  })
})

test_that("plot_spatial_atom works with LoadingsHandle", {
  skip_if_not_installed("Matrix")

  # Create explicit LoadingsHandle
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
    id = "test-spatial-plot-loadings",
    dim = as.integer(dim(load_mat)),
    kind = "explicit",
    spec = list(matrix = load_mat),
    label = "test-loadings"
  )

  mask_arr <- array(TRUE, dim = c(2, 2, 1))

  result <- plot_spatial_atom(l_handle, mask_arr, idx = 1L)

  expect_true(is.array(result))
  expect_equal(dim(result), dim(mask_arr))
  # Values should match first column of load_mat
  expect_equal(as.vector(result), c(1.0, 0.5, 0.1, 0.0))
})

test_that("plot_spatial_atom with LoadingsHandle extracts correct component", {
  skip_if_not_installed("Matrix")

  # Create explicit LoadingsHandle with multiple columns
  load_mat <- Matrix::Matrix(
    rbind(
      c(1.0, 10.0, 100.0),
      c(2.0, 20.0, 200.0),
      c(3.0, 30.0, 300.0),
      c(4.0, 40.0, 400.0)
    ),
    sparse = FALSE
  )

  l_handle <- new("LoadingsHandle",
    id = "test-spatial-plot-loadings-multi",
    dim = as.integer(dim(load_mat)),
    kind = "explicit",
    spec = list(matrix = load_mat),
    label = "test-loadings-multi"
  )

  mask_arr <- array(TRUE, dim = c(2, 2, 1))

  # Test idx = 2
  result <- plot_spatial_atom(l_handle, mask_arr, idx = 2L)
  expect_equal(as.vector(result), c(10.0, 20.0, 30.0, 40.0))

  # Test idx = 3
  result3 <- plot_spatial_atom(l_handle, mask_arr, idx = 3L)
  expect_equal(as.vector(result3), c(100.0, 200.0, 300.0, 400.0))
})

test_that("plot_spatial_atom preserves zero values in loadings", {
  mask_arr <- array(TRUE, dim = c(2, 2, 1))
  # Loadings with explicit zeros
  loadings <- matrix(c(0, 1, 0, 2), nrow = 4, ncol = 1)

  result <- plot_spatial_atom(loadings, mask_arr, idx = 1L)

  expect_equal(result[1, 1, 1], 0)
  expect_equal(result[2, 1, 1], 1)
  expect_equal(result[1, 2, 1], 0)
  expect_equal(result[2, 2, 1], 2)
})

test_that("plot_spatial_atom works with partial mask", {
  # Not all voxels in the array are part of the mask
  mask_arr <- array(FALSE, dim = c(3, 3, 1))
  mask_arr[1, 1, 1] <- TRUE
  mask_arr[3, 3, 1] <- TRUE  # Only 2 voxels masked

  loadings <- matrix(c(5, 10), nrow = 2, ncol = 1)

  result <- plot_spatial_atom(loadings, mask_arr, idx = 1L)

  expect_equal(dim(result), c(3, 3, 1))
  expect_equal(result[1, 1, 1], 5)
  expect_equal(result[3, 3, 1], 10)
  # All other positions should be 0
  expect_equal(sum(result == 0), 7)
})

test_that("plot_spatial_atom works with default main parameter", {
  mask_arr <- array(TRUE, dim = c(2, 2, 1))
  loadings <- matrix(1:4, nrow = 4, ncol = 1)

  # Should not error when main is not provided (default NULL)
  result <- plot_spatial_atom(loadings, mask_arr, idx = 1L)
  expect_true(is.array(result))
})
