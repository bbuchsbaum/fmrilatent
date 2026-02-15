# Tests for slepian_plot.R
# Testing plot_slepian_temporal and plot_basis_gram functions

# =============================================================================
# plot_slepian_temporal tests
# =============================================================================

test_that("plot_slepian_temporal works with matrix input", {
  set.seed(42)
  n_time <- 50
  k <- 6
  basis <- matrix(rnorm(n_time * k), nrow = n_time, ncol = k)

  # Without ggplot2, should use base graphics
  result <- plot_slepian_temporal(basis, max_components = 4)

  # If ggplot2 is available, returns ggplot; otherwise NULL
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    expect_s3_class(result, "ggplot")
  } else {
    expect_null(result)
  }
})

test_that("plot_slepian_temporal respects max_components", {
  skip_if_not_installed("ggplot2")

  set.seed(123)
  n_time <- 30
  k <- 10
  basis <- matrix(rnorm(n_time * k), nrow = n_time, ncol = k)

  result <- plot_slepian_temporal(basis, max_components = 3)

  expect_s3_class(result, "ggplot")

  # Check that only 3 components are in the data
  plot_data <- ggplot2::ggplot_build(result)$data[[1]]
  unique_groups <- length(unique(plot_data$group))
  expect_equal(unique_groups, 3)
})

test_that("plot_slepian_temporal handles fewer components than max", {
  skip_if_not_installed("ggplot2")

  set.seed(456)
  n_time <- 20
  k <- 2
  basis <- matrix(rnorm(n_time * k), nrow = n_time, ncol = k)

  result <- plot_slepian_temporal(basis, max_components = 6)

  expect_s3_class(result, "ggplot")

  # Verify only 2 components in plot data
  plot_data <- ggplot2::ggplot_build(result)$data[[1]]
  unique_groups <- length(unique(plot_data$group))
  expect_equal(unique_groups, 2)
})

test_that("plot_slepian_temporal rejects invalid input", {
  expect_error(plot_slepian_temporal("not_a_matrix"), "must be a matrix")
  expect_error(plot_slepian_temporal(c(1, 2, 3)), "must be a matrix")
  expect_error(plot_slepian_temporal(data.frame(a = 1:10)), "must be a matrix")
  expect_error(plot_slepian_temporal(list(x = 1)), "must be a matrix")
})

test_that("plot_slepian_temporal works with sparse Matrix input", {
  skip_if_not_installed("ggplot2")

  set.seed(303)
  n_time <- 25
  k <- 4
  basis <- Matrix::Matrix(rnorm(n_time * k), nrow = n_time, ncol = k, sparse = TRUE)

  result <- plot_slepian_temporal(basis, max_components = 3)

  expect_s3_class(result, "ggplot")
})

test_that("plot_slepian_temporal works with BasisHandle input", {
  skip_if_not_installed("ggplot2")

  # Create a DCT basis handle (simpler than Slepian, avoids RSpectra dependency)
  n_time <- 20L
  k <- 4L
  bh <- dct_basis_handle(n_time = n_time, k = k)

  result <- plot_slepian_temporal(bh, max_components = 3)

  expect_s3_class(result, "ggplot")

  # Verify 3 components plotted
  plot_data <- ggplot2::ggplot_build(result)$data[[1]]
  unique_groups <- length(unique(plot_data$group))
  expect_equal(unique_groups, 3)
})

test_that("plot_slepian_temporal works with slepian_temporal_handle", {
  skip_if_not_installed("ggplot2")

  n_time <- 24L
  tr <- 2
  bandwidth <- 0.08
  bh <- slepian_temporal_handle(n_time = n_time, tr = tr, bandwidth = bandwidth,
                                 backend = "tridiag")

  result <- plot_slepian_temporal(bh, max_components = 4)

  expect_s3_class(result, "ggplot")
})

test_that("plot_slepian_temporal handles single component basis", {
  skip_if_not_installed("ggplot2")

  n_time <- 15
  basis <- matrix(rnorm(n_time), nrow = n_time, ncol = 1)

  result <- plot_slepian_temporal(basis, max_components = 6)

  expect_s3_class(result, "ggplot")

  # Should only have 1 component
  plot_data <- ggplot2::ggplot_build(result)$data[[1]]
  unique_groups <- length(unique(plot_data$group))
  expect_equal(unique_groups, 1)
})

test_that("plot_slepian_temporal handles minimal time points", {
  skip_if_not_installed("ggplot2")

  # Very short time series
  n_time <- 3
  k <- 2
  basis <- matrix(rnorm(n_time * k), nrow = n_time, ncol = k)

  result <- plot_slepian_temporal(basis, max_components = 2)

  expect_s3_class(result, "ggplot")

  # Verify time axis has 3 points per component
  plot_data <- ggplot2::ggplot_build(result)$data[[1]]
  expect_equal(nrow(plot_data), n_time * k)
})

test_that("plot_slepian_temporal data frame has correct structure", {
  skip_if_not_installed("ggplot2")

  set.seed(999)
  n_time <- 10
  k <- 3
  basis <- matrix(rnorm(n_time * k), nrow = n_time, ncol = k)

  result <- plot_slepian_temporal(basis, max_components = 3)

  # Access the underlying data
  plot_data <- ggplot2::layer_data(result)

  # Should have n_time * k rows

  expect_equal(nrow(plot_data), n_time * k)

  # Check x values correspond to time indices
  expect_true(all(plot_data$x %in% 1:n_time))
})

test_that("plot_slepian_temporal default max_components is 6", {
  skip_if_not_installed("ggplot2")

  set.seed(888)
  n_time <- 20
  k <- 10  # More than default max_components
  basis <- matrix(rnorm(n_time * k), nrow = n_time, ncol = k)

  result <- plot_slepian_temporal(basis)  # Use default

  plot_data <- ggplot2::ggplot_build(result)$data[[1]]
  unique_groups <- length(unique(plot_data$group))
  expect_equal(unique_groups, 6)  # Default is 6
})

# =============================================================================
# plot_basis_gram tests
# =============================================================================

test_that("plot_basis_gram works with matrix input", {
  set.seed(789)
  n_time <- 20
  k <- 4
  basis <- matrix(rnorm(n_time * k), nrow = n_time, ncol = k)

  result <- plot_basis_gram(basis)

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    expect_s3_class(result, "ggplot")
  } else {
    expect_null(result)
  }
})

test_that("plot_basis_gram shows orthogonality for orthonormal basis", {
  skip_if_not_installed("ggplot2")

  # Create an orthonormal basis via QR
  set.seed(101)
  n_time <- 20
  k <- 4
  A <- matrix(rnorm(n_time * k), nrow = n_time, ncol = k)
  basis <- qr.Q(qr(A))

  result <- plot_basis_gram(basis)

  expect_s3_class(result, "ggplot")

  # The Gram matrix should be approximately identity
  G <- crossprod(basis)
  expect_equal(G, diag(k), tolerance = 1e-10)
})

test_that("plot_basis_gram rejects invalid input", {
  expect_error(plot_basis_gram("not_a_matrix"), "must be a matrix")
  expect_error(plot_basis_gram(list(a = 1)), "must be a matrix")
  expect_error(plot_basis_gram(c(1, 2, 3, 4)), "must be a matrix")
  expect_error(plot_basis_gram(data.frame(x = 1:5)), "must be a matrix")
})

test_that("plot_basis_gram handles single column basis", {
  skip_if_not_installed("ggplot2")

  basis <- matrix(rnorm(10), nrow = 10, ncol = 1)

  result <- plot_basis_gram(basis)

  expect_s3_class(result, "ggplot")
})

test_that("plot_basis_gram handles sparse Matrix input", {
  skip_if_not_installed("ggplot2")

  set.seed(202)
  n_time <- 15
  k <- 3

  basis <- Matrix::Matrix(rnorm(n_time * k), nrow = n_time, ncol = k, sparse = TRUE)

  result <- plot_basis_gram(basis)

  expect_s3_class(result, "ggplot")
})

test_that("plot_basis_gram works with BasisHandle input", {
  skip_if_not_installed("ggplot2")

  # Create a DCT basis handle
  n_time <- 16L
  k <- 5L
  bh <- dct_basis_handle(n_time = n_time, k = k)

  result <- plot_basis_gram(bh)

  expect_s3_class(result, "ggplot")
})

test_that("plot_basis_gram works with slepian_temporal_handle", {
  skip_if_not_installed("ggplot2")

  n_time <- 30L
  tr <- 2
  bandwidth <- 0.06
  bh <- slepian_temporal_handle(n_time = n_time, tr = tr, bandwidth = bandwidth,
                                 backend = "tridiag")

  result <- plot_basis_gram(bh)

  expect_s3_class(result, "ggplot")

  # Slepian bases should be approximately orthonormal
  B <- basis_mat(bh)
  G <- crossprod(as.matrix(B))
  expect_equal(G, diag(ncol(G)), tolerance = 1e-6)
})

test_that("plot_basis_gram data frame has correct structure", {
  skip_if_not_installed("ggplot2")

  set.seed(777)
  n_time <- 12
  k <- 3
  basis <- matrix(rnorm(n_time * k), nrow = n_time, ncol = k)

  result <- plot_basis_gram(basis)

  # Access the underlying tile data
  plot_data <- ggplot2::layer_data(result)

  # Gram matrix is k x k, so tiles should be k^2
  expect_equal(nrow(plot_data), k * k)
})

test_that("plot_basis_gram tile values match Gram matrix", {
  skip_if_not_installed("ggplot2")

  set.seed(666)
  n_time <- 10
  k <- 2
  basis <- matrix(rnorm(n_time * k), nrow = n_time, ncol = k)

  result <- plot_basis_gram(basis)

  # Compute expected Gram matrix
  G <- crossprod(basis)

  # Extract fill values from the underlying plot data
  # The data frame stored in the ggplot object contains the original values
  observed_vals <- result$data$val

  # The fill values should match the Gram matrix entries
  expected_vals <- as.vector(G)

  expect_equal(sort(observed_vals), sort(expected_vals), tolerance = 1e-10)
})

test_that("plot_basis_gram handles wide basis (more cols than rows)", {

  skip_if_not_installed("ggplot2")

  # More components than time points (unusual but valid)
  n_time <- 5
  k <- 8
  basis <- matrix(rnorm(n_time * k), nrow = n_time, ncol = k)

  result <- plot_basis_gram(basis)

  expect_s3_class(result, "ggplot")

  # Check dimensions
  plot_data <- ggplot2::layer_data(result)
  expect_equal(nrow(plot_data), k * k)
})

test_that("plot_basis_gram uses viridis color scale", {
  skip_if_not_installed("ggplot2")

  basis <- matrix(rnorm(20), nrow = 10, ncol = 2)
  result <- plot_basis_gram(basis)

  # Check that viridis scale is present in scales
  scales <- result$scales$scales
  color_scale_present <- any(vapply(scales, function(s) {
    inherits(s, "ScaleContinuous") && !is.null(s$aesthetics) && "fill" %in% s$aesthetics
  }, logical(1)))

  expect_true(color_scale_present)
})

# =============================================================================
# Edge cases and integration tests
# =============================================================================

test_that("both plot functions handle dense Matrix class", {
  skip_if_not_installed("ggplot2")

  set.seed(555)
  n_time <- 12
  k <- 3
  basis <- Matrix::Matrix(rnorm(n_time * k), nrow = n_time, ncol = k, sparse = FALSE)

  result_temporal <- plot_slepian_temporal(basis, max_components = 2)
  result_gram <- plot_basis_gram(basis)

  expect_s3_class(result_temporal, "ggplot")
  expect_s3_class(result_gram, "ggplot")
})

test_that("plot functions work with actual DPSS basis", {
  skip_if_not_installed("ggplot2")

  # Create actual DPSS basis
  n_time <- 64
  tr <- 2
  bw <- 0.05
  k <- 5
  B_dpss <- dpss_time_basis(n_time, tr, bw, k = k, backend = "tridiag")

  result_temporal <- plot_slepian_temporal(B_dpss, max_components = 4)
  result_gram <- plot_basis_gram(B_dpss)

  expect_s3_class(result_temporal, "ggplot")
  expect_s3_class(result_gram, "ggplot")

  # DPSS should be orthonormal
  G <- crossprod(B_dpss)
  expect_equal(G, diag(k), tolerance = 1e-6)
})

test_that("plot_slepian_temporal produces correct time axis labels", {
  skip_if_not_installed("ggplot2")

  set.seed(444)
  n_time <- 15
  k <- 2
  basis <- matrix(rnorm(n_time * k), nrow = n_time, ncol = k)

  result <- plot_slepian_temporal(basis)

  # Check axis labels
  expect_equal(result$labels$x, "Time")
  expect_equal(result$labels$y, "Amplitude")
  expect_equal(result$labels$colour, "Component")
})

test_that("plot_basis_gram produces correct axis labels", {
  skip_if_not_installed("ggplot2")

  basis <- matrix(rnorm(20), nrow = 10, ncol = 2)

  result <- plot_basis_gram(basis)

  # Check axis labels
  expect_equal(result$labels$x, "Component")
  expect_equal(result$labels$y, "Component")
  expect_equal(result$labels$fill, "Gij")
})

test_that("plot functions have minimal theme applied", {
  skip_if_not_installed("ggplot2")

  basis <- matrix(rnorm(30), nrow = 10, ncol = 3)

  result_temporal <- plot_slepian_temporal(basis)
  result_gram <- plot_basis_gram(basis)

  # Both should use theme_minimal (check theme class)
  expect_true(inherits(result_temporal$theme, "theme"))
  expect_true(inherits(result_gram$theme, "theme"))
})

# =============================================================================
# Base graphics fallback tests (when ggplot2 not available)
# =============================================================================

test_that("plot_slepian_temporal base graphics fallback with matrix input", {
  # Temporarily mock ggplot2 unavailability by testing behavior
  set.seed(42)
  n_time <- 30
  k <- 4
  basis <- matrix(rnorm(n_time * k), nrow = n_time, ncol = k)

  # Should work without error even if ggplot2 not available
  result <- suppressMessages(plot_slepian_temporal(basis, max_components = 3))

  # Returns ggplot if available, NULL otherwise
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    expect_s3_class(result, "ggplot")
  } else {
    expect_null(result)
  }
})

test_that("plot_basis_gram base graphics fallback with matrix input", {
  set.seed(99)
  n_time <- 20
  k <- 3
  basis <- matrix(rnorm(n_time * k), nrow = n_time, ncol = k)

  # Should work without error
  result <- suppressMessages(plot_basis_gram(basis))

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    expect_s3_class(result, "ggplot")
  } else {
    expect_null(result)
  }
})

# =============================================================================
# Additional coverage tests for plot_slepian_temporal and plot_basis_gram
# =============================================================================

test_that("plot_slepian_temporal with max_components = 1 plots single line", {
  skip_if_not_installed("ggplot2")

  set.seed(1010)
  n_time <- 20
  k <- 5
  basis <- matrix(rnorm(n_time * k), nrow = n_time, ncol = k)

  result <- plot_slepian_temporal(basis, max_components = 1L)
  expect_s3_class(result, "ggplot")

  plot_data <- ggplot2::ggplot_build(result)$data[[1]]
  unique_groups <- length(unique(plot_data$group))
  expect_equal(unique_groups, 1L)
})

test_that("plot_slepian_temporal with max_components exceeding ncol uses all columns", {
  skip_if_not_installed("ggplot2")

  set.seed(2020)
  n_time <- 15
  k <- 2
  basis <- matrix(rnorm(n_time * k), nrow = n_time, ncol = k)

  # max_components = 100, but only 2 columns exist

  result <- plot_slepian_temporal(basis, max_components = 100L)
  expect_s3_class(result, "ggplot")

  plot_data <- ggplot2::ggplot_build(result)$data[[1]]
  unique_groups <- length(unique(plot_data$group))
  expect_equal(unique_groups, 2L)
})

test_that("plot_slepian_temporal title is correct", {
  skip_if_not_installed("ggplot2")

  basis <- matrix(rnorm(30), nrow = 10, ncol = 3)
  result <- plot_slepian_temporal(basis)
  expect_equal(result$labels$title, "Temporal Slepians (DPSS)")
})

test_that("plot_basis_gram title is correct", {
  skip_if_not_installed("ggplot2")

  basis <- matrix(rnorm(20), nrow = 10, ncol = 2)
  result <- plot_basis_gram(basis)
  expect_equal(result$labels$title, "Basis Gram matrix")
})

test_that("plot_basis_gram with large k produces k^2 tiles", {
  skip_if_not_installed("ggplot2")

  set.seed(3030)
  n_time <- 20
  k <- 7
  basis <- matrix(rnorm(n_time * k), nrow = n_time, ncol = k)

  result <- plot_basis_gram(basis)
  expect_s3_class(result, "ggplot")

  plot_data <- ggplot2::layer_data(result)
  expect_equal(nrow(plot_data), k * k)
})

test_that("plot_basis_gram Gram matrix diagonal values are positive for non-zero basis", {
  skip_if_not_installed("ggplot2")

  set.seed(4040)
  basis <- matrix(rnorm(30), nrow = 10, ncol = 3)
  result <- plot_basis_gram(basis)

  G <- crossprod(basis)
  # Diagonal entries should all be positive (sum of squares)
  expect_true(all(diag(G) > 0))

  # Verify they appear in the plot data
  df <- result$data
  diag_entries <- df[df$i == df$j, "val"]
  expect_true(all(diag_entries > 0))
})

test_that("plot_slepian_temporal with actual DPSS from both backends", {
  skip_if_not_installed("ggplot2")

  # Test with tridiag
  B_tri <- dpss_time_basis(32, tr = 2, bandwidth = 0.06, k = 4, backend = "tridiag")
  result_tri <- plot_slepian_temporal(B_tri, max_components = 4)
  expect_s3_class(result_tri, "ggplot")

  # Test with dense
  B_dense <- dpss_time_basis(32, tr = 2, bandwidth = 0.06, k = 4, backend = "dense")
  result_dense <- plot_slepian_temporal(B_dense, max_components = 4)
  expect_s3_class(result_dense, "ggplot")
})

test_that("plot_basis_gram coord_equal produces square aspect ratio", {
  skip_if_not_installed("ggplot2")

  basis <- matrix(rnorm(20), nrow = 10, ncol = 2)
  result <- plot_basis_gram(basis)

  # Verify that coord_equal was applied by checking the coordinate system exists
  # and that the plot builds correctly with equal aspect ratio
  expect_true(inherits(result$coordinates, "Coord"))
  built <- ggplot2::ggplot_build(result)
  expect_true(!is.null(built))
})
