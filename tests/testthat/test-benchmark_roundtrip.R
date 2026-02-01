# Tests for benchmark_roundtrip.R

test_that("benchmark_roundtrip requires bench package", {
  skip_if_not_installed("bench")
  # If we get here, bench is installed, so just verify we can call the function
  # with minimal settings
  result <- benchmark_roundtrip(
    mask_dims = c(4, 4, 2),
    n_time = 3L,
    methods = "slepian_space",
    iterations = 1L
  )
  expect_s3_class(result, "data.frame")
})

test_that("benchmark_roundtrip returns correct structure for slepian_space", {
  skip_if_not_installed("bench")

  result <- benchmark_roundtrip(
    mask_dims = c(4, 4, 2),
    n_time = 3L,
    methods = "slepian_space",
    iterations = 1L
  )

  # Check structure

  expect_s3_class(result, "data.frame")
  expect_true(all(c("method", "median_ms", "itr", "rmse") %in% names(result)))
  expect_equal(nrow(result), 1L)
  expect_equal(result$method, "slepian_space")
  expect_equal(result$itr, 1L)

  # Check values are sensible

  expect_true(is.numeric(result$median_ms))
  expect_true(result$median_ms >= 0)
  expect_true(is.numeric(result$rmse))
  expect_true(result$rmse >= 0)
})

test_that("benchmark_roundtrip works with hrbf method", {
  skip_if_not_installed("bench")

  result <- benchmark_roundtrip(
    mask_dims = c(4, 4, 2),
    n_time = 3L,
    methods = "hrbf",
    iterations = 1L
  )

  expect_s3_class(result, "data.frame")
  expect_equal(result$method, "hrbf")
  expect_true(is.numeric(result$rmse))
})

test_that("benchmark_roundtrip works with wavelet_active method", {
  skip_if_not_installed("bench")

  result <- benchmark_roundtrip(
    mask_dims = c(4, 4, 2),
    n_time = 4L,
    methods = "wavelet_active",
    iterations = 1L
  )


  expect_s3_class(result, "data.frame")
  expect_equal(result$method, "wavelet_active")
  expect_true(is.numeric(result$rmse))
})

test_that("benchmark_roundtrip works with bspline_hrbf_st method", {
  skip_if_not_installed("bench")

  result <- benchmark_roundtrip(
    mask_dims = c(4, 4, 2),
    n_time = 5L,
    methods = "bspline_hrbf_st",
    iterations = 1L
  )

  expect_s3_class(result, "data.frame")
  expect_equal(result$method, "bspline_hrbf_st")
  expect_true(is.numeric(result$rmse))
})

test_that("benchmark_roundtrip handles multiple methods", {
  skip_if_not_installed("bench")

  result <- benchmark_roundtrip(
    mask_dims = c(4, 4, 2),
    n_time = 5L,
    methods = c("slepian_space", "hrbf"),
    iterations = 1L
  )


  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2L)
  expect_equal(result$method, c("slepian_space", "hrbf"))
})

test_that("benchmark_roundtrip errors on unknown method", {
  skip_if_not_installed("bench")

  expect_error(
    benchmark_roundtrip(
      mask_dims = c(4, 4, 2),
      n_time = 3L,
      methods = "unknown_method",
      iterations = 1L
    ),
    "Unknown method"
  )
})

test_that("plot_benchmark_roundtrip returns ggplot when ggplot2 available", {
  skip_if_not_installed("bench")
  skip_if_not_installed("ggplot2")

  result <- benchmark_roundtrip(
    mask_dims = c(4, 4, 2),
    n_time = 3L,
    methods = "slepian_space",
    iterations = 1L
  )

  p <- plot_benchmark_roundtrip(result)
  expect_s3_class(p, "ggplot")
})

test_that("plot_benchmark_roundtrip works without ggplot2", {
  skip_if_not_installed("bench")

  # Create a mock data frame (no need to run actual benchmark)
  df <- data.frame(
    method = "slepian_space",
    median_ms = 100.0,
    itr = 1L,
    rmse = 0.5,
    stringsAsFactors = FALSE
  )

  # If ggplot2 is available, we can't really test "without" it in this session,

  # but we can verify the function handles the df correctly
  result <- plot_benchmark_roundtrip(df)

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    expect_s3_class(result, "ggplot")
  } else {
    # When ggplot2 is not available, it returns the df invisibly
    expect_equal(result, df)
  }
})

test_that("plot_benchmark_roundtrip handles multiple methods", {
  skip_if_not_installed("bench")
  skip_if_not_installed("ggplot2")

  result <- benchmark_roundtrip(
    mask_dims = c(4, 4, 2),
    n_time = 5L,
    methods = c("slepian_space", "hrbf"),
    iterations = 1L
  )

  p <- plot_benchmark_roundtrip(result)
  expect_s3_class(p, "ggplot")

  # Verify the plot has the expected aesthetics
  expect_true("colour" %in% names(p$mapping) || "color" %in% names(p$mapping))
})
