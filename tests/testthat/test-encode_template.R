library(testthat)

test_that("symmetric matrix factor fallback preserves the source matrix", {
  mat <- matrix(c(1, 1, 1, 1), nrow = 2)

  factor <- .symmetric_matrix_factor(mat)

  expect_equal(factor, t(factor), tolerance = 1e-12)
  expect_equal(t(factor) %*% factor, mat, tolerance = 1e-12)
})

test_that("symmetric matrix factor rejects non-square matrices before symmetrizing", {
  expect_error(
    .symmetric_matrix_factor(matrix(seq_len(6), nrow = 2)),
    class = "fmrilatent_error_dim"
  )
})

test_that("metric transforms use pseudoinverse fallback for rank-deficient factors", {
  raw_loadings <- matrix(c(1, 1, 1, 1), nrow = 2)
  payload <- .template_coordinate_payload(
    raw_loadings = raw_loadings,
    default_measure = "unit"
  )
  transform <- payload$analysis_transform
  z_analysis <- matrix(c(2, 2), ncol = 1)

  expect_equal(transform$matrix %*% transform$to_raw(z_analysis), z_analysis, tolerance = 1e-12)
  expect_equal(
    crossprod(payload$analysis_loadings),
    transform$matrix %*% transform$inverse_matrix,
    tolerance = 1e-12
  )
  expect_equal(
    .transform_quadratic_form(diag(2), transform),
    t(transform$inverse_matrix) %*% transform$inverse_matrix,
    tolerance = 1e-12
  )
})
