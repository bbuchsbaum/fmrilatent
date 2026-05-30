# Structured-condition coverage for haar_wavelet.R
# (bd-01KSQQNYTQCW36WNEYPYMX8188). Bare stop() calls were converted to classed
# .encoder_cli_abort() conditions. Triggers below are dependency-free.

library(testthat)

test_that("haar_wavelet_forward rejects a non-3D mask", {
  expect_error(
    haar_wavelet_forward(matrix(0, nrow = 3, ncol = 4), array(TRUE, dim = c(2, 2))),
    class = "fmrilatent_error_invalid_mask"
  )
})

test_that("precompute_haar_scalings rejects a non-positive levels argument", {
  expect_error(
    fmrilatent:::precompute_haar_scalings(array(TRUE, dim = c(2, 2, 2)), levels = 0L),
    class = "fmrilatent_error_invalid_argument"
  )
})

test_that("lna_forward_lift_matrix rejects a non-matrix data argument", {
  expect_error(
    lna_forward_lift_matrix("not a matrix", logical(8), c(2L, 2L, 2L), 1L, list()),
    class = "fmrilatent_error_invalid_type"
  )
})

test_that("get_roi_detail_indices flags an ROI/mask dimension mismatch", {
  expect_error(
    fmrilatent:::get_roi_detail_indices(
      array(TRUE, dim = c(2, 2, 2)),
      array(TRUE, dim = c(4, 4, 4)),
      levels = 1L
    ),
    class = "fmrilatent_error_dimension_mismatch"
  )
})

test_that("as_haar_latent rejects a non-Haar object (class + message preserved)", {
  expect_error(as_haar_latent(list()), class = "fmrilatent_error_invalid_type")
  # The original message text survives the cli_abort conversion.
  expect_error(as_haar_latent(list()), "not Haar implicit latent")
})
