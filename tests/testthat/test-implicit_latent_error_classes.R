# Structured-condition coverage for implicit_latent.R
# (bd-01KSQQNYTQCW36WNEYPYMX8188). Bare stop() calls were converted to classed
# .encoder_cli_abort() conditions. Triggers below are dependency-free.

library(testthat)

good_decoder <- function(time_idx = NULL, roi_mask = NULL, levels_keep = NULL, ...) NULL

test_that("implicit_latent requires meta$family", {
  expect_error(
    implicit_latent(coeff = 1, decoder = good_decoder, meta = list()),
    class = "fmrilatent_error_invalid_argument"
  )
})

test_that("implicit_latent rejects a non-function decoder", {
  expect_error(
    implicit_latent(coeff = 1, decoder = 42, meta = list(family = "f"),
                    mask = array(TRUE, dim = c(2, 2, 1))),
    class = "fmrilatent_error_invalid_argument"
  )
})

test_that("implicit_latent rejects a decoder missing required formals", {
  expect_error(
    implicit_latent(coeff = 1, decoder = function(x) x, meta = list(family = "f"),
                    mask = array(TRUE, dim = c(2, 2, 1))),
    class = "fmrilatent_error_invalid_argument"
  )
})

test_that("implicit_latent requires a mask or support", {
  expect_error(
    implicit_latent(coeff = 1, decoder = good_decoder, meta = list(family = "f")),
    class = "fmrilatent_error_missing_argument"
  )
})

test_that("as_implicit_latent has no default coercion method", {
  expect_error(
    as_implicit_latent(42),
    class = "fmrilatent_error_unsupported_operation"
  )
})

test_that("wrap_decoded_volume flags a support-cardinality mismatch", {
  mask_arr <- array(c(TRUE, TRUE, FALSE, FALSE), dim = c(2, 2, 1))  # cardinality 2
  expect_error(
    fmrilatent:::.wrap_decoded_volume(c(1, 2, 3), mask_arr),
    class = "fmrilatent_error_dimension_mismatch"
  )
})

test_that("mask() is unsupported for a non-volumetric implicit latent", {
  lat <- implicit_latent(coeff = 1, decoder = good_decoder, meta = list(family = "f"),
                         support = 1:3)
  expect_error(mask(lat), class = "fmrilatent_error_unsupported_operation")
})
