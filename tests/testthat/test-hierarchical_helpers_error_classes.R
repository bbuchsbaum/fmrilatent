# Structured-condition coverage for hierarchical_helpers.R
# (bd-01KSQQNYTQCW36WNEYPYMX8188). Bare stop()/warning() calls were converted to
# classed .encoder_cli_abort()/.encoder_cli_warn() conditions. Triggers are
# dependency-free.

library(testthat)

test_that("spectral_ward_hclust guards its W argument with classed conditions", {
  expect_error(spectral_ward_hclust("not_matrix"), class = "fmrilatent_error_invalid_type")
  expect_error(spectral_ward_hclust(matrix(1, nrow = 2, ncol = 3)),
               class = "fmrilatent_error_dimension_mismatch")
})

test_that("validate_nested_parcellations rejects an empty level list", {
  expect_error(validate_nested_parcellations(list()),
               class = "fmrilatent_error_invalid_argument")
  expect_error(validate_nested_parcellations(list()), "non-empty list")
})

test_that("validate_nested_parcellations flags unequal level lengths", {
  expect_error(
    validate_nested_parcellations(list(c(1L, 1L), c(1L, 2L, 3L))),
    class = "fmrilatent_error_dimension_mismatch"
  )
})

test_that("validate_nested_parcellations flags non-nested levels", {
  # child id 1 (positions 1,3) maps to two different parents (1 and 2)
  expect_error(
    validate_nested_parcellations(list(c(1L, 1L, 2L, 2L), c(1L, 2L, 1L, 2L))),
    class = "fmrilatent_error_not_nested"
  )
})

test_that("cut_hclust_nested rejects a non-hclust object", {
  expect_error(cut_hclust_nested(42, c(2L, 4L)),
               class = "fmrilatent_error_invalid_type")
})

test_that("cut_hclust_nested warns (classed) when k_levels are unsorted", {
  set.seed(1)
  hc <- stats::hclust(stats::dist(matrix(rnorm(20), nrow = 10)))
  expect_warning(
    cut_hclust_nested(hc, c(4L, 2L)),
    class = "fmrilatent_warning_unsorted_levels"
  )
})
