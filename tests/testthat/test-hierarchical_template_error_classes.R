# Structured-condition coverage for hierarchical_template.R
# (bd-01KSQQNYTQCW36WNEYPYMX8188). Bare stop() calls were converted to classed
# .encoder_cli_abort() conditions. These triggers fire during validation, before
# any RSpectra/rgsp work, so they are dependency-free.

library(testthat)

small_mask <- function() array(TRUE, dim = c(2, 2, 2))  # 8 voxels

test_that("build_hierarchical_template rejects an empty parcellation list", {
  expect_error(
    build_hierarchical_template(small_mask(), parcellations = list(),
                                k_per_level = integer(0)),
    class = "fmrilatent_error_invalid_argument"
  )
  expect_error(
    build_hierarchical_template(small_mask(), parcellations = list(),
                                k_per_level = integer(0)),
    "non-empty list"
  )
})

test_that("build_hierarchical_template requires k_per_level to match level count", {
  expect_error(
    build_hierarchical_template(small_mask(), parcellations = list(rep(1L, 8)),
                                k_per_level = c(2L, 2L)),
    class = "fmrilatent_error_invalid_argument"
  )
})

test_that("build_hierarchical_template flags a parcellation/mask length mismatch", {
  expect_error(
    build_hierarchical_template(small_mask(), parcellations = list(rep(1L, 4)),
                                k_per_level = 2L),
    class = "fmrilatent_error_dimension_mismatch"
  )
})

test_that("encode_hierarchical and project_hierarchical reject non-templates", {
  expect_error(
    encode_hierarchical(matrix(0, 2, 8), template = 42),
    class = "fmrilatent_error_invalid_type"
  )
  expect_error(
    project_hierarchical(42, matrix(0, 2, 8)),
    class = "fmrilatent_error_invalid_type"
  )
})

test_that("load_template rejects an RDS without a supported template", {
  f <- tempfile(fileext = ".rds")
  on.exit(unlink(f), add = TRUE)
  saveRDS(list(not = "a template"), f)
  expect_error(load_template(f), class = "fmrilatent_error_invalid_type")
})
