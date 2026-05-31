library(testthat)

test_that("build_dct_basis rejects non-positive counts before sequence construction", {
  expect_error(
    build_dct_basis(n_time = 8L, k = 0L),
    class = "fmrilatent_error_invalid_count"
  )
  expect_error(
    build_dct_basis(n_time = 0L, k = 1L),
    class = "fmrilatent_error_invalid_count"
  )
  expect_error(
    build_dct_basis(n_time = 8L, k = -1L),
    class = "fmrilatent_error_invalid_count"
  )
})

test_that("build_dct_basis rejects non-integer counts", {
  expect_error(
    build_dct_basis(n_time = 8.5, k = 2L),
    class = "fmrilatent_error_invalid_count"
  )
  expect_error(
    build_dct_basis(n_time = 8L, k = 2.5),
    class = "fmrilatent_error_invalid_count"
  )
})

test_that("dct_basis_handle rejects invalid dimensions immediately", {
  expect_error(
    dct_basis_handle(n_time = 0L, k = 1L),
    class = "fmrilatent_error_invalid_count"
  )
  expect_error(
    dct_basis_handle(n_time = 4L, k = 0L),
    class = "fmrilatent_error_invalid_count"
  )
  expect_error(
    dct_basis_handle(n_time = 4L, k = 5L),
    class = "fmrilatent_error_invalid_count"
  )
})
