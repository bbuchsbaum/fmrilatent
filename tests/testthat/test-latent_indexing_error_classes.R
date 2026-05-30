# Structured-condition coverage for latent_indexing.R
# (bd-01KSQQNYTQCW36WNEYPYMX8188). Bare stop() calls in the indexing methods
# were converted to classed .encoder_cli_abort() conditions. Triggers are
# dependency-free.

library(testthat)

make_lvec <- function(nx = 3, ny = 3, nz = 2, nt = 5, k = 2) {
  mask_arr <- array(TRUE, dim = c(nx, ny, nz))
  mask_vol <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(c(nx, ny, nz)))
  space <- neuroim2::NeuroSpace(c(nx, ny, nz, nt))
  set.seed(42)
  basis <- Matrix::Matrix(matrix(rnorm(nt * k), nt, k), sparse = FALSE)
  loadings <- Matrix::Matrix(matrix(rnorm(sum(mask_arr) * k), sum(mask_arr), k), sparse = FALSE)
  LatentNeuroVec(basis = basis, loadings = loadings, space = space,
                 mask = mask_vol, offset = numeric(0), label = "idx-test")
}

test_that("[[ guards index validity with classed conditions", {
  lvec <- make_lvec()
  expect_error(lvec[[c(1, 2)]], class = "fmrilatent_error_invalid_index")  # not scalar
  expect_error(lvec[[1.5]], class = "fmrilatent_error_invalid_index")      # not integer
  expect_error(lvec[[100]], class = "fmrilatent_error_invalid_index")      # out of range
  # message text preserved for the existing substring matchers
  expect_error(lvec[[c(1, 2)]], "single finite number")
})

test_that("linear_access guards index bounds with classed conditions", {
  lvec <- make_lvec()
  expect_error(linear_access(lvec, NA_integer_), class = "fmrilatent_error_invalid_index")
  expect_error(linear_access(lvec, 1e9), class = "fmrilatent_error_invalid_index")
  expect_error(linear_access(lvec, 1e9), "out of bounds")
})

test_that("matricized_access rejects a non-2-column matrix index", {
  lvec <- make_lvec()
  expect_error(
    matricized_access(lvec, matrix(1, nrow = 1, ncol = 3)),
    class = "fmrilatent_error_invalid_argument"
  )
})

test_that("[ guards subscript range with classed conditions", {
  lvec <- make_lvec()
  expect_error(lvec[1, 1, 1, 99], class = "fmrilatent_error_invalid_index")
})
