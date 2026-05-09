# Tests for the show() method on LatentNeuroVec.

skip_if_no_neuroim2 <- function() {
  testthat::skip_if_not_installed("neuroim2")
}

build_small_lvec <- function(label = "show-test", basis_kind = c("dense", "handle")) {
  basis_kind <- match.arg(basis_kind)
  set.seed(42)

  n_time <- 8L
  k      <- 3L
  mask   <- neuroim2::LogicalNeuroVol(
    array(TRUE, dim = c(3L, 3L, 3L)),
    neuroim2::NeuroSpace(c(3L, 3L, 3L))
  )
  n_vox  <- sum(as.array(mask))
  space  <- neuroim2::NeuroSpace(c(3L, 3L, 3L, n_time))

  basis <- if (basis_kind == "dense") {
    Matrix::Matrix(matrix(rnorm(n_time * k), n_time, k), sparse = FALSE)
  } else {
    dct_basis_handle(n_time = n_time, k = k, id = paste0("show-test-handle-", label))
  }

  loadings <- Matrix::Matrix(matrix(rnorm(n_vox * k), n_vox, k), sparse = FALSE)
  LatentNeuroVec(basis = basis, loadings = loadings, space = space, mask = mask, label = label)
}

test_that("show() prints a sensible summary for a dense LatentNeuroVec", {
  skip_if_no_neuroim2()

  lvec <- build_small_lvec(label = "show-dense", basis_kind = "dense")
  out <- capture.output(methods::show(lvec))
  text <- paste(out, collapse = "\n")

  expect_match(text, "LatentNeuroVec Object")
  expect_match(text, "Label:.*show-dense")
  expect_match(text, "Spatial:.*3 x 3 x 3")
  expect_match(text, "Temporal:.*8")
  expect_match(text, "Number:.*3")
  # First-component preview must be numeric, not the dgeMatrix repr.
  expect_match(text, "First component, first 5 coeffs:")
  expect_false(grepl("S4 class", text, fixed = TRUE))
  expect_false(grepl("dgeMatrix", text, fixed = TRUE))
  expect_match(text, "Non-zero voxels:.*27")
})

test_that("show() reports handle kind for BasisHandle backed objects", {
  skip_if_no_neuroim2()

  lvec <- build_small_lvec(label = "show-handle", basis_kind = "handle")
  out  <- capture.output(methods::show(lvec))
  text <- paste(out, collapse = "\n")

  expect_match(text, "Basis:.*handle \\(kind=dct, lazy\\)")
  expect_match(text, "Loadings:.*[0-9]")  # dense loadings still report a size
  expect_false(grepl("S4 class", text, fixed = TRUE))
})

test_that("show() suppresses label line when label is empty", {
  skip_if_no_neuroim2()

  lvec <- build_small_lvec(label = "", basis_kind = "dense")
  out  <- capture.output(methods::show(lvec))
  text <- paste(out, collapse = "\n")

  expect_false(grepl("Label:", text, fixed = TRUE))
})

test_that("show() preview adapts when fewer than 5 timepoints are available", {
  skip_if_no_neuroim2()

  set.seed(7)
  n_time <- 3L
  k      <- 2L
  mask   <- neuroim2::LogicalNeuroVol(
    array(TRUE, dim = c(2L, 2L, 2L)),
    neuroim2::NeuroSpace(c(2L, 2L, 2L))
  )
  n_vox <- sum(as.array(mask))
  space <- neuroim2::NeuroSpace(c(2L, 2L, 2L, n_time))
  B <- Matrix::Matrix(matrix(rnorm(n_time * k), n_time, k), sparse = FALSE)
  L <- Matrix::Matrix(matrix(rnorm(n_vox  * k), n_vox,  k), sparse = FALSE)
  lvec <- LatentNeuroVec(B, L, space, mask, label = "tiny")

  out  <- capture.output(methods::show(lvec))
  text <- paste(out, collapse = "\n")

  expect_match(text, "First component, first 3 coeffs:")
  # Trailing ellipsis only when more rows exist than were previewed.
  expect_false(grepl("First component, first 3 coeffs:.*\\.\\.\\.", text))
})
