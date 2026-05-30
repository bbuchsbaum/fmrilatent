library(testthat)
library(Matrix)

.skip_if_no_neurosurf_bval <- function() {
  skip_if(
    any(commandArgs() == "-f"),
    "neurosurf surface examples abort under this local R/BATCH setup"
  )
  skip_if_not_installed("neurosurf")
}

.make_surface_bval <- function(hemi = c("left", "right"), basis = NULL) {
  hemi <- match.arg(hemi)
  .skip_if_no_neurosurf_bval()
  geom <- neurosurf::example_surface_geometry()
  geom@hemi <- hemi
  if (is.null(basis)) {
    basis <- Matrix(matrix(c(
      1, 0,
      0, 1,
      1, 1
    ), nrow = 3, byrow = TRUE), sparse = FALSE)
  }
  loadings <- Matrix(matrix(c(
    1, 0,
    0, 1,
    1, 1
  ), nrow = 3, byrow = TRUE), sparse = FALSE)
  LatentNeuroSurfaceVector(
    basis = basis,
    loadings = loadings,
    geometry = geom,
    support = 1:3,
    offset = c(0, 0.5, -0.5),
    meta = list(family = paste0("surface_", hemi))
  )
}

# A basis with the SAME dimensions (3 x 2) as the default but DIFFERENT values.
.different_basis_bval <- function() {
  Matrix(matrix(c(
    2, 0,
    0, 3,
    1, 4
  ), nrow = 3, byrow = TRUE), sparse = FALSE)
}

test_that("BilatLatentNeuroSurfaceVector accepts a matched left/right basis pair", {
  left <- .make_surface_bval("left")
  right <- .make_surface_bval("right")
  expect_error(
    BilatLatentNeuroSurfaceVector(left, right),
    NA
  )
})

test_that("BilatLatentNeuroSurfaceVector rejects same-dimension but different bases", {
  left <- .make_surface_bval("left")
  right <- .make_surface_bval("right", basis = .different_basis_bval())

  # Sanity: the bases share dimensions but differ in value, so the old
  # dimension-only check would have silently accepted this pair.
  expect_identical(dim(as.matrix(basis(left))), dim(as.matrix(basis(right))))
  expect_false(isTRUE(all.equal(as.matrix(basis(left)), as.matrix(basis(right)))))

  expect_error(
    BilatLatentNeuroSurfaceVector(left, right),
    "basis matrices must be equal"
  )
  expect_error(
    BilatLatentNeuroSurfaceVector(left, right),
    class = "fmrilatent_error_dim"
  )
})

test_that("setValidity catches a directly-constructed mismatched bilateral object", {
  left <- .make_surface_bval("left")
  right <- .make_surface_bval("right", basis = .different_basis_bval())

  # Bypass the constructor; validity must still reject on basis inequality.
  expect_error(
    methods::new("BilatLatentNeuroSurfaceVector", left = left, right = right,
                 label = "", meta = list()),
    "basis matrices must be equal"
  )

  # validObject() on a forcibly-built object should also fail.
  obj <- methods::new("BilatLatentNeuroSurfaceVector", left = left, right = left,
                      label = "", meta = list())
  obj@right <- right
  expect_error(methods::validObject(obj), "basis matrices must be equal")
})

test_that("BlockLatentNeuroVector still rejects same-dim mismatched block bases", {
  .skip_if_no_neurosurf_bval()
  left <- .make_surface_bval("left")
  right_bad <- .make_surface_bval("right", basis = .different_basis_bval())

  # Both surface bases are 3 x 2; only their values differ.
  expect_identical(
    dim(as.matrix(basis(left))),
    dim(as.matrix(basis(right_bad)))
  )

  expect_error(
    BlockLatentNeuroVector(list(a = left, b = right_bad)),
    "basis matrix must match the shared block basis"
  )
})
