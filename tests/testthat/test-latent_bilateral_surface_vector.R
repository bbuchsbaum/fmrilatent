library(testthat)
library(Matrix)

.skip_if_no_neurosurf_blsv <- function() {
  skip_if_not_installed("neurosurf")
}

.make_unilateral_blsv <- function(hemi = c("left", "right")) {
  hemi <- match.arg(hemi)
  .skip_if_no_neurosurf_blsv()
  geom <- neurosurf::example_surface_geometry()
  geom@hemi <- hemi
  basis <- Matrix(matrix(c(
    1, 0,
    0, 1,
    1, 1
  ), nrow = 3, byrow = TRUE), sparse = FALSE)
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

test_that("BilatLatentNeuroSurfaceVector reconstructs bilateral outputs", {
  left <- .make_unilateral_blsv("left")
  right <- .make_unilateral_blsv("right")
  bilat <- BilatLatentNeuroSurfaceVector(left, right, meta = list(family = "bilat_test"))

  mat <- reconstruct_matrix(bilat)

  expect_true(is_explicit_latent(bilat))
  expect_equal(latent_meta(bilat)$family, "bilat_test")
  expect_equal(ncol(mat), 6L)
  expect_equal(mat, cbind(reconstruct_matrix(left), reconstruct_matrix(right)), tolerance = 1e-8)

  wrapped <- wrap_decoded(bilat, mat)
  expect_s4_class(wrapped, "BilatNeuroSurfaceVector")
  expect_equal(as.matrix(as.matrix(wrapped)), t(mat), tolerance = 1e-8)
})

test_that("BilatLatentNeuroSurfaceVector supports decode and ROI splitting", {
  left <- .make_unilateral_blsv("left")
  right <- .make_unilateral_blsv("right")
  bilat <- BilatLatentNeuroSurfaceVector(left, right)
  L <- as.matrix(loadings(bilat))
  gamma <- c(2, -1)
  Sigma <- matrix(c(2, 0.25,
                    0.25, 1), nrow = 2, byrow = TRUE)

  expect_equal(
    decode_coefficients(bilat, gamma, space = "native"),
    as.vector(L %*% matrix(gamma, ncol = 1)),
    tolerance = 1e-8
  )
  expect_equal(
    decode_covariance(bilat, Sigma, space = "native", diag_only = FALSE),
    L %*% Sigma %*% t(L),
    tolerance = 1e-8
  )

  roi <- list(left = c(TRUE, FALSE, TRUE), right = c(FALSE, TRUE, FALSE))
  mat <- reconstruct_matrix(bilat)
  expect_equal(
    reconstruct_matrix(bilat, roi_mask = roi),
    cbind(mat[, c(1, 3), drop = FALSE], mat[, 5, drop = FALSE]),
    tolerance = 1e-8
  )
})
