library(testthat)
library(Matrix)

.skip_if_no_neurosurf_lsv <- function() {
  skip_if_not_installed("neurosurf")
}

.make_surface_geometry_lsv <- function() {
  .skip_if_no_neurosurf_lsv()
  neurosurf::example_surface_geometry()
}

.make_latent_surface_vector_lsv <- function() {
  geom <- .make_surface_geometry_lsv()
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
    offset = c(0.5, 0, -0.5),
    meta = list(family = "surface_explicit_test")
  )
}

test_that("LatentNeuroSurfaceVector reconstructs and wraps surface outputs", {
  lv <- .make_latent_surface_vector_lsv()
  mat <- as.matrix(lv)

  expect_true(is_explicit_latent(lv))
  expect_equal(latent_meta(lv)$family, "surface_explicit_test")
  expect_true(methods::is(latent_domain(lv), "SurfaceGeometry"))
  expect_equal(as.integer(latent_support(lv)), 1:3)
  expect_equal(reconstruct_matrix(lv), mat)

  wrapped_matrix <- wrap_decoded(lv, mat)
  wrapped_vector <- wrap_decoded(lv, mat[1, ])
  expect_s4_class(wrapped_matrix, "NeuroSurfaceVector")
  expect_s4_class(wrapped_vector, "NeuroSurface")
  expect_equal(as.integer(neuroim2::indices(wrapped_matrix)), 1:3)
  expect_equal(as.matrix(wrapped_matrix), t(mat), tolerance = 1e-8)
})

test_that("LatentNeuroSurfaceVector supports coefficient decoding and covariance pushforward", {
  lv <- .make_latent_surface_vector_lsv()
  L <- as.matrix(loadings(lv))
  gamma <- c(2, -1)
  Sigma <- matrix(c(2, 0.5,
                    0.5, 1), nrow = 2, byrow = TRUE)

  expect_equal(coef_time(lv), as.matrix(basis(lv)))
  expect_equal(coef_metric(lv, "analysis"), diag(2))
  expect_equal(
    decode_coefficients(lv, gamma, space = "native"),
    as.vector(L %*% matrix(gamma, ncol = 1)),
    tolerance = 1e-8
  )
  expect_equal(
    decode_covariance(lv, Sigma, space = "native", diag_only = FALSE),
    L %*% Sigma %*% t(L),
    tolerance = 1e-8
  )
})

test_that("LatentNeuroSurfaceVector respects support-level ROI masks", {
  lv <- .make_latent_surface_vector_lsv()
  mat <- reconstruct_matrix(lv)
  roi <- c(TRUE, FALSE, TRUE)

  expect_equal(
    reconstruct_matrix(lv, roi_mask = roi),
    mat[, roi, drop = FALSE],
    tolerance = 1e-8
  )
  expect_error(reconstruct_array(lv), "not defined for surface latent objects")
})
