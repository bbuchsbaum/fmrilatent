library(testthat)
library(Matrix)

.skip_if_no_neurosurf_sawpt <- function() {
  skip_if_not_installed("neurosurf")
}

.make_surface_geom_sawpt <- function() {
  .skip_if_no_neurosurf_sawpt()
  neurosurf::example_surface_geometry()
}

.make_surface_identity_operator_sawpt <- function(geometry, support) {
  .linear_map_from_matrix(
    diag(length(support)),
    source_domain_id = "surface_template",
    target_domain_id = "surface_native",
    provenance = list(
      target_domain = geometry,
      target_support = support
    )
  )
}

test_that("surface AWPT template builds from a mesh Laplacian", {
  geom <- .make_surface_geom_sawpt()
  tmpl <- awpt_surface_basis_template(
    geometry = geom,
    basis_spec = basis_awpt_wavelet(scales = c(0.5, 1)),
    support = 1:4
  )

  expect_true(is_awpt_template(tmpl))
  expect_true(is_template(tmpl))
  expect_true(inherits(tmpl, "SurfaceBasisTemplate"))
  expect_equal(as.integer(template_support(tmpl)), 1:4)
  expect_equal(template_domain(tmpl), geom)
  expect_null(template_measure(tmpl))
  expect_equal(nrow(template_loadings(tmpl)), 4L)
  expect_gt(ncol(template_loadings(tmpl)), 0L)
  expect_equal(dim(template_roughness(tmpl)),
               c(template_rank(tmpl), template_rank(tmpl)))
})

test_that("encode_awpt works on a surface target domain", {
  geom <- .make_surface_geom_sawpt()
  tmpl <- awpt_surface_basis_template(
    geometry = geom,
    basis_spec = basis_awpt_wavelet(scales = c(0.5, 1)),
    support = 1:4
  )
  B <- .materialize_linear_map(basis_decoder(tmpl))
  op <- .make_surface_identity_operator_sawpt(geom, 1:4)

  Z <- matrix(seq_len(2 * ncol(B)), nrow = 2)
  X <- Z %*% t(B)

  fit <- encode_awpt(
    x = X,
    basis_asset = tmpl,
    observation_operator = op,
    domain = geom,
    support = 1:4,
    spatial_lambda = 0,
    temporal_lambda = 0,
    center = FALSE
  )

  expect_true(is_transport_latent(fit))
  expect_equal(latent_domain(fit), geom)
  expect_equal(as.integer(latent_support(fit)), 1:4)
  expect_equal(reconstruct_matrix(fit), X, tolerance = 1e-6)
  expect_s4_class(wrap_decoded(fit, reconstruct_matrix(fit)), "NeuroSurfaceVector")
  expect_equal(
    decode_coefficients(fit, coef_time(fit)[1, ], space = "native"),
    as.vector(X[1, ]),
    tolerance = 1e-6
  )
})

test_that("surface transport ROI subsetting works with support-level logical masks", {
  geom <- .make_surface_geom_sawpt()
  tmpl <- awpt_surface_basis_template(
    geometry = geom,
    basis_spec = basis_awpt_wavelet(scales = c(1)),
    support = 1:4
  )
  B <- .materialize_linear_map(basis_decoder(tmpl))
  op <- .make_surface_identity_operator_sawpt(geom, 1:4)
  Z <- matrix(c(1, -1, 0, 2), nrow = 1)
  X <- Z %*% t(B)

  fit <- encode_awpt(
    x = X,
    basis_asset = tmpl,
    observation_operator = op,
    domain = geom,
    support = 1:4,
    spatial_lambda = 0,
    temporal_lambda = 0,
    center = FALSE
  )

  roi <- c(TRUE, FALSE, TRUE, FALSE)
  expect_equal(
    reconstruct_matrix(fit, roi_mask = roi),
    X[, roi, drop = FALSE],
    tolerance = 1e-6
  )
})

test_that("surface AWPT template_project honors center = TRUE", {
  geom <- .make_surface_geom_sawpt()
  tmpl <- awpt_surface_basis_template(
    geometry = geom,
    basis_spec = basis_awpt_wavelet(scales = c(1)),
    support = 1:4,
    center = TRUE
  )

  X <- matrix(
    c(1, 3, 5, 7,
      2, 4, 6, 8),
    nrow = 2,
    byrow = TRUE
  )
  proj <- template_project(tmpl, X)
  expected <- .template_projection_payload(
    data = X,
    raw_loadings = template_loadings(tmpl),
    measure = template_measure(tmpl),
    center = TRUE,
    default_measure = "null"
  )$coefficients

  expect_equal(proj$offset, colMeans(X))
  expect_equal(as.matrix(proj$coefficients), as.matrix(expected), tolerance = 1e-6)
})
