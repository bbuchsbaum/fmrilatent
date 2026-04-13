library(testthat)
library(Matrix)

.skip_if_no_neurosurf_st <- function() {
  skip_if_not_installed("neurosurf")
}

.make_surface_geometry_st <- function() {
  .skip_if_no_neurosurf_st()
  neurosurf::example_surface_geometry()
}

.make_surface_template_st <- function() {
  geom <- .make_surface_geometry_st()
  loadings <- Matrix(matrix(
    c(1, 0,
      0, 1,
      1, 1),
    nrow = 3, byrow = TRUE
  ), sparse = FALSE)
  surface_basis_template(
    geometry = geom,
    loadings = loadings,
    support = 1:3,
    roughness = diag(2),
    measure = c(1, 2, 3),
    label = "surface_test",
    meta = list(family = "surface_test")
  )
}

test_that("SurfaceBasisTemplate satisfies the template protocol", {
  tmpl <- .make_surface_template_st()
  payload <- .template_coordinate_payload(
    raw_loadings = template_loadings(tmpl),
    measure = template_measure(tmpl),
    default_measure = "null"
  )

  expect_true(is_surface_template(tmpl))
  expect_true(is_template(tmpl))
  expect_equal(template_rank(tmpl), 2L)
  expect_equal(as.integer(template_support(tmpl)), 1:3)
  expect_type(template_measure(tmpl), "double")
  expect_true(methods::is(template_domain(tmpl), "SurfaceGeometry"))
  expect_equal(as.matrix(template_loadings(tmpl)),
               matrix(c(1, 0,
                        0, 1,
                        1, 1), nrow = 3, byrow = TRUE))
  expect_equal(as.matrix(template_roughness(tmpl, "raw")), diag(2), tolerance = 1e-8)
  expect_equal(
    as.matrix(template_roughness(tmpl, "analysis")),
    .transform_quadratic_form(diag(2), payload$analysis_transform, coordinates = "analysis"),
    tolerance = 1e-8
  )
  expect_error(template_mask(tmpl), "not defined")
})

test_that("SurfaceBasisTemplate decoder and projection agree with direct linear algebra", {
  tmpl <- .make_surface_template_st()
  payload <- .template_coordinate_payload(
    raw_loadings = template_loadings(tmpl),
    measure = template_measure(tmpl),
    default_measure = "null"
  )
  B <- as.matrix(payload$analysis_loadings)
  X <- matrix(c(
    1, 2, 3,
    4, 5, 6
  ), nrow = 2, byrow = TRUE)

  dec <- basis_decoder(tmpl)
  expect_equal(.materialize_linear_map(dec), B, tolerance = 1e-8)

  proj <- template_project(tmpl, X)
  expected <- .template_projection_payload(
    data = X,
    raw_loadings = template_loadings(tmpl),
    measure = template_measure(tmpl),
    default_measure = "null"
  )$coefficients

  expect_equal(as.matrix(proj$coefficients), as.matrix(expected), tolerance = 1e-8)
  expect_equal(proj$offset, numeric(0))
})

test_that("ImplicitLatent can wrap decoded surface outputs", {
  geom <- .make_surface_geometry_st()
  support <- 1:3
  full <- matrix(c(
    1, 2, 3,
    4, 5, 6
  ), nrow = 2, byrow = TRUE)

  il <- implicit_latent(
    coeff = list(full = full),
    decoder = function(time_idx = NULL, roi_mask = NULL, levels_keep = NULL, ...) {
      rec <- full
      if (!is.null(time_idx)) {
        rec <- rec[time_idx, , drop = FALSE]
      }
      rec
    },
    meta = list(family = "surface_implicit_test"),
    domain = geom,
    support = support
  )

  wrapped_matrix <- wrap_decoded(il, reconstruct_matrix(il))
  wrapped_vector <- wrap_decoded(il, c(10, 20, 30))

  expect_equal(latent_domain(il), geom)
  expect_equal(as.integer(latent_support(il)), support)
  expect_s4_class(wrapped_matrix, "NeuroSurfaceVector")
  expect_s4_class(wrapped_vector, "NeuroSurface")
  expect_equal(as.integer(neuroim2::indices(wrapped_matrix)), support)
  expect_equal(as.integer(neuroim2::indices(wrapped_vector)), support)
  expect_equal(as.matrix(wrapped_matrix), t(full), tolerance = 1e-8)
})
