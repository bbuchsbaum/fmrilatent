surface_contracts_enabled <- function() {
  check_path <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  any(commandArgs() == "-f") ||
    grepl("[.]Rcheck($|/)", check_path) ||
    identical(Sys.getenv("FMRILATENT_RUN_SURFACE_CONTRACT_TESTS"), "true")
}

skip_unless_surface_contracts_enabled <- function() {
  skip_if_not(
    surface_contracts_enabled(),
    "surface contract doubles run only in R CMD check or when explicitly enabled"
  )
}

fake_surface_geometry <- function(n_nodes = 4L) {
  if (!methods::isClass("SurfaceGeometry")) {
    methods::setClass("SurfaceGeometry", slots = c(label = "character"))
  }
  slots <- names(methods::getClass("SurfaceGeometry")@slots)
  if (all(c("mesh", "graph", "hemi", "label", "surf_to_world") %in% slots)) {
    mesh <- structure(
      list(
        vb = matrix(0, nrow = 4L, ncol = n_nodes),
        it = matrix(c(1L, 2L, 3L, 1L, 3L, 4L), nrow = 3L)
      ),
      class = c("mesh3d", "shape3d")
    )
    graph <- getExportedValue("igraph", "make_empty_graph")(
      n = n_nodes,
      directed = FALSE
    )
    geom <- methods::new(
      "SurfaceGeometry",
      mesh = mesh,
      graph = graph,
      hemi = "left",
      label = "fmrilatent-test",
      surf_to_world = diag(4)
    )
  } else {
    geom <- methods::new("SurfaceGeometry", label = "fmrilatent-test")
  }
  attr(geom, "fmrilatent.surface_n_nodes") <- as.integer(n_nodes)
  geom
}

make_contract_surface_latent <- function() {
  LatentNeuroSurfaceVector(
    basis = matrix(c(
      1, 0,
      0, 1,
      1, 1
    ), nrow = 3, byrow = TRUE),
    loadings = matrix(c(
      1, 0,
      0, 1,
      1, 1
    ), nrow = 3, byrow = TRUE),
    geometry = fake_surface_geometry(4L),
    support = 1:3,
    offset = c(0.5, 0, -0.5),
    meta = list(family = "surface_contract")
  )
}

test_that("surface latent contracts run without loading neurosurf in R CMD check", {
  skip_unless_surface_contracts_enabled()
  lv <- make_contract_surface_latent()

  expected <- matrix(c(
    1.5, 0.0, 0.5,
    0.5, 1.0, 0.5,
    1.5, 1.0, 1.5
  ), nrow = 3, byrow = TRUE)

  expect_s4_class(lv, "LatentNeuroSurfaceVector")
  expect_true(is_explicit_latent(lv))
  expect_equal(as.integer(latent_support(lv)), 1:3)
  expect_equal(reconstruct_matrix(lv), expected, tolerance = 1e-8)
  expect_equal(
    reconstruct_matrix(lv, roi_mask = c(TRUE, FALSE, TRUE, FALSE)),
    expected[, c(1, 3), drop = FALSE],
    tolerance = 1e-8
  )
})

test_that("surface support validation uses the domain node count helper", {
  skip_unless_surface_contracts_enabled()
  expect_error(
    LatentNeuroSurfaceVector(
      basis = matrix(1, nrow = 2, ncol = 1),
      loadings = matrix(1, nrow = 1, ncol = 1),
      geometry = fake_surface_geometry(4L),
      support = 5L
    ),
    class = "fmrilatent_error_invalid_support"
  )
})

test_that("surface template contracts run without loading neurosurf in R CMD check", {
  skip_unless_surface_contracts_enabled()
  tmpl <- surface_basis_template(
    geometry = fake_surface_geometry(4L),
    loadings = matrix(c(
      1, 0,
      0, 1,
      1, 1
    ), nrow = 3, byrow = TRUE),
    support = 1:3,
    measure = c(1, 2, 3),
    roughness = diag(2),
    meta = list(family = "surface_contract_template")
  )

  X <- matrix(c(
    1, 2, 3,
    4, 5, 6
  ), nrow = 2, byrow = TRUE)
  proj <- template_project(tmpl, X)

  expect_true(is_surface_template(tmpl))
  expect_equal(template_rank(tmpl), 2L)
  expect_equal(as.integer(template_support(tmpl)), 1:3)
  expect_equal(nrow(as.matrix(proj$coefficients)), 2L)
  expect_equal(ncol(as.matrix(proj$coefficients)), 2L)
  expect_equal(.materialize_linear_map(basis_decoder(tmpl)),
               as.matrix(fmrilatent:::.template_coordinate_payload(
                 raw_loadings = template_loadings(tmpl),
                 measure = template_measure(tmpl),
                 default_measure = "null"
               )$analysis_loadings),
               tolerance = 1e-8)
})
