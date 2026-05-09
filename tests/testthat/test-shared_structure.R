library(testthat)
library(Matrix)
library(neuroim2)

make_shared_structure_template <- function() {
  mask_vol <- LogicalNeuroVol(array(TRUE, dim = c(2, 2, 1)), NeuroSpace(c(2, 2, 1)))
  loadings <- Matrix(matrix(
    c(1, 0,
      0, 1,
      1, 1,
      0, 1),
    nrow = 4, byrow = TRUE
  ), sparse = FALSE)
  gram_factor <- Matrix::Cholesky(Matrix::crossprod(loadings), perm = TRUE)
  reduction <- make_cluster_reduction(mask_vol, c(1L, 1L, 2L, 2L))
  structure(
    list(
      loadings = loadings,
      gram_factor = gram_factor,
      reduction = reduction,
      basis_spec = basis_slepian(k = 2),
      center = FALSE,
      meta = list(
        family = "spec_slepian",
        k = 2L,
        ridge = 1e-8,
        label_map = list(A = 1L, B = 2L),
        cluster_map = list(`1` = 1:2, `2` = 3:4)
      )
    ),
    class = "ParcelBasisTemplate"
  )
}

test_that("shared component contracts validate dimensions, support, and measures", {
  loadings <- Matrix(matrix(c(1, 0, 0, 1, 1, 1), nrow = 3, byrow = TRUE), sparse = FALSE)
  contract <- shared_component_contract(
    loadings,
    family = "test_basis",
    domain_id = "toy_domain",
    support = c("a", "b", "c"),
    measure = c(1, 2, 1)
  )

  expect_s3_class(contract, "SharedComponentContract")
  expect_equal(contract$n_features, 3L)
  expect_equal(contract$n_components, 2L)
  expect_equal(contract$domain_id, "toy_domain")
  expect_true(nzchar(contract$digest))
  expect_identical(validate_shared_component_contract(contract), contract)

  expect_error(
    shared_component_contract(loadings, support = c("a", "b")),
    "support cardinality 2 does not match loadings rows 3"
  )
  expect_false(validate_shared_component_contract(structure(list(), class = "SharedComponentContract"), error = FALSE))
})

test_that("shared references stay in-session and materialize lazily", {
  shared_reference_clear()
  counter <- new.env(parent = emptyenv())
  counter$n <- 0L
  ref <- shared_reference(
    id = "toy-ref",
    kind = "loadings",
    materialize = function() {
      counter$n <- counter$n + 1L
      matrix(1:4, nrow = 2)
    }
  )

  expect_true(is_shared_reference(ref))
  expect_false(shared_reference_info(ref)$has_cached_value)
  expect_equal(resolve_shared_reference(ref), matrix(1:4, nrow = 2))
  expect_equal(resolve_shared_reference(ref), matrix(1:4, nrow = 2))
  expect_equal(counter$n, 1L)
  expect_true(shared_reference_info(ref)$has_cached_value)
  expect_equal(shared_reference_info(ref)$persistence, "session")
})

test_that("group-plus-delta loadings materialize and feed component contracts", {
  group <- Matrix(matrix(c(1, 0, 0, 1, 1, 1), nrow = 3, byrow = TRUE), sparse = FALSE)
  delta <- Matrix(matrix(c(0, 1, 1, 0, 0, -1), nrow = 3, byrow = TRUE), sparse = FALSE)
  gd <- group_delta_loadings(group, delta, scale = 0.5)

  expect_s3_class(gd, "GroupDeltaLoadings")
  expect_equal(as.matrix(gd), as.matrix(group + 0.5 * delta))

  contract <- shared_component_contract(gd, family = "group_delta", support = 1:3)
  expect_equal(contract$n_features, 3L)
  expect_equal(contract$n_components, 2L)
  expect_error(group_delta_loadings(group, delta[1:2, ]), "identical dimensions")
})

test_that("shared temporal specs and event dictionaries are executable descriptors", {
  dct <- shared_temporal_spec("dct", n_time = 8L, rank = 3L)
  B <- materialize_shared_temporal_spec(dct)
  expect_equal(dim(B), c(8L, 3L))
  expect_equal(crossprod(B), diag(3), tolerance = 1e-10)

  custom_basis <- matrix(c(1, 0, 0, 1, 1, 1), nrow = 3, byrow = TRUE)
  custom <- shared_temporal_spec("custom", basis = custom_basis)
  expect_equal(materialize_shared_temporal_spec(custom), custom_basis)

  dictionary <- shared_event_dictionary(
    shapes = list(impulse = 1, box2 = c(1, 1)),
    n_time = 5L
  )
  events <- data.frame(
    atom = c(1L, 2L),
    time = c(2L, 4L),
    amplitude = c(3, 2),
    shape_id = c("impulse", "box2")
  )
  rendered <- render_shared_events(dictionary, events, n_atoms = 2L)
  expected <- matrix(0, nrow = 2, ncol = 5)
  expected[1, 2] <- 3
  expected[2, 4:5] <- 2

  expect_equal(as.matrix(rendered), expected)
  expect_error(render_shared_events(dictionary, transform(events, shape_id = "missing"), 2L), "unknown event shape_id")
})

test_that("template protocol and neuroarchive handoff enforce the ownership boundary", {
  tmpl <- make_shared_structure_template()
  manifest <- validate_template_protocol(tmpl)
  contract <- shared_component_contract(template_loadings(tmpl), support = template_support(tmpl))
  handoff <- neuroarchive_handoff_contract(
    components = list(contract),
    templates = list(tmpl),
    references = list(shared_reference(template_loadings(tmpl), kind = "loadings")),
    meta = list(scope = "unit-test")
  )

  expect_equal(manifest$rank, 2L)
  expect_equal(manifest$n_features, 4L)
  expect_s3_class(handoff, "NeuroarchiveHandoffContract")
  expect_true("projection and reconstruction behavior" %in% handoff$fmrilatent_owns)
  expect_true("HDF5 artifacts" %in% handoff$neuroarchive_owns)
  expect_false(handoff$persistent)
  expect_null(handoff$archive_locator)
  expect_identical(validate_neuroarchive_handoff_contract(handoff), handoff)

  bad_handoff <- handoff
  bad_handoff$archive_locator <- "archive.h5"
  expect_error(validate_neuroarchive_handoff_contract(bad_handoff), "persistent archive locators")
})
