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

make_shared_structure_awpt_template <- function() {
  mask_vol <- LogicalNeuroVol(array(TRUE, dim = c(3, 1, 1)), NeuroSpace(c(3, 1, 1)))
  reduction <- make_cluster_reduction(mask_vol, c(1L, 1L, 2L))
  loadings <- Matrix(matrix(c(1, 0, 0, 1, 1, 1), nrow = 3, byrow = TRUE), sparse = FALSE)
  structure(
    list(
      loadings = loadings,
      gram_factor = Matrix::Cholesky(Matrix::crossprod(loadings), perm = TRUE),
      roughness = Matrix::Diagonal(2),
      reduction = reduction,
      basis_spec = basis_awpt_wavelet(scales = c(1, 2)),
      atoms = data.frame(
        col_id = 1:2,
        cluster_id = c(1L, 2L),
        scale = c(1, 2),
        scale_index = 1:2,
        roughness_weight = c(1, 1)
      ),
      center = FALSE,
      meta = list(family = "awpt_wavelet", k = 2L)
    ),
    class = "AWPTBasisTemplate"
  )
}

make_shared_structure_hierarchical_template <- function() {
  mask_vol <- LogicalNeuroVol(array(TRUE, dim = c(3, 1, 1)), NeuroSpace(c(3, 1, 1)))
  loadings <- Matrix(matrix(c(1, 0, 0, 1, 1, 1), nrow = 3, byrow = TRUE), sparse = FALSE)
  new(
    "HierarchicalBasisTemplate",
    mask = mask_vol,
    space = NeuroSpace(c(3, 1, 1, 1)),
    levels = list(c(1L, 1L, 2L)),
    parents = list(),
    loadings = loadings,
    gram_factor = Matrix::Cholesky(Matrix::crossprod(loadings), perm = TRUE),
    atoms = data.frame(
      col_id = 1:2,
      level = 1L,
      parcel_id = c(1L, 2L),
      parent_id = NA_integer_,
      mode = "test",
      label = c("h1", "h2")
    ),
    meta = list(family = "hierarchical", k_per_level = 2L)
  )
}

test_that("matrix-like finite checks preserve sparse representation", {
  sparse_ok <- Matrix::sparseMatrix(i = 1L, j = 1L, x = 1, dims = c(100L, 100L))
  sparse_bad <- sparse_ok
  sparse_bad@x[1L] <- Inf

  expect_true(fmrilatent:::.matrix_like_all_finite(sparse_ok))
  expect_false(fmrilatent:::.matrix_like_all_finite(sparse_bad))
  expect_s4_class(fmrilatent:::.shared_matrix_like(sparse_ok, "sparse_ok"), "sparseMatrix")
  expect_error(
    fmrilatent:::.shared_matrix_like(sparse_bad, "sparse_bad"),
    class = "fmrilatent_error_invalid_value"
  )
})

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

test_that("shared component contracts preserve dimensional invariants across random shapes", {
  set.seed(40)
  for (iter in seq_len(10L)) {
    n_features <- sample(3:8, 1L)
    n_components <- sample(1:min(4L, n_features), 1L)
    loadings <- Matrix(
      matrix(rnorm(n_features * n_components), nrow = n_features),
      sparse = iter %% 2L == 0L
    )
    support <- seq_len(n_features)
    measure <- runif(n_features, min = 0.1, max = 2)

    contract <- shared_component_contract(
      loadings,
      support = support,
      measure = measure,
      family = paste0("random_", iter)
    )

    expect_equal(contract$n_features, n_features)
    expect_equal(contract$n_components, n_components)
    expect_identical(validate_shared_component_contract(contract), contract)
  }
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

test_that("shared references can bypass caching and clear session state", {
  shared_reference_clear()
  counter <- new.env(parent = emptyenv())
  counter$n <- 0L
  ref <- shared_reference(
    id = "uncached-ref",
    kind = "loadings",
    materialize = function() {
      counter$n <- counter$n + 1L
      matrix(counter$n, nrow = 1, ncol = 1)
    }
  )

  expect_equal(resolve_shared_reference(ref, cache = FALSE), matrix(1, 1, 1))
  expect_equal(resolve_shared_reference(ref, cache = FALSE), matrix(2, 1, 1))
  expect_false(shared_reference_info(ref)$has_cached_value)
  expect_equal(resolve_shared_reference(ref, cache = TRUE), matrix(3, 1, 1))
  expect_true(shared_reference_info(ref)$has_cached_value)
  shared_reference_clear()
  expect_false(shared_reference_info(ref)$has_cached_value)
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

test_that("group-plus-delta materialization is affine in scale", {
  group <- Matrix(matrix(c(1, 2, 3, 4), nrow = 2), sparse = FALSE)
  delta <- Matrix(matrix(c(0.5, -0.25, 1, -2), nrow = 2), sparse = TRUE)
  gd1 <- group_delta_loadings(group, delta, scale = 1)
  gd2 <- group_delta_loadings(group, delta, scale = 2)

  expect_equal(
    as.matrix(gd2) - as.matrix(gd1),
    as.matrix(delta),
    tolerance = 1e-12
  )
})

test_that("shared temporal specs and event dictionaries are executable descriptors", {
  dct <- shared_temporal_spec("dct", n_time = 8L, rank = 3L)
  B <- materialize_shared_temporal_spec(dct)
  expect_equal(dim(B), c(8L, 3L))
  expect_equal(crossprod(B), diag(3), tolerance = 1e-10)

  slepian <- shared_temporal_spec(
    "slepian",
    n_time = 8L,
    rank = 2L,
    params = list(tr = 1, bandwidth = 0.2)
  )
  B_slepian <- materialize_shared_temporal_spec(slepian)
  expect_equal(dim(B_slepian), c(8L, 2L))
  expect_equal(crossprod(B_slepian), diag(2), tolerance = 1e-10)

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

test_that("shared event rendering is additive under event decomposition", {
  dictionary <- shared_event_dictionary(
    shapes = list(impulse = 1, ramp = c(1, 0.5)),
    n_time = 4L
  )
  event_a <- data.frame(atom = 1L, time = 2L, amplitude = 3, shape_id = "ramp")
  event_b <- data.frame(atom = 1L, time = 3L, amplitude = 2, shape_id = "impulse")
  together <- rbind(event_a, event_b)

  expect_equal(
    as.matrix(render_shared_events(dictionary, together, n_atoms = 1L)),
    as.matrix(render_shared_events(dictionary, event_a, n_atoms = 1L)) +
      as.matrix(render_shared_events(dictionary, event_b, n_atoms = 1L))
  )
  expect_error(
    render_shared_events(dictionary, data.frame(atom = 2L, time = 1L, amplitude = 1), n_atoms = 1L),
    "out-of-range"
  )
})

test_that("shared event rendering handles many events without changing values", {
  dictionary <- shared_event_dictionary(
    shapes = list(impulse = 1, box3 = c(1, 1, 1)),
    n_time = 30L
  )
  events <- data.frame(
    atom = rep(1:3, length.out = 120L),
    time = rep(1:30, length.out = 120L),
    amplitude = seq_len(120L) / 10,
    shape_id = rep(c("impulse", "box3"), length.out = 120L)
  )

  rendered <- render_shared_events(dictionary, events, n_atoms = 3L)
  expected <- matrix(0, nrow = 3L, ncol = 30L)
  for (idx in seq_len(nrow(events))) {
    shape <- dictionary$shapes[[events$shape_id[[idx]]]]
    cols <- seq.int(events$time[[idx]], length.out = length(shape))
    keep <- cols <= ncol(expected)
    expected[events$atom[[idx]], cols[keep]] <-
      expected[events$atom[[idx]], cols[keep]] + events$amplitude[[idx]] * shape[keep]
  }

  expect_equal(as.matrix(rendered), expected)
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
  expect_identical(handoff$components[[1]], contract)
  expect_s3_class(handoff$templates[[1]], "ParcelBasisTemplate")
  expect_equal(template_loadings(handoff$templates[[1]]), template_loadings(tmpl))
  expect_equal(handoff$template_manifests[[1]]$rank, manifest$rank)
  expect_identical(validate_neuroarchive_handoff_contract(handoff), handoff)

  bad_handoff <- handoff
  bad_handoff$archive_locator <- "archive.h5"
  expect_error(validate_neuroarchive_handoff_contract(bad_handoff), "persistent archive locators")
})

test_that("shared template protocol covers parcel, AWPT, and hierarchical families", {
  templates <- list(
    parcel = make_shared_structure_template(),
    awpt = make_shared_structure_awpt_template(),
    hierarchical = make_shared_structure_hierarchical_template()
  )
  manifests <- lapply(templates, validate_template_protocol)

  expect_equal(vapply(manifests, `[[`, integer(1), "rank"), c(parcel = 2L, awpt = 2L, hierarchical = 2L))
  expect_equal(vapply(manifests, `[[`, integer(1), "n_features"), c(parcel = 4L, awpt = 3L, hierarchical = 3L))

  contracts <- Map(
    function(template, family) {
      shared_component_contract(
        template_loadings(template),
        family = family,
        support = template_support(template)
      )
    },
    templates,
    names(templates)
  )
  expect_true(all(vapply(contracts, inherits, logical(1), "SharedComponentContract")))
})

test_that("template and handoff validators reject wrong-contract objects", {
  tmpl <- make_shared_structure_template()
  bad_tmpl <- tmpl
  bad_tmpl$decoder_map <- as_portable_linear_map(matrix(1, nrow = 3L, ncol = 2L))

  expect_error(
    validate_template_protocol(bad_tmpl),
    "dimensions do not match template loadings"
  )
  expect_error(
    neuroarchive_handoff_contract(references = list(list(id = "not-shared"))),
    "SharedReference"
  )

  bad_component <- shared_component_contract(template_loadings(tmpl), support = template_support(tmpl))
  bad_component$contract_version <- 2L
  expect_error(
    neuroarchive_handoff_contract(components = list(bad_component)),
    "unsupported shared component contract_version"
  )
})
