make_identity_volume <- function(n_time = 3L, component_order = 1:2,
                                 origin = c(0, 0, 0),
                                 support_order = 1:4,
                                 analysis_transform = NULL) {
  mask <- neuroim2::LogicalNeuroVol(
    array(TRUE, dim = c(2, 2, 1)),
    neuroim2::NeuroSpace(c(2, 2, 1), origin = origin)
  )
  scores <- matrix(seq_len(n_time * 2L), nrow = n_time, ncol = 2L)
  loadings <- matrix(c(
    1, 0,
    0, 1,
    1, 1,
    0, 1
  ), nrow = 4, byrow = TRUE)[, component_order, drop = FALSE]
  meta <- list(family = "identity_fixture")
  if (!is.null(analysis_transform)) {
    meta$analysis_transform <- analysis_transform
  }
  value <- LatentNeuroVec(
    basis = Matrix::Matrix(scores[, component_order, drop = FALSE]),
    loadings = Matrix::Matrix(loadings),
    space = neuroim2::NeuroSpace(c(2, 2, 1, n_time), origin = origin),
    mask = mask,
    meta = meta
  )
  value@map@indices <- as.integer(support_order)
  value
}

make_identity_surface <- function(hemi = c("left", "right"),
                                  support_order = 1:3) {
  hemi <- match.arg(hemi)
  geometry <- neurosurf::example_surface_geometry()
  geometry@hemi <- hemi
  LatentNeuroSurfaceVector(
    basis = Matrix::Matrix(diag(2)),
    loadings = Matrix::Matrix(matrix(c(1, 0, 0, 1, 1, 1), nrow = 3,
      byrow = TRUE
    )),
    geometry = geometry,
    support = support_order
  )
}

make_identity_transport_mask <- function() {
  neuroim2::LogicalNeuroVol(
    array(TRUE, dim = c(2, 2, 1)),
    neuroim2::NeuroSpace(c(2, 2, 1))
  )
}

make_identity_transport_template <- function() {
  mask <- make_identity_transport_mask()
  loadings <- Matrix::Matrix(matrix(c(
    1, 0,
    0, 1,
    1, 1,
    0, 1
  ), nrow = 4, byrow = TRUE))
  structure(
    list(
      loadings = loadings,
      gram_factor = Matrix::Cholesky(Matrix::crossprod(loadings), perm = TRUE),
      reduction = make_cluster_reduction(mask, c(1L, 1L, 2L, 2L)),
      basis_spec = basis_slepian(k = 2),
      center = FALSE,
      meta = list(family = "identity_transport_template", k = 2L)
    ),
    class = "ParcelBasisTemplate"
  )
}

make_identity_transport_operator <- function() {
  .linear_map_from_matrix(
    diag(4),
    source_domain_id = "identity_template_field",
    target_domain_id = "identity_native_field",
    provenance = list(
      operator = "identity_transport",
      target_mask = make_identity_transport_mask()
    )
  )
}

test_that("explicit identity excludes observation values and count", {
  anchor <- make_identity_volume(n_time = 3L)
  task <- make_identity_volume(n_time = 7L)
  task@basis[] <- rev(seq_along(task@basis))

  anchor_id <- latent_space_id(anchor)
  task_id <- latent_space_id(task)

  expect_s3_class(anchor_id, "fmrilatent_latent_space_id")
  expect_identical(anchor_id$contract, "fmrilatent.latent-space-id.v1")
  expect_identical(anchor_id$coordinate_id, task_id$coordinate_id)
  expect_identical(anchor_id$decoder_domain_id, task_id$decoder_domain_id)
  expect_identical(anchor_id$decoder_id, task_id$decoder_id)
  expect_identical(anchor_id$dimension, 2L)
})

test_that("explicit identity separates coordinate and decoder-domain layers", {
  baseline <- make_identity_volume()
  reordered_components <- make_identity_volume(component_order = 2:1)
  changed_domain <- make_identity_volume(origin = c(5, 0, 0))
  changed_support <- make_identity_volume(support_order = 4:1)
  transform_matrix <- diag(c(2, 0.5), nrow = 2)
  changed_transform <- make_identity_volume(
    analysis_transform = list(
      type = "linear",
      dim = 2L,
      matrix = transform_matrix,
      to_analysis = function(x) transform_matrix %*% as.matrix(x),
      to_raw = function(x) solve(transform_matrix, as.matrix(x))
    )
  )

  base_id <- latent_space_id(baseline)
  component_id <- latent_space_id(reordered_components)
  domain_id <- latent_space_id(changed_domain)
  support_id <- latent_space_id(changed_support)
  transform_id <- latent_space_id(changed_transform)

  expect_false(identical(base_id$coordinate_id, component_id$coordinate_id))
  expect_identical(
    base_id$decoder_domain_id,
    component_id$decoder_domain_id
  )
  expect_false(identical(base_id$coordinate_id, transform_id$coordinate_id))
  expect_identical(base_id$decoder_domain_id, transform_id$decoder_domain_id)

  expect_identical(base_id$coordinate_id, domain_id$coordinate_id)
  expect_false(identical(base_id$decoder_domain_id, domain_id$decoder_domain_id))
  expect_identical(base_id$coordinate_id, support_id$coordinate_id)
  expect_false(identical(base_id$decoder_domain_id, support_id$decoder_domain_id))
})

test_that("analysis and raw coordinates have explicit distinct identities", {
  latent <- make_identity_volume()

  analysis_id <- latent_space_id(latent, coordinates = "analysis")
  raw_id <- latent_space_id(latent, coordinates = "raw")

  expect_false(identical(analysis_id$coordinate_id, raw_id$coordinate_id))
  expect_identical(
    analysis_id$decoder_domain_id,
    raw_id$decoder_domain_id
  )
  expect_false(identical(analysis_id$decoder_id, raw_id$decoder_id))
})

test_that("identity survives serialization and fresh observation scores", {
  original <- make_identity_volume()
  recreated <- unserialize(serialize(original, NULL, version = 3L))
  recreated@basis[] <- recreated@basis[] * -10

  expect_identical(latent_space_id(original), latent_space_id(recreated))
})

test_that("handle-backed identity remains lazy", {
  kind <- "test_latent_space_id_no_materialize"
  register_handle_kind(kind, function(handle) {
    stop("identity materialization forbidden", call. = FALSE)
  }, type = "loadings")
  handle <- new(
    "LoadingsHandle",
    id = "latent-space-id-no-materialize",
    dim = as.integer(c(4L, 2L)),
    kind = kind,
    spec = list(derivation = "throwing-identity-receipt"),
    label = "throwing identity handle"
  )
  latent <- make_identity_volume()
  latent@loadings <- handle

  expect_no_error(identity <- latent_space_id(latent))
  expect_match(identity$coordinate_id, "^latent-coordinate:")
  expect_match(identity$decoder_domain_id, "^latent-decoder-domain:")
  expect_error(
    decoder(latent)$forward(c(1, 0)),
    "identity materialization forbidden"
  )
})

test_that("surface identity separates geometry and support from coordinates", {
  skip_if_not_installed("neurosurf")
  baseline <- make_identity_surface("left")
  reordered <- make_identity_surface("left", support_order = 3:1)

  base_id <- latent_space_id(baseline)
  reordered_id <- latent_space_id(reordered)
  expect_identical(base_id$coordinate_id, reordered_id$coordinate_id)
  expect_false(identical(
    base_id$decoder_domain_id,
    reordered_id$decoder_domain_id
  ))
})

test_that("bilateral identity derives from ordered hemisphere identities", {
  skip_if_not_installed("neurosurf")
  left <- make_identity_surface("left")
  right <- make_identity_surface("right")
  baseline <- BilatLatentNeuroSurfaceVector(left, right)
  recreated <- unserialize(serialize(baseline, NULL, version = 3L))
  reordered_right <- make_identity_surface("right", support_order = 3:1)
  changed <- BilatLatentNeuroSurfaceVector(left, reordered_right)

  base_id <- latent_space_id(baseline)
  changed_id <- latent_space_id(changed)

  expect_named(base_id$child_ids, c("left", "right"))
  expect_identical(base_id, latent_space_id(recreated))
  expect_identical(base_id$coordinate_id, changed_id$coordinate_id)
  expect_false(identical(
    base_id$decoder_domain_id,
    changed_id$decoder_domain_id
  ))
})

test_that("composite identity derives from named ordered children", {
  left <- make_identity_volume()
  right <- make_identity_volume(origin = c(2, 0, 0))
  baseline <- BlockLatentNeuroVector(list(left = left, right = right))
  recreated <- unserialize(serialize(baseline, NULL, version = 3L))
  reordered <- BlockLatentNeuroVector(list(right = right, left = left))

  base_id <- latent_space_id(baseline)
  recreated_id <- latent_space_id(recreated)
  reordered_id <- latent_space_id(reordered)

  expect_named(base_id$child_ids, c("left", "right"))
  expect_identical(base_id, recreated_id)
  expect_false(identical(base_id$coordinate_id, reordered_id$coordinate_id))
  expect_false(identical(
    base_id$decoder_domain_id,
    reordered_id$decoder_domain_id
  ))
})

test_that("transport identity excludes observations and separates decoder spaces", {
  template <- make_identity_transport_template()
  operator <- make_identity_transport_operator()
  decoder_matrix <- .materialize_linear_map(operator) %*%
    .materialize_linear_map(basis_decoder(template))
  anchor_scores <- matrix(c(1, 0, 0, 1, 1, 1), nrow = 3, byrow = TRUE)
  task_scores <- matrix(c(4, -2, 3, 5), nrow = 2, byrow = TRUE)
  anchor <- encode_transport(
    anchor_scores %*% t(decoder_matrix),
    basis_asset = template,
    field_operator = operator,
    center = FALSE,
    lambda = 0
  )
  task <- encode_transport(
    task_scores %*% t(decoder_matrix),
    basis_asset = template,
    field_operator = operator,
    center = FALSE,
    lambda = 0
  )

  anchor_native <- latent_space_id(anchor, space = "native")
  task_native <- latent_space_id(task, space = "native")
  anchor_template <- latent_space_id(anchor, space = "template")

  expect_identical(anchor_native, task_native)
  expect_identical(
    anchor_native$coordinate_id,
    anchor_template$coordinate_id
  )
  expect_false(identical(
    anchor_native$decoder_domain_id,
    anchor_template$decoder_domain_id
  ))
  expect_false(identical(anchor_native$decoder_id, anchor_template$decoder_id))
})

test_that("plain implicit latents fail the supported identity contract", {
  latent <- implicit_latent(
    coeff = matrix(1:4, nrow = 2),
    decoder = function(time_idx = NULL, roi_mask = NULL,
                       levels_keep = NULL) matrix(1, nrow = 2, ncol = 2),
    meta = list(family = "plain_implicit"),
    mask = array(TRUE, dim = c(2, 1, 1))
  )

  expect_error(
    latent_space_id(latent),
    class = "fmrilatent_error_unsupported_operation"
  )
})
