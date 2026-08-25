make_declared_units <- function(response_scaling = "percent_signal_change",
                                sign_convention = "arbitrary_component_sign") {
  latent_units_record(
    response_scaling = response_scaling,
    coefficient_units = "percent_signal_change",
    loading_normalization = "unit_l2_columns",
    loading_metric = "euclidean_sample_inner_product",
    analysis_coordinate_metric = "euclidean",
    sign_convention = sign_convention
  )
}

make_units_mask <- function() {
  neuroim2::LogicalNeuroVol(
    array(TRUE, dim = c(2, 2, 1)),
    neuroim2::NeuroSpace(c(2, 2, 1))
  )
}

make_units_volume <- function(units = NULL, origin = c(0, 0, 0), n_time = 3L) {
  mask <- neuroim2::LogicalNeuroVol(
    array(TRUE, dim = c(2, 2, 1)),
    neuroim2::NeuroSpace(c(2, 2, 1), origin = origin)
  )
  meta <- list(family = "units_fixture")
  if (!is.null(units)) meta$latent_units <- units
  LatentNeuroVec(
    basis = Matrix::Matrix(matrix(seq_len(n_time * 2L), nrow = n_time, ncol = 2)),
    loadings = Matrix::Matrix(matrix(c(
      1, 0,
      0, 1,
      1, 1,
      0, 1
    ), nrow = 4, byrow = TRUE)),
    space = neuroim2::NeuroSpace(c(2, 2, 1, n_time), origin = origin),
    mask = mask,
    meta = meta
  )
}

make_units_transport_template <- function(mask) {
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
      meta = list(family = "units_transport_template", k = 2L)
    ),
    class = "ParcelBasisTemplate"
  )
}

test_that("declared units records are deterministic and immutable", {
  units <- make_declared_units()
  recreated <- make_declared_units()
  annotated <- latent_units_record(
    response_scaling = "percent_signal_change",
    coefficient_units = "percent_signal_change",
    loading_normalization = "unit_l2_columns",
    loading_metric = "euclidean_sample_inner_product",
    analysis_coordinate_metric = "euclidean",
    sign_convention = "arbitrary_component_sign",
    notes = "reader-facing note"
  )

  expect_s3_class(units, "fmrilatent_units")
  expect_identical(units, recreated)
  expect_identical(units$status, "declared")
  expect_match(units$compatibility_id, "^latent-units:")
  expect_identical(units$compatibility_id, annotated$compatibility_id)
  expect_false(identical(units$integrity_id, annotated$integrity_id))
  expect_error(units$response_scaling <- "raw_signal", "immutable")
  expect_error(units[["response_scaling"]] <- "raw_signal", "immutable")
  expect_error(units["response_scaling"] <- "raw_signal", "immutable")
  expect_error(names(units) <- rev(names(units)), "immutable")

  tampered <- unclass(units)
  tampered$response_scaling <- "raw_signal"
  class(tampered) <- c("fmrilatent_units", "list")
  latent <- make_units_volume()
  latent@meta$latent_units <- tampered
  expect_error(latent_units(latent), "integrity check")
})

test_that("legacy latent objects report undeclared units without guessing", {
  explicit <- make_units_volume()
  implicit <- implicit_latent(
    coeff = matrix(1:4, nrow = 2),
    decoder = function(time_idx = NULL, roi_mask = NULL,
                       levels_keep = NULL) matrix(1, 2, 2),
    meta = list(family = "legacy_implicit"),
    mask = array(TRUE, dim = c(2, 1, 1))
  )

  explicit_units <- latent_units(explicit)
  implicit_units <- latent_units(implicit)
  expect_identical(explicit_units$status, "undeclared")
  expect_identical(implicit_units$status, "undeclared")
  expect_identical(
    explicit_units$compatibility_id,
    "latent-units:undeclared"
  )
  expect_true(all(is.na(unlist(
    explicit_units[c(
      "response_scaling", "coefficient_units", "loading_normalization",
      "loading_metric", "analysis_coordinate_metric", "sign_convention"
    )]
  ))))
})

test_that("standard encode captures caller-declared units", {
  units <- make_declared_units()
  mask <- make_units_mask()
  values <- matrix(seq_len(20), nrow = 5, ncol = 4)

  encoded <- encode(
    values,
    spec_time_dct(k = 2),
    mask = mask,
    materialize = "matrix",
    units = units
  )
  legacy <- encode(
    values,
    spec_time_dct(k = 2),
    mask = mask,
    materialize = "matrix"
  )

  expect_identical(latent_units(encoded), units)
  expect_identical(latent_meta(encoded)$latent_units, units)
  expect_identical(latent_units(legacy)$status, "undeclared")
  expect_error(
    encode(
      values,
      spec_time_dct(k = 2),
      mask = mask,
      materialize = "matrix",
      units = list(response_scaling = "guessed")
    ),
    "latent_units_record"
  )
})

test_that("transport encode captures the same units contract", {
  units <- make_declared_units()
  mask <- make_units_mask()
  template <- make_units_transport_template(mask)
  field_operator <- as_portable_linear_map(
    diag(4),
    source_domain_id = "units_template_field",
    target_domain_id = "units_native_field",
    target_support = mask,
    provenance = list(operator = "units_identity")
  )
  decoder_matrix <- field_operator$materialize() %*%
    basis_decoder(template)$materialize()
  observations <- matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE) %*%
    t(decoder_matrix)

  encoded <- encode_transport(
    observations,
    basis_asset = template,
    field_operator = field_operator,
    center = FALSE,
    lambda = 0,
    units = units
  )

  expect_identical(latent_units(encoded), units)
  expect_identical(latent_meta(encoded)$latent_units, units)
  expect_true("units" %in% names(formals(encode_operator)))
  expect_true("units" %in% names(formals(encode_transport)))
  expect_true("units" %in% names(formals(encode_awpt)))
})

test_that("composite units derive only from compatible children", {
  units <- make_declared_units()
  other_units <- make_declared_units(response_scaling = "raw_signal")
  compatible <- BlockLatentNeuroVector(list(
    first = make_units_volume(units),
    second = make_units_volume(units, origin = c(3, 0, 0))
  ))
  incompatible <- BlockLatentNeuroVector(list(
    first = make_units_volume(units),
    second = make_units_volume(other_units, origin = c(3, 0, 0))
  ))
  legacy <- BlockLatentNeuroVector(list(
    first = make_units_volume(),
    second = make_units_volume(origin = c(3, 0, 0))
  ))

  expect_identical(latent_units(compatible), units)
  expect_identical(latent_units(incompatible)$status, "incompatible")
  expect_false(identical(
    latent_units(incompatible)$child_compatibility_ids[[1L]],
    latent_units(incompatible)$child_compatibility_ids[[2L]]
  ))
  expect_identical(latent_units(legacy)$status, "undeclared")
})

test_that("units compatibility is observation-independent and serializable", {
  units <- make_declared_units()
  first <- make_units_volume(units)
  second <- make_units_volume(units, n_time = 7L)
  second@basis[] <- rev(seq_along(second@basis))
  recreated <- unserialize(serialize(first, NULL, version = 3L))

  expect_identical(latent_units(first), latent_units(second))
  expect_identical(latent_units(first), latent_units(recreated))
})
