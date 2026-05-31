# Temporal encode_spec methods

#' @include encode.R
NULL

#' @exportS3Method
encode_spec.spec_time_slepian <- function(x, spec, mask, reduction, materialize, label, ...) {
  materialize <- .resolve_materialize(materialize, context = "encode_spec.spec_time_slepian")
  n_time <- nrow(x)
  k_requested <- spec$k %||% .slepian_default_k(n_time, spec$tr, spec$bandwidth)
  k_use <- max(1L, min(as.integer(k_requested), n_time))
  if (materialize == "matrix") {
    basis <- dpss_time_basis(n_time, tr = spec$tr, bandwidth = spec$bandwidth, k = k_use, backend = spec$backend)
  } else {
    basis <- slepian_temporal_handle(n_time = n_time, tr = spec$tr,
      bandwidth = spec$bandwidth, k = k_use, backend = spec$backend)
  }
  basis_matrix <- basis_mat(basis)
  loadings <- .temporal_loadings_from_basis(
    x,
    basis_matrix,
    context = "encode_spec.spec_time_slepian"
  )
  spc <- .space_with_time_from_mask(mask, n_time, "encode_spec.spec_time_slepian")
  meta <- list(
    family = "time_slepian",
    k = ncol(basis_matrix),
    tr = spec$tr,
    bandwidth = spec$bandwidth,
    backend = spec$backend
  )
  LatentNeuroVec(basis = basis, loadings = loadings, space = spc, mask = mask,
                 offset = numeric(0), label = label, meta = meta)
}

#' @exportS3Method
encode_spec.spec_time_dct <- function(x, spec, mask, reduction, materialize, label, ...) {
  materialize <- .resolve_materialize(materialize, context = "encode_spec.spec_time_dct")
  n_time <- nrow(x)
  k_use <- .validate_auto_positive_count(spec$k, n_time, default = n_time,
                                         name = "spec_time_dct$k",
                                         n_name = "n_time")
  if (materialize == "matrix") {
    basis <- build_dct_basis(n_time, k = k_use, norm = spec$norm)
  } else {
    basis <- dct_basis_handle(n_time = n_time, k = k_use, norm = spec$norm)
  }
  loadings <- .temporal_loadings_from_basis(
    x,
    basis_mat(basis),
    context = "encode_spec.spec_time_dct"
  )
  spc <- .space_with_time_from_mask(mask, n_time, "encode_spec.spec_time_dct")
  meta <- list(
    family = "time_dct",
    k = k_use,
    norm = spec$norm
  )
  LatentNeuroVec(basis = basis, loadings = loadings, space = spc, mask = mask,
                 offset = numeric(0), label = label, meta = meta)
}

#' @exportS3Method
encode_spec.spec_time_bspline <- function(x, spec, mask, reduction, materialize, label, ...) {
  materialize <- .resolve_materialize(materialize, context = "encode_spec.spec_time_bspline")
  n_time <- nrow(x)
  k_use <- .validate_auto_positive_count(spec$k, n_time, default = min(5L, n_time),
                                         name = "spec_time_bspline$k",
                                         n_name = "n_time")
  if (materialize == "matrix") {
    basis <- build_bspline_basis(
      n_time = n_time,
      k = k_use,
      degree = spec$degree,
      include_intercept = spec$include_intercept,
      orthonormalize = spec$orthonormalize
    )
  } else {
    basis <- bspline_basis_handle(
      n_time = n_time,
      k = k_use,
      degree = spec$degree,
      include_intercept = spec$include_intercept,
      orthonormalize = spec$orthonormalize,
      id = NULL,
      label = NULL
    )
  }
  loadings <- .temporal_loadings_from_basis(
    x,
    basis_mat(basis),
    context = "encode_spec.spec_time_bspline"
  )
  spc <- .space_with_time_from_mask(mask, n_time, "encode_spec.spec_time_bspline")
  meta <- list(
    family = "time_bspline",
    k = k_use,
    degree = spec$degree,
    include_intercept = spec$include_intercept,
    orthonormalize = spec$orthonormalize
  )
  LatentNeuroVec(basis = basis, loadings = loadings, space = spc, mask = mask,
                 offset = numeric(0), label = label, meta = meta)
}
