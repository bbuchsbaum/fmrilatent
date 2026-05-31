# Spatial encode_spec methods

#' @include encode.R
NULL

#' @exportS3Method
encode_spec.spec_space_slepian <- function(x, spec, mask, reduction, materialize, label, ...) {
  materialize <- .resolve_materialize(materialize, context = "encode_spec.spec_space_slepian")
  mask_arr <- .extract_mask_array(mask, "encode_spec.spec_space_slepian")
  if (is.null(reduction)) reduction <- make_cluster_reduction(mask, seq_len(sum(mask_arr)))
  k_use <- spec$k %||% 3L
  k_use <- .validate_positive_count(k_use, "spec_space_slepian$k")
  if (materialize == "matrix") {
    loadings <- lift(reduction, basis_slepian(k = k_use), k_neighbors = spec$k_neighbors)
  } else {
    loadings <- slepian_spatial_loadings_handle(
      reduction,
      basis_spec = basis_slepian(k = k_use),
      data = NULL,
      k_neighbors = spec$k_neighbors,
      label = "slepian-spatial"
    )
  }
  # Multiply in Matrix space so a sparse loadings dictionary is never densified;
  # only the small time x k product is realized as a dense matrix.
  basis <- as.matrix(x %*% loadings_mat(loadings))
  spc <- .space_with_time_from_mask(mask, nrow(x), "encode_spec.spec_space_slepian")
  meta <- list(
    family = "space_slepian",
    k = k_use,
    k_neighbors = spec$k_neighbors
  )
  LatentNeuroVec(basis = basis, loadings = loadings, space = spc, mask = mask,
                 offset = numeric(0), label = label, meta = meta)
}

#' @exportS3Method
encode_spec.spec_space_heat <- function(x, spec, mask, reduction, materialize, label, ...) {
  materialize <- .resolve_materialize(materialize, context = "encode_spec.spec_space_heat")
  mask_arr <- .extract_mask_array(mask, "encode_spec.spec_space_heat")
  if (is.null(reduction)) reduction <- make_cluster_reduction(mask, seq_len(sum(mask_arr)))
  spec_hw <- basis_heat_wavelet(scales = spec$scales, order = spec$order, threshold = spec$threshold)
  if (materialize == "matrix") {
    loadings <- lift(reduction, spec_hw, k_neighbors = spec$k_neighbors)
  } else {
    loadings <- heat_wavelet_loadings_handle(
      reduction,
      basis_spec = spec_hw,
      data = NULL,
      k_neighbors = spec$k_neighbors,
      label = "heat-wavelet"
    )
  }
  # Multiply in Matrix space so a sparse loadings dictionary is never densified;
  # only the small time x k product is realized as a dense matrix.
  basis <- as.matrix(x %*% loadings_mat(loadings))
  spc <- .space_with_time_from_mask(mask, nrow(x), "encode_spec.spec_space_heat")
  meta <- list(
    family = "space_heat",
    scales = spec$scales,
    order = spec$order,
    sparsify_eps = spec$sparsify_eps %||% spec$threshold,
    threshold = spec$threshold,
    k_neighbors = spec$k_neighbors
  )
  LatentNeuroVec(basis = basis, loadings = loadings, space = spc, mask = mask,
                 offset = numeric(0), label = label, meta = meta)
}

#' @exportS3Method
encode_spec.spec_space_hrbf <- function(x, spec, mask, reduction, materialize, label, ...) {
  materialize <- .resolve_materialize(materialize, context = "encode_spec.spec_space_hrbf")
  .extract_mask_array(mask, "encode_spec.spec_space_hrbf")
  params <- spec$params %||% list()
  B_atoms <- hrbf_generate_basis(params, mask) # atoms x vox
  loadings_mat_hrbf <- Matrix::t(Matrix::Matrix(B_atoms, sparse = TRUE)) # vox x atoms
  coeff <- as.matrix(Matrix::tcrossprod(as.matrix(x), B_atoms))
  basis <- Matrix::Matrix(coeff, sparse = FALSE) # time x atoms
  loadings <- .loadings_for_materialize(
    loadings_mat_hrbf,
    materialize,
    id_prefix = "hrbf-loadings",
    label = "hrbf-loadings"
  )
  spc <- .space_with_time_from_mask(mask, nrow(x), "encode_spec.spec_space_hrbf")
  meta <- list(
    family = "space_hrbf",
    params = params
  )
  LatentNeuroVec(basis = basis, loadings = loadings, space = spc, mask = mask,
                 offset = numeric(0), label = label, meta = meta)
}

#' @exportS3Method
encode_spec.spec_space_pca <- function(x, spec, mask, reduction, materialize, label, ...) {
  materialize <- .resolve_materialize(materialize, context = "encode_spec.spec_space_pca")
  mask_arr <- .extract_mask_array(mask, "encode_spec.spec_space_pca")
  n_time <- nrow(x)
  n_vox <- sum(mask_arr)

  if (is.null(reduction)) {
    reduction <- make_cluster_reduction(mask, rep.int(1L, n_vox))
  }
  k_use <- spec$k %||% 3L
  k_use <- .validate_positive_count(k_use, "spec_space_pca$k")

  offset <- .encode_center(x, center = isTRUE(spec$center),
                           context = "encode_spec.spec_space_pca")$offset

  loadings <- lift(
    reduction,
    basis_pca(k = k_use, whiten = FALSE),
    data = x,
    center = isTRUE(spec$center),
    scale = isTRUE(spec$scale),
    offset = if (length(offset) > 0) offset else NULL,
    backend = spec$backend %||% "auto",
    ...
  )

  basis <- x %*% loadings
  # Centering is owned by lift.basis_pca: it returns the per-voxel mean projected
  # onto the PCA loadings as the `fmrilatent.mean_scores` attribute. We subtract
  # that mean contribution here rather than re-deriving it from `offset`, so the
  # mean-removal lives in a single place.
  mean_scores <- attr(loadings, "fmrilatent.mean_scores")
  if (length(mean_scores) > 0) {
    basis <- basis - matrix(1, nrow = n_time, ncol = 1) %*%
      matrix(mean_scores, nrow = 1)
  }

  if (isTRUE(spec$whiten)) {
    d <- attr(loadings, "fmrilatent.singular_values")
    if (is.null(d) || length(d) != ncol(loadings)) {
      .encoder_cli_abort(
        "PCA whitening requested, but singular values were not returned by lift().",
        class = "fmrilatent_error_invalid_pca_whitening"
      )
    }
    if (any(!is.finite(d)) || any(d <= 0)) {
      .encoder_cli_abort(
        "PCA whitening requires strictly positive finite singular values.",
        class = "fmrilatent_error_invalid_pca_whitening"
      )
    }
    scale_fac <- sqrt(max(1, n_time - 1))
    basis <- sweep(as.matrix(basis), 2, d, "/") * scale_fac
    loadings <- loadings %*% Matrix::Diagonal(x = d / scale_fac)
  }
  loadings <- .loadings_for_materialize(
    loadings,
    materialize,
    id_prefix = "pca-loadings",
    label = "pca-loadings"
  )

  spc <- .space_with_time_from_mask(mask, n_time, "encode_spec.spec_space_pca")
  meta <- list(
    family = "pca_spatial",
    k = k_use,
    center = isTRUE(spec$center),
    scale = isTRUE(spec$scale),
    whiten = isTRUE(spec$whiten),
    backend = spec$backend %||% "auto"
  )

  LatentNeuroVec(
    basis = Matrix::Matrix(basis, sparse = FALSE),
    loadings = loadings,
    space = spc,
    mask = mask,
    offset = offset,
    label = label,
    meta = meta
  )
}

#' @exportS3Method
encode_spec.spec_space_wavelet_active <- function(x, spec, mask, reduction, materialize, label, ...) {
  materialize <- .resolve_materialize(
    materialize,
    supported = "matrix",
    default = "matrix",
    context = "encode_spec.spec_space_wavelet_active"
  )
  result <- wavelet_active_latent(
    X = x,
    mask = mask,
    levels_space = spec$levels_space,
    levels_time = spec$levels_time,
    threshold = spec$threshold
  )
  if (!is(result, "LatentNeuroVec") && !inherits(result, "ImplicitLatent")) {
    .encoder_cli_abort(
      paste0(
        "encode_spec.spec_space_wavelet_active expected wavelet_active_latent() ",
        "to return a LatentNeuroVec or ImplicitLatent, got: ",
        paste(class(result), collapse = ", ")
      ),
      class = "fmrilatent_error_invalid_type"
    )
  }
  result
}
