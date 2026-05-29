# Spatial encode_spec methods

#' @include encode.R
NULL

#' @exportS3Method
encode_spec.spec_space_slepian <- function(x, spec, mask, reduction, materialize, label, ...) {
  materialize <- .resolve_materialize(materialize, context = "encode_spec.spec_space_slepian")
  mask_arr <- .extract_mask_array(mask, "encode_spec.spec_space_slepian")
  if (is.null(reduction)) reduction <- make_cluster_reduction(mask, seq_len(sum(mask_arr)))
  if (materialize == "matrix") {
    loadings <- lift(reduction, basis_slepian(k = spec$k), k_neighbors = spec$k_neighbors)
  } else {
    loadings <- slepian_spatial_loadings_handle(
      reduction,
      basis_spec = basis_slepian(k = spec$k),
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
    k = spec$k,
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
    threshold = spec$threshold,
    k_neighbors = spec$k_neighbors
  )
  LatentNeuroVec(basis = basis, loadings = loadings, space = spc, mask = mask,
                 offset = numeric(0), label = label, meta = meta)
}

#' @exportS3Method
encode_spec.spec_space_hrbf <- function(x, spec, mask, reduction, materialize, label, ...) {
  materialize <- .resolve_materialize(materialize, context = "encode_spec.spec_space_hrbf")
  mask_arr <- .extract_mask_array(mask, "encode_spec.spec_space_hrbf")
  params <- spec$params %||% list()
  B_atoms <- hrbf_generate_basis(params, mask) # atoms x vox
  loadings_mat_hrbf <- Matrix::t(Matrix::Matrix(B_atoms, sparse = TRUE)) # vox x atoms
  gram <- as.matrix(B_atoms %*% Matrix::t(B_atoms))
  rhs <- t(as.matrix(x) %*% Matrix::t(B_atoms))
  coeff <- t(.robust_gram_solve(gram, rhs))
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

  offset <- numeric(0)
  if (isTRUE(spec$center)) {
    offset <- colMeans(x)
  }

  loadings <- lift(
    reduction,
    basis_pca(k = spec$k, whiten = isTRUE(spec$whiten)),
    data = x,
    center = isTRUE(spec$center),
    offset = if (length(offset) > 0) offset else NULL,
    backend = spec$backend %||% "auto",
    ...
  )

  basis <- x %*% loadings
  if (length(offset) > 0) {
    mu_scores <- as.matrix(crossprod(offset, loadings))
    basis <- basis - matrix(1, nrow = n_time, ncol = 1) %*% mu_scores
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
    k = spec$k,
    center = isTRUE(spec$center),
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
    supported = c("matrix", "handle"),
    default = "matrix",
    context = "encode_spec.spec_space_wavelet_active"
  )
  mask_arr <- .extract_mask_array(mask, "encode_spec.spec_space_wavelet_active")
  wavelet_active_latent(
    X = x,
    mask = mask,
    levels_space = spec$levels_space,
    levels_time = spec$levels_time,
    threshold = spec$threshold
  )
}
