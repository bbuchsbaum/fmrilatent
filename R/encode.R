# High-level encoding API and specs

# Robust Gram matrix solve with ridge regularization
# Adds a small ridge to the diagonal if solve fails
.robust_gram_solve <- function(gram, rhs, ridge = 1e-8) {
  gram_mat <- as.matrix(gram)
  result <- tryCatch(
    solve(gram_mat, rhs),
    error = function(e) NULL
  )
  if (is.null(result)) {
    # Add ridge regularization and retry
    diag(gram_mat) <- diag(gram_mat) + ridge
    result <- tryCatch(
      solve(gram_mat, rhs),
      error = function(e) {
        # Final fallback: pseudoinverse via SVD
        svd_g <- svd(gram_mat)
        tol <- max(dim(gram_mat)) * max(svd_g$d) * .Machine$double.eps
        pos <- svd_g$d > tol
        if (any(pos)) {
          svd_g$v[, pos, drop = FALSE] %*% (t(svd_g$u[, pos, drop = FALSE]) / svd_g$d[pos]) %*% rhs
        } else {
          matrix(0, nrow = ncol(gram_mat), ncol = ncol(rhs))
        }
      }
    )
  }
  result
}

# --- Spec constructors --------------------------------------------------------

#' Temporal Slepian/DPSS spec
#'
#' @param tr Repetition time (seconds).
#' @param bandwidth Half-bandwidth in Hz (default 0.1).
#' @param k Optional number of components (default floor(2*NW)-1).
#' @param backend Backend to use ("tridiag" or "dense").
#' @return A `spec_time_slepian` object for `encode()` / `spec_st()`.
#' @export
spec_time_slepian <- function(tr, bandwidth = 0.1, k = NULL, backend = c("tridiag", "dense")) {
  backend <- match.arg(backend)
  structure(list(tr = tr, bandwidth = bandwidth, k = k, backend = backend), class = "spec_time_slepian")
}

#' Temporal DCT spec
#'
#' @param k Components.
#' @param norm Normalization ("ortho" or "none").
#' @return A `spec_time_dct` object.
#' @export
spec_time_dct <- function(k, norm = c("ortho", "none")) {
  norm <- match.arg(norm)
  structure(list(k = k, norm = norm), class = "spec_time_dct")
}

#' Temporal B-spline spec
#'
#' @param k Components (df).
#' @param degree Spline degree (default 3).
#' @param include_intercept Logical include intercept.
#' @param orthonormalize Logical orthonormalize columns (default TRUE).
#' @return A `spec_time_bspline` object.
#' @export
spec_time_bspline <- function(k, degree = 3L, include_intercept = FALSE, orthonormalize = TRUE) {
  structure(list(k = k, degree = degree, include_intercept = include_intercept,
                 orthonormalize = orthonormalize),
            class = "spec_time_bspline")
}

#' Spatial Slepian spec
#'
#' @param k Components per cluster.
#' @param k_neighbors k-NN graph parameter.
#' @return A `spec_space_slepian` object.
#' @export
spec_space_slepian <- function(k = 3L, k_neighbors = 6L) {
  structure(list(k = k, k_neighbors = k_neighbors), class = "spec_space_slepian")
}

#' Spatial PCA spec (cluster-local)
#'
#' Computes PCA eigenvectors within each cluster/parcel specified by a
#' `ClusterReduction` and returns a block-sparse spatial dictionary.
#'
#' @param k Components per cluster.
#' @param center Logical; center voxels before PCA (default TRUE). When TRUE,
#'   voxel means are stored in `LatentNeuroVec@offset`.
#' @param whiten Logical; if TRUE, return whitened scores (unit-variance) and
#'   rescaled loadings such that reconstruction is unchanged.
#' @param backend SVD backend: "auto" (default), "svds" (RSpectra), or "svd" (base).
#' @return A `spec_space_pca` object.
#' @export
spec_space_pca <- function(k = 3L, center = TRUE, whiten = FALSE,
                           backend = c("auto", "svds", "svd")) {
  backend <- match.arg(backend)
  structure(
    list(
      k = as.integer(k),
      center = isTRUE(center),
      whiten = isTRUE(whiten),
      backend = backend
    ),
    class = "spec_space_pca"
  )
}

#' Spatial heat-wavelet spec (graph diffusion)
#'
#' @param scales Heat scales.
#' @param order Polynomial order.
#' @param threshold Threshold for small coefficients.
#' @param k_neighbors k-NN graph parameter.
#' @return A `spec_space_heat` object.
#' @export
spec_space_heat <- function(scales = c(1, 2, 4, 8), order = 30L, threshold = 1e-6, k_neighbors = 6L) {
  structure(list(scales = scales, order = order, threshold = threshold, k_neighbors = k_neighbors),
            class = "spec_space_heat")
}

#' Spatial HRBF spec
#'
#' @param params HRBF parameter list (sigma0, levels, radius_factor, kernel_type, seed).
#' @return A `spec_space_hrbf` object.
#' @export
spec_space_hrbf <- function(params = list()) {
  structure(list(params = params), class = "spec_space_hrbf")
}

#' Spatial wavelet (active pencil) spec
#'
#' @param levels_space Spatial lifting levels.
#' @param levels_time Optional time lifting levels.
#' @param threshold Threshold after space transform.
#' @return A `spec_space_wavelet_active` object.
#' @export
spec_space_wavelet_active <- function(levels_space = 2L, levels_time = 0L, threshold = 0) {
  structure(list(levels_space = levels_space, levels_time = levels_time, threshold = threshold),
            class = "spec_space_wavelet_active")
}
#' Spatiotemporal spec (separable)
#'
#' @param time Temporal spec (`spec_time_*`).
#' @param space Spatial spec (`spec_space_*`).
#' @param core_mode Reserved (currently "auto").
#' @return A `spec_st` object.
#' @export
spec_st <- function(time, space, core_mode = c("auto", "explicit")) {
  core_mode <- match.arg(core_mode)
  structure(list(time = time, space = space, core_mode = core_mode), class = "spec_st")
}

# --- Encode generic -----------------------------------------------------------

#' Encode data into a latent representation using a spec
#'
#' @param x Matrix (time x voxels in mask order).
#' @param spec Spec object created by `spec_time_*`, `spec_space_*`, or `spec_st`.
#' @param mask LogicalNeuroVol or logical array (required for spatial pieces).
#' @param reduction Optional GraphReduction (for spatial specs).
#' @param materialize "handle", "matrix", or "auto" (default "handle").
#' @param label Optional label.
#' @param ... Additional arguments passed to methods.
#' @return A `LatentNeuroVec` (explicit bases) or `ImplicitLatent` (separable cases).
#' @export
encode <- function(x, spec, mask, reduction = NULL,
                    materialize = c("handle", "auto", "matrix"),
                    label = "", ...) {
  UseMethod("encode")
}

#' @export
encode.default <- function(x, spec, mask, reduction = NULL,
                           materialize = c("handle", "auto", "matrix"),
                           label = "", ...) {
  stop("No encode method for class: ", paste(class(x), collapse = ","), call. = FALSE)
}

#' @export
encode.matrix <- function(x, spec, mask, reduction = NULL,
                          materialize = c("handle", "auto", "matrix"),
                          label = "", ...) {
  materialize <- match.arg(materialize)
  encode_spec(
    x, spec,
    mask = mask,
    reduction = reduction,
    materialize = materialize,
    label = label,
    ...
  )
}

#' @export
encode.NeuroVec <- function(x, spec, mask, reduction = NULL,
                            materialize = c("handle", "auto", "matrix"),
                            label = "", ...) {
  materialize <- match.arg(materialize)
  X <- t(neuroim2::series(x, mask != 0))  # series returns voxels x time, transpose to time x voxels
  encode(X, spec, mask = mask, reduction = reduction,
         materialize = materialize, label = label, ...)
}

# --- Factory helper -----------------------------------------------------------

#' Simple factory to build a spec and encode in one call
#'
#' @param family One of: "dct_time", "slepian_time", "slepian_space", "heat_space", "slepian_st".
#' @param x Data matrix (time x voxels).
#' @param mask Mask (required for spatial families).
#' @param reduction Optional GraphReduction for spatial specs.
#' @param ... Passed to spec constructors and encode().
#' @param materialize "handle", "matrix", or "auto" (default "handle").
#' @param label Optional label for the resulting object.
#' @return A LatentNeuroVec or ImplicitLatent object.
#' @export
latent_factory <- function(family, x, mask, reduction = NULL, ..., materialize = "handle", label = "") {
  family <- match.arg(family, c("dct_time", "slepian_time", "slepian_space", "pca_space", "parcel_space", "heat_space", "slepian_st", "bspline_hrbf_st", "wavelet_active", "hierarchical"))
  spec <- switch(
    family,
    dct_time = spec_time_dct(...),
    slepian_time = spec_time_slepian(...),
    slepian_space = spec_space_slepian(...),
    pca_space = spec_space_pca(...),
    heat_space = spec_space_heat(...),
    bspline_hrbf_st = {
      args <- list(...)
      time_spec <- args$time %||% spec_time_bspline(k = args$k_time %||% args$k %||% 5L,
                                                   degree = args$degree %||% 3L,
                                                   include_intercept = args$include_intercept %||% FALSE,
                                                   orthonormalize = TRUE)
      space_spec <- args$space %||% spec_space_hrbf(params = args$params %||% list())
      spec_st(time = time_spec, space = space_spec)
    },
    parcel_space = spec_space_parcel(...),
    wavelet_active = spec_space_wavelet_active(...),
    hierarchical = spec_hierarchical_template(...),
    slepian_st = {
      args <- list(...)
      time_spec <- args$time %||% spec_time_slepian(tr = args$tr, bandwidth = args$bandwidth %||% 0.1, k = args$k_time %||% NULL)
      space_spec <- args$space %||% spec_space_slepian(k = args$k_space %||% 3L, k_neighbors = args$k_neighbors %||% 6L)
      spec_st(time = time_spec, space = space_spec)
    }
  )
  encode(x, spec, mask = mask, reduction = reduction, materialize = materialize, label = label)
}
#' Dispatch encoding based on spec type
#'
#' @param x Data matrix.
#' @param spec Spec object.
#' @param ... Additional arguments passed to methods.
#' @return Encoded representation.
#' @export
encode_spec <- function(x, spec, ...) UseMethod("encode_spec", spec)

#' @exportS3Method
encode_spec.default <- function(x, spec, ...) {
  stop("Unknown spec class: ", paste(class(spec), collapse = ","), call. = FALSE)
}

# Hierarchical template spec -------------------------------------------------

#' Create hierarchical template spec
#' @param template HierarchicalBasisTemplate object
#' @param template_file Path to saved template file
#' @return Spec object of class spec_hierarchical
#' @export
spec_hierarchical_template <- function(template = NULL, template_file = NULL) {
  if (is.null(template) && is.null(template_file)) {
    stop("Provide either a template object or template_file")
  }
  structure(list(template = template, template_file = template_file), class = "spec_hierarchical")
}

#' @exportS3Method
encode_spec.spec_hierarchical <- function(x, spec, mask, materialize, label, ...) {
  tmpl <- spec$template
  if (is.null(tmpl)) {
    if (is.null(spec$template_file)) stop("template_file is required when template is NULL")
    tmpl <- load_hierarchical_template(spec$template_file)
  }
  encode_hierarchical(x, tmpl, label = label)
}

#' @exportS3Method
encode_spec.spec_time_slepian <- function(x, spec, mask, materialize, label, ...) {
  n_time <- nrow(x)
  k_use <- spec$k %||% max(1L, floor(2 * n_time * spec$bandwidth * spec$tr) - 1L)
  if (materialize == "matrix") {
    basis <- dpss_time_basis(n_time, tr = spec$tr, bandwidth = spec$bandwidth, k = k_use, backend = spec$backend)
  } else {
    basis <- slepian_temporal_handle(n_time = n_time, tr = spec$tr, bandwidth = spec$bandwidth, k = k_use, backend = spec$backend)
  }
  loadings <- Matrix::Matrix(crossprod(x, basis), sparse = FALSE)
  spc <- neuroim2::NeuroSpace(c(dim(as.array(mask)), n_time))
  meta <- list(
    family = "time_slepian",
    k = k_use,
    tr = spec$tr,
    bandwidth = spec$bandwidth,
    backend = spec$backend
  )
  LatentNeuroVec(basis = basis, loadings = loadings, space = spc, mask = mask,
                 offset = numeric(0), label = label, meta = meta)
}

#' @exportS3Method
encode_spec.spec_time_dct <- function(x, spec, mask, materialize, label, ...) {
  n_time <- nrow(x)
  k_use <- spec$k
  if (materialize == "matrix") {
    basis <- build_dct_basis(n_time, k = k_use, norm = spec$norm)
  } else {
    basis <- dct_basis_handle(n_time = n_time, k = k_use, norm = spec$norm)
  }
  loadings <- Matrix::Matrix(crossprod(x, basis), sparse = FALSE)
  spc <- neuroim2::NeuroSpace(c(dim(as.array(mask)), n_time))
  meta <- list(
    family = "time_dct",
    k = k_use,
    norm = spec$norm
  )
  LatentNeuroVec(basis = basis, loadings = loadings, space = spc, mask = mask,
                 offset = numeric(0), label = label, meta = meta)
}

#' @exportS3Method
encode_spec.spec_time_bspline <- function(x, spec, mask, materialize, label, ...) {
  n_time <- nrow(x)
  k_use <- spec$k
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
      id = NULL,
      label = NULL
    )
    basis@spec$orthonormalize <- spec$orthonormalize
  }
  loadings <- Matrix::Matrix(crossprod(x, basis_mat(basis)), sparse = FALSE)
  spc <- neuroim2::NeuroSpace(c(dim(as.array(mask)), n_time))
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

#' @exportS3Method
encode_spec.spec_space_slepian <- function(x, spec, mask, reduction, materialize, label, ...) {
  mask_arr <- .mask_to_array(mask, "encode_spec.spec_space_slepian")
  if (is.null(reduction)) reduction <- make_cluster_reduction(mask, seq_len(sum(mask_arr)))
  if (materialize == "matrix") {
    loadings <- lift(reduction, basis_slepian(k = spec$k), k_neighbors = spec$k_neighbors)
  } else {
    loadings <- slepian_spatial_loadings_handle(reduction, basis_spec = basis_slepian(k = spec$k), data = NULL, label = "slepian-spatial")
  }
  basis <- as.matrix(x) %*% as.matrix(loadings_mat(loadings))
  spc <- neuroim2::NeuroSpace(c(dim(mask_arr), nrow(x)))
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
  mask_arr <- .mask_to_array(mask, "encode_spec.spec_space_heat")
  if (is.null(reduction)) reduction <- make_cluster_reduction(mask, seq_len(sum(mask_arr)))
  spec_hw <- basis_heat_wavelet(scales = spec$scales, order = spec$order, threshold = spec$threshold)
  if (materialize == "matrix") {
    loadings <- lift(reduction, spec_hw, k_neighbors = spec$k_neighbors)
  } else {
    loadings <- heat_wavelet_loadings_handle(reduction, basis_spec = spec_hw, data = NULL, label = "heat-wavelet")
  }
  basis <- as.matrix(x) %*% as.matrix(loadings_mat(loadings))
  spc <- neuroim2::NeuroSpace(c(dim(mask_arr), nrow(x)))
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
  mask_arr <- .mask_to_array(mask, "encode_spec.spec_space_hrbf")
  params <- spec$params %||% list()
  B_atoms <- hrbf_generate_basis(params, mask) # atoms x vox
  loadings <- Matrix::t(Matrix::Matrix(B_atoms, sparse = TRUE)) # vox x atoms
  gram <- as.matrix(B_atoms %*% Matrix::t(B_atoms))
  rhs <- t(as.matrix(x) %*% Matrix::t(B_atoms))
  coeff <- t(.robust_gram_solve(gram, rhs))
  basis <- Matrix::Matrix(coeff, sparse = FALSE) # time x atoms
  spc <- neuroim2::NeuroSpace(c(dim(mask_arr), nrow(x)))
  meta <- list(
    family = "space_hrbf",
    params = params
  )
  LatentNeuroVec(basis = basis, loadings = loadings, space = spc, mask = mask,
                 offset = numeric(0), label = label, meta = meta)
}

#' @exportS3Method
encode_spec.spec_space_pca <- function(x, spec, mask, reduction, materialize, label, ...) {
  mask_arr <- .mask_to_array(mask, "encode_spec.spec_space_pca")
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
      stop("PCA whitening requested, but singular values were not returned by lift().", call. = FALSE)
    }
    if (any(!is.finite(d)) || any(d <= 0)) {
      stop("PCA whitening requires strictly positive finite singular values.", call. = FALSE)
    }
    scale_fac <- sqrt(max(1, n_time - 1))
    basis <- sweep(as.matrix(basis), 2, d, "/") * scale_fac
    loadings <- loadings %*% Matrix::Diagonal(x = d / scale_fac)
  }

  spc <- neuroim2::NeuroSpace(c(dim(mask_arr), n_time))
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
  mask_arr <- .mask_to_array(mask, "encode_spec.spec_space_wavelet_active")
  wavelet_active_latent(
    X = x,
    mask = mask,
    levels_space = spec$levels_space,
    levels_time = spec$levels_time,
    threshold = spec$threshold
  )
}

#' @exportS3Method
encode_spec.spec_st <- function(x, spec, mask, reduction, materialize, label, ...) {
  mask_arr <- .mask_to_array(mask, "encode_spec.spec_st")
  n_time <- nrow(x)

  if (inherits(spec$time, "spec_time_slepian")) {
    k_time <- spec$time$k %||% max(1L, floor(2 * n_time * spec$time$bandwidth * spec$time$tr) - 1L)
    B_t <- basis_mat(slepian_temporal_handle(n_time = n_time, tr = spec$time$tr,
                                             bandwidth = spec$time$bandwidth, k = k_time,
                                             backend = spec$time$backend))
  } else if (inherits(spec$time, "spec_time_bspline")) {
    k_time <- spec$time$k
    B_t <- as.matrix(build_bspline_basis(
      n_time = n_time,
      k = k_time,
      degree = spec$time$degree,
      include_intercept = spec$time$include_intercept,
      orthonormalize = spec$time$orthonormalize
    ))
    if (!spec$time$orthonormalize) {
      q <- qr(B_t)
      B_t <- qr.Q(q)
    }
  } else {
    stop("Unsupported spec_st$time class: ", class(spec$time)[1])
  }

  if (is.null(reduction)) reduction <- make_cluster_reduction(mask, seq_len(sum(mask_arr)))
  if (inherits(spec$space, "spec_space_slepian")) {
    L_s <- as.matrix(lift(reduction, basis_slepian(k = spec$space$k), k_neighbors = spec$space$k_neighbors))
    B_atoms <- t(L_s)
  } else if (inherits(spec$space, "spec_space_hrbf")) {
    params <- spec$space$params %||% list()
    B_atoms <- as.matrix(hrbf_generate_basis(params, mask))
    L_s <- t(B_atoms)
  } else if (inherits(spec$space, "spec_space_wavelet_active")) {
    stop("spec_st with space = spec_space_wavelet_active not yet supported; use spec_space_wavelet_active directly via encode().", call. = FALSE)
  } else {
    stop("Unsupported spec_st$space class: ", class(spec$space)[1])
  }

  gram <- B_atoms %*% t(B_atoms)
  rhs_st <- t(as.matrix(x) %*% t(B_atoms))
  C_atoms <- t(.robust_gram_solve(as.matrix(gram), rhs_st))

  core <- crossprod(B_t, C_atoms)

  decoder <- function(time_idx = NULL, roi_mask = NULL, ...) {
    t_sel <- if (is.null(time_idx)) seq_len(n_time) else as.integer(time_idx)
    B_sel <- B_t[t_sel, , drop = FALSE]
    rec <- B_sel %*% core %*% t(L_s)
    if (!is.null(roi_mask)) {
      global_idx <- which(as.logical(mask_arr))
      roi_global <- which(as.logical(roi_mask))
      col_keep <- which(global_idx %in% roi_global)
      rec <- rec[, col_keep, drop = FALSE]
    }
    rec
  }

  meta <- list(
    family = "st_separable",
    time = class(spec$time)[1],
    space = class(spec$space)[1],
    k_time = ncol(B_t),
    k_space = ncol(L_s),
    label = label
  )

  implicit_latent(
    coeff = list(core = core, B_t = B_t, L_s = L_s),
    decoder = decoder,
    meta = meta,
    mask = mask_arr
  )
}
