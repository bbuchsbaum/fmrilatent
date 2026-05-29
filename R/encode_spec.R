# Spec constructors for the encode() pipeline.
#
# These S3 list-with-class objects are consumed by encode_spec() S3 methods
# defined in R/encode.R (and per-family in R/awpt.R, R/parcel_basis.R, etc.).
# Validators (.validate_positive_scalar, .validate_positive_count,
# .validate_flag_scalar, .validate_nonnegative_scalar,
# .validate_nonnegative_count, .validate_hrbf_params) live in R/encode.R and
# are available because R/encode.R is loaded before this file via Collate.

#' @include encode.R
NULL

#' Temporal Slepian/DPSS spec
#'
#' @param tr Repetition time (seconds).
#' @param bandwidth Half-bandwidth in Hz (default 0.1).
#' @param k Optional number of components (default floor(2*NW)-1).
#' @param backend Backend to use. Only "tridiag" is currently supported.
#' @return A `spec_time_slepian` object for `encode()` / `spec_st()`.
#' @export
spec_time_slepian <- function(tr, bandwidth = 0.1, k = NULL, backend = c("tridiag", "dense")) {
  backend <- match.arg(backend)
  if (identical(backend, "dense")) {
    .encoder_cli_abort(
      "backend = \"dense\" is disabled because it can return a different DPSS subspace under eigenvalue degeneracy; use backend = \"tridiag\".",
      class = "fmrilatent_error_unsupported_dpss_backend"
    )
  }
  tr <- .validate_positive_scalar(tr, "tr")
  bandwidth <- .validate_positive_scalar(bandwidth, "bandwidth")
  if (!is.null(k)) {
    k <- .validate_positive_count(k, "k")
  }
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
  k <- .validate_positive_count(k, "k")
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
  k <- .validate_positive_count(k, "k")
  degree <- .validate_positive_count(degree, "degree")
  include_intercept <- .validate_flag_scalar(include_intercept, "include_intercept")
  orthonormalize <- .validate_flag_scalar(orthonormalize, "orthonormalize")
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
  k <- .validate_positive_count(k, "k")
  k_neighbors <- .validate_positive_count(k_neighbors, "k_neighbors")
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
  k <- .validate_positive_count(k, "k")
  center <- .validate_flag_scalar(center, "center")
  whiten <- .validate_flag_scalar(whiten, "whiten")
  structure(
    list(
      k = k,
      center = center,
      whiten = whiten,
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
  if (!is.numeric(scales) || length(scales) == 0L || any(!is.finite(scales)) || any(scales <= 0)) {
    stop("scales must be a non-empty numeric vector of positive finite values.", call. = FALSE)
  }
  order <- .validate_positive_count(order, "order")
  threshold <- .validate_nonnegative_scalar(threshold, "threshold")
  k_neighbors <- .validate_positive_count(k_neighbors, "k_neighbors")
  structure(list(scales = as.numeric(scales), order = order, threshold = threshold, k_neighbors = k_neighbors),
            class = "spec_space_heat")
}

#' Spatial HRBF spec
#'
#' @param params HRBF parameter list (sigma0, levels, radius_factor, kernel_type, seed).
#' @return A `spec_space_hrbf` object.
#' @export
spec_space_hrbf <- function(params = list()) {
  structure(list(params = .validate_hrbf_params(params)), class = "spec_space_hrbf")
}

#' Spatial wavelet (active pencil) spec
#'
#' @param levels_space Spatial lifting levels.
#' @param levels_time Optional time lifting levels.
#' @param threshold Threshold after space transform.
#' @return A `spec_space_wavelet_active` object.
#' @export
spec_space_wavelet_active <- function(levels_space = 2L, levels_time = 0L, threshold = 0) {
  levels_space <- .validate_positive_count(levels_space, "levels_space")
  levels_time <- .validate_nonnegative_count(levels_time, "levels_time")
  threshold <- .validate_nonnegative_scalar(threshold, "threshold")
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
