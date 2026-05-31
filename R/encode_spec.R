# Spec constructors for the encode() pipeline.
#
# These S3 list-with-class objects are consumed by encode_spec() S3 methods
# defined in R/encode.R (and per-family in R/awpt.R, R/parcel_basis.R, etc.).
# Validators (.validate_positive_scalar, .validate_positive_count,
# .validate_optional_positive_count, .validate_flag_scalar,
# .validate_nonnegative_scalar, .validate_nonnegative_count,
# .validate_hrbf_params) live in R/encoder_validators.R;
# .encoder_cli_abort lives in R/encoder_utilities.R.

#' @include encode.R encoder_utilities.R encoder_validators.R
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
#' @param k Optional number of components. If `NULL`, the encoder uses all
#'   available DCT columns (`k = n_time`) at encode time.
#' @param norm Normalization ("ortho" or "none").
#' @return A `spec_time_dct` object.
#' @export
spec_time_dct <- function(k = NULL, norm = c("ortho", "none")) {
  norm <- match.arg(norm)
  k <- .validate_optional_positive_count(k, "k")
  structure(list(k = k, norm = norm), class = "spec_time_dct")
}

#' Temporal B-spline spec
#'
#' @param k Optional number of components (df). If `NULL`, the encoder uses
#'   `min(5, n_time)` at encode time.
#' @param degree Spline degree (default 3).
#' @param include_intercept Logical include intercept.
#' @param orthonormalize Logical orthonormalize columns (default TRUE).
#' @return A `spec_time_bspline` object.
#' @export
spec_time_bspline <- function(k = NULL, degree = 3L, include_intercept = FALSE, orthonormalize = TRUE) {
  k <- .validate_optional_positive_count(k, "k")
  degree <- .validate_positive_count(degree, "degree")
  include_intercept <- .validate_flag_scalar(include_intercept, "include_intercept")
  orthonormalize <- .validate_flag_scalar(orthonormalize, "orthonormalize")
  structure(list(k = k, degree = degree, include_intercept = include_intercept,
                 orthonormalize = orthonormalize),
            class = "spec_time_bspline")
}

#' Spatial Slepian spec
#'
#' @param k Optional components per cluster. If `NULL`, encoders use their
#'   family default (`3L`).
#' @param k_neighbors k-NN graph parameter.
#' @return A `spec_space_slepian` object.
#' @export
spec_space_slepian <- function(k = NULL, k_neighbors = 6L) {
  k <- .validate_optional_positive_count(k, "k")
  k_neighbors <- .validate_positive_count(k_neighbors, "k_neighbors")
  structure(list(k = k, k_neighbors = k_neighbors), class = "spec_space_slepian")
}

#' Spatial PCA spec (cluster-local)
#'
#' Computes PCA eigenvectors within each cluster/parcel specified by a
#' `ClusterReduction` and returns a block-sparse spatial dictionary.
#'
#' @param k Optional components per cluster. If `NULL`, encoders use their
#'   family default (`3L`).
#' @param center Logical; center voxels before PCA (default TRUE). When TRUE,
#'   voxel means are stored in `LatentNeuroVec@offset`.
#' @param scale Logical; scale voxels before PCA (default FALSE). Passed
#'   through to `lift(ClusterReduction, spec_pca)`.
#' @param whiten Logical; if TRUE, return whitened scores (unit-variance) and
#'   rescaled loadings such that reconstruction is unchanged.
#' @param backend SVD backend: "auto" (default), "svds" (RSpectra), or "svd" (base).
#' @return A `spec_space_pca` object.
#' @export
spec_space_pca <- function(k = NULL, center = TRUE, whiten = FALSE,
                           backend = c("auto", "svds", "svd"), scale = FALSE) {
  backend <- match.arg(backend)
  k <- .validate_optional_positive_count(k, "k")
  center <- .validate_flag_scalar(center, "center")
  scale <- .validate_flag_scalar(scale, "scale")
  whiten <- .validate_flag_scalar(whiten, "whiten")
  structure(
    list(
      k = k,
      center = center,
      scale = scale,
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
#' @param sparsify_eps Non-negative threshold for small heat-wavelet
#'   coefficients. This is stored as `threshold` for compatibility.
#' @param threshold Deprecated alias for `sparsify_eps`.
#' @param k_neighbors k-NN graph parameter.
#' @return A `spec_space_heat` object.
#' @export
spec_space_heat <- function(scales = c(1, 2, 4, 8), order = 30L,
                            threshold = NULL, k_neighbors = 6L,
                            sparsify_eps = NULL) {
  if (!is.numeric(scales) || length(scales) == 0L || any(!is.finite(scales)) || any(scales <= 0)) {
    .encoder_cli_abort(
      "scales must be a non-empty numeric vector of positive finite values.",
      class = "fmrilatent_error_invalid_scalar"
    )
  }
  order <- .validate_positive_count(order, "order")
  if (!is.null(sparsify_eps) && !is.null(threshold) && !identical(sparsify_eps, threshold)) {
    .encoder_cli_abort(
      "Specify only one of sparsify_eps or threshold.",
      class = "fmrilatent_error_invalid_argument"
    )
  }
  sparsify_eps <- sparsify_eps %||% threshold %||% 1e-6
  sparsify_eps <- .validate_nonnegative_scalar(sparsify_eps, "sparsify_eps")
  k_neighbors <- .validate_positive_count(k_neighbors, "k_neighbors")
  structure(
    list(
      scales = as.numeric(scales),
      order = order,
      sparsify_eps = sparsify_eps,
      threshold = sparsify_eps,
      k_neighbors = k_neighbors
    ),
    class = "spec_space_heat"
  )
}

#' Spatial HRBF spec
#'
#' @param sigma0 Base kernel width. `NULL` uses the HRBF default (`6`).
#' @param levels Number of dyadic resolution levels. `NULL` uses the HRBF
#'   default (`3L`).
#' @param radius_factor Atom support radius as a multiple of level width.
#'   `NULL` uses the HRBF default (`2.5`).
#' @param num_extra_fine_levels Optional extra fine levels. `NULL` uses the
#'   HRBF default (`0L`).
#' @param kernel_type Radial kernel family. `NULL` uses `"gaussian"`.
#' @param kernel_type_fine_levels Kernel family for extra fine levels. `NULL`
#'   uses the HRBF default.
#' @param seed Random seed. `NULL` uses the HRBF default (`1L`).
#' @param params Optional named list for advanced or compatibility HRBF
#'   parameters. Explicit formals override entries with the same names.
#' @param ... Additional named HRBF parameters for advanced use.
#' @return A `spec_space_hrbf` object.
#' @export
spec_space_hrbf <- function(params = list(), sigma0 = NULL, levels = NULL,
                            radius_factor = NULL, num_extra_fine_levels = NULL,
                            kernel_type = NULL, kernel_type_fine_levels = NULL,
                            seed = NULL, ...) {
  if (!is.list(params)) {
    .encoder_cli_abort("params must be a list.", class = "fmrilatent_error_invalid_params")
  }
  extra <- list(...)
  if (length(extra) > 0L &&
      (is.null(names(extra)) || any(!nzchar(names(extra))))) {
    .encoder_cli_abort(
      "... must contain only named HRBF parameters.",
      class = "fmrilatent_error_invalid_params"
    )
  }
  params <- c(params, extra)
  explicit <- list(
    sigma0 = sigma0,
    levels = levels,
    radius_factor = radius_factor,
    num_extra_fine_levels = num_extra_fine_levels,
    kernel_type = kernel_type,
    kernel_type_fine_levels = kernel_type_fine_levels,
    seed = seed
  )
  explicit <- explicit[!vapply(explicit, is.null, logical(1))]
  params[names(explicit)] <- explicit
  structure(list(params = .validate_hrbf_params(params)), class = "spec_space_hrbf")
}

#' Spatial wavelet (active pencil) spec
#'
#' @param levels_space Spatial lifting levels.
#' @param levels_time Optional time lifting levels.
#' @param select_threshold Non-negative active-set coefficient threshold after
#'   the spatial transform. This is stored as `threshold` for compatibility.
#' @param threshold Deprecated alias for `select_threshold`.
#' @return A `spec_space_wavelet_active` object.
#' @export
spec_space_wavelet_active <- function(levels_space = 2L, levels_time = 0L,
                                      select_threshold = NULL, threshold = NULL) {
  levels_space <- .validate_positive_count(levels_space, "levels_space")
  levels_time <- .validate_nonnegative_count(levels_time, "levels_time")
  if (!is.null(select_threshold) && !is.null(threshold) &&
      !identical(select_threshold, threshold)) {
    .encoder_cli_abort(
      "Specify only one of select_threshold or threshold.",
      class = "fmrilatent_error_invalid_argument"
    )
  }
  select_threshold <- select_threshold %||% threshold %||% 0
  select_threshold <- .validate_nonnegative_scalar(select_threshold, "select_threshold")
  structure(
    list(
      levels_space = levels_space,
      levels_time = levels_time,
      select_threshold = select_threshold,
      threshold = select_threshold
    ),
    class = "spec_space_wavelet_active"
  )
}

#' Spatiotemporal spec (separable)
#'
#' @param time Temporal spec (`spec_time_*`).
#' @param space Spatial spec (`spec_space_*`).
#' @param core_mode Controls the returned representation. `"auto"` (the
#'   default) always returns an `ImplicitLatent` (decoder-only separable
#'   form). `"explicit"` returns an explicit [LatentNeuroVec] (the separable
#'   `basis %*% t(loadings)` pair) when both the time and space factors are
#'   materializable, which holds for every currently supported `spec_st`
#'   combination; if the factors are not both materializable it falls back to
#'   `ImplicitLatent` with a `fmrilatent_warning_explicit_core_unavailable`
#'   warning. The reconstruction is identical either way.
#' @return A `spec_st` object.
#' @export
spec_st <- function(time, space, core_mode = c("auto", "explicit")) {
  core_mode <- match.arg(core_mode)
  structure(list(time = time, space = space, core_mode = core_mode), class = "spec_st")
}
