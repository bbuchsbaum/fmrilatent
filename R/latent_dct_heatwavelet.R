# Heat-wavelet template constructor (basis is a zero placeholder).

#' Create a template LatentNeuroVec with heat-wavelet spatial loadings
#'
#' Builds a \code{LatentNeuroVec} whose loadings are a heat-wavelet
#' \code{LoadingsHandle} and whose basis is a placeholder zero matrix of
#' \code{n_time x k} (where \code{k} is determined by the heat-wavelet
#' loadings, not by the caller). The caller is expected to overwrite
#' \code{lvec@basis} with real coefficients (e.g. fitted from data) before
#' the object represents a valid factorization.
#'
#' Despite the name, no DCT basis is constructed: the "dct" reference
#' predates the spec/encoder pipeline and is retained for API
#' compatibility. For an encoded DCT-temporal + heat-wavelet-spatial
#' pipeline see \code{spec_st(time = spec_time_dct(...), space = spec_space_heat(...))}
#' passed to \code{\link{encode}}.
#'
#' @param n_time        Number of time points.
#' @param k_time        Optional ignored legacy argument. The number of
#'   components is determined by the heat-wavelet loadings. Supplying a
#'   non-NULL value warns with class `fmrilatent_warning_deprecated`.
#' @param mask          LogicalNeuroVol or logical array mask (3D).
#' @param cluster_map   Optional integer vector mapping voxels (mask order) to clusters.
#' @param reduction     Graph reduction object; if NULL, built via
#'   `make_cluster_reduction(mask, cluster_map)` with default one-cluster-per-voxel map.
#' @param hw_basis_spec Heat-wavelet basis spec; defaults to `basis_heat_wavelet()`.
#' @param offset        Optional voxel-wise offset (length n_vox).
#' @param label         Optional label.
#'
#' @return A \code{LatentNeuroVec} with placeholder basis matrix.
#' @export
latent_dct_heatwavelet <- function(n_time,
                                   k_time = NULL,
                                   mask,
                                   cluster_map   = NULL,
                                   reduction     = NULL,
                                   hw_basis_spec = NULL,
                                   offset        = numeric(0),
                                   label         = "DCT + heat-wavelet") {
  if (missing(n_time)) {
    .encoder_cli_abort("n_time must be a positive integer.",
                       class = "fmrilatent_error_invalid_count")
  }
  n_time <- .validate_positive_count(n_time, "n_time")
  if (!is.null(k_time)) {
    .encoder_cli_warn(
      "k_time is ignored by latent_dct_heatwavelet(); omit it in new code.",
      class = "fmrilatent_warning_deprecated",
      call = rlang::caller_env()
    )
  }

  # Build mask volume and array first

  mask_arr <- .mask_to_array(mask, "latent_dct_heatwavelet")
  mask_vol <- if (inherits(mask, "LogicalNeuroVol")) {
    mask
  } else {
    neuroim2::LogicalNeuroVol(mask_arr, neuroim2::NeuroSpace(dim(mask_arr)))
  }

  # Spatial reduction + heat-wavelet spec
  if (is.null(hw_basis_spec)) {
    hw_basis_spec <- basis_heat_wavelet()
  }
  if (is.null(reduction)) {
    if (is.null(cluster_map)) {
      cluster_map <- seq_len(sum(mask_arr))
    }
    reduction <- make_cluster_reduction(mask, cluster_map)
  }

  # Create loadings handle first - this determines the number of components

  L_handle <- heat_wavelet_loadings_handle(
    reduction  = reduction,
    basis_spec = hw_basis_spec,
    data       = NULL,
    label      = "heat-wavelet"
  )

  # Get dimensions from loadings
  n_vox <- .latent_loadings_dim(L_handle)[1L]
  k_components <- .latent_loadings_dim(L_handle)[2L]

  # Validate offset
  if (!length(offset)) {
    offset <- numeric(n_vox)
  } else if (length(offset) != n_vox) {
    .encoder_cli_abort(
      paste0("length(offset) (", length(offset),
             ") must equal n_vox (", n_vox, ")."),
      class = "fmrilatent_error_dim",
      call = rlang::caller_env()
    )
  }

  # Create explicit basis matrix with matching number of columns

  # Initialized to zeros - users should populate with actual coefficients
  basis0 <- Matrix::Matrix(
    0,
    nrow = as.integer(n_time),
    ncol = k_components,
    sparse = FALSE
  )

  # Build 4D space with any geometry carried by a LogicalNeuroVol mask.
  spc <- .space_with_time_from_mask(mask, n_time, "latent_dct_heatwavelet")

  LatentNeuroVec(
    basis    = basis0,
    loadings = L_handle,
    space    = spc,
    mask     = mask_vol,
    offset   = offset,
    label    = label,
    meta     = list(
      time_basis   = "template",
      time_k       = k_components,
      spatial_dict = "heat_wavelet"
    )
  )
}
